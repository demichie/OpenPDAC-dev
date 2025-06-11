# Salvalo come vtk_to_ascii_all_vars.py
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LightSource
from linecache import getline

import matplotlib.pyplot as plt

import shutil
import os
import re
import sys
import numpy as np
import pyvista as pv
import netCDF4 as nc
from scipy.interpolate import griddata
import json
import glob
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process VTK data into raster maps.")
    parser.add_argument("--lx", type=float, default=None, help="Lower-left x-coordinate of the grid (optional).")
    parser.add_argument("--ly", type=float, default=None, help="Lower-left y-coordinate of the grid (optional).")
    parser.add_argument("--res", type=float, default=None, help="Grid resolution (optional).")
    parser.add_argument("--nrows", type=int, default=None, help="Number of rows in the output grid (optional).")
    parser.add_argument("--ncols", type=int, default=None, help="Number of columns in the output grid (optional).")
    parser.add_argument("-np", type=int, default=1, help="Number of processes to use for parallel execution.")
    return parser.parse_args()

def validate_args(args):
    # Verifica solo i parametri relativi alla griglia, ignora --np
    provided = {
        "lx": args.lx is not None,
        "ly": args.ly is not None,
        "res": args.res is not None,
        "nrows": args.nrows is not None,
        "ncols": args.ncols is not None
    }

    total_provided = sum(provided.values())

    if total_provided == 0:
        return "default"  # caso 1: nessun parametro fornito

    if total_provided == 1 and provided["res"]:
        return "res_only"  # caso 2: solo --res fornito

    if all(provided.values()):
        return "full"  # caso 3: tutti i parametri forniti

    raise ValueError(
        "Invalid combination of grid arguments.\n"
        "Accepted usage:\n"
        "1. No grid arguments\n"
        "2. Only --res\n"
        "3. All of --lx, --ly, --res, --nrows, --ncols\n"
        "(--np is always optional)"
    )


# Function to extract floating point number from directory name
def extract_float(directory_name):
    try:
        return float(directory_name)
    except ValueError:
        return None


def find_and_create_links(root_folder):
    vtk_files = {}
    top_dir = os.getcwd()
    os.chdir(root_folder)
    current_dir = os.getcwd()

    # Sort subfolders
    directories = [d for d in os.listdir() if os.path.isdir(d)]
    directories.sort(key=extract_float)

    # Iterate over subfolders
    for subdir in directories:
        subdir_path = os.path.join(current_dir, subdir)
        if os.path.isdir(subdir_path):
            # Iterate over files in subfolder
            for file in os.listdir(subdir_path):
                if file.endswith('.vtk'):
                    base_name = os.path.splitext(file)[0]
                    vtk_files.setdefault(base_name, []).append(
                        (subdir, os.path.join(subdir_path, file)))

    # Create symbolic links
    for base_name, files in vtk_files.items():
        index = 0
        for subdir, file_path in files:
            new_name = f"{base_name}_{index:03d}.vtk"
            new_link_path = os.path.join(current_dir, new_name)
            try:
                os.symlink(file_path, new_link_path)
                print(f"Created symlink: {new_link_path} -> {file_path}")
            except Exception:
                shutil.copy(file_path, new_link_path)
                print(f"Created copy: {new_link_path} -> {file_path}")

            index += 1

    os.chdir(top_dir)


def extract_surface_names_from_foam_dict(file_path):
    """
    Extracts surface names defined within the 'surfaces (...);' block
    in an OpenFOAM dictionary file.

    Args:
        file_path (str): The path to the OpenFOAM dictionary file.

    Returns:
        list: A list of strings containing the extracted surface names,
              or an empty list if the block is not found or is empty,
              or None if the file cannot be read.
    """
    surface_names = []
    # State machine states:
    # 0: Searching for the 'surfaces' keyword at the start of a block
    # 1: Found 'surfaces', searching for the opening parenthesis '('
    # 2: Found '(', inside the surfaces block, searching for names or ')'
    state = 0
    potential_name = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                cleaned_line = line.strip()
                # Remove potential end-of-line comments starting with //
                if '//' in cleaned_line:
                    cleaned_line = cleaned_line.split('//', 1)[0].strip()

                # Ignore empty lines or lines that are purely comments
                if not cleaned_line or cleaned_line.startswith(
                        '/*') or cleaned_line.endswith('*/'):
                    continue

                # --- State Machine Logic ---
                if state == 0:
                    # Look for the start of the surfaces definition block
                    if cleaned_line == "surfaces":
                        state = 1
                        potential_name = None  # Reset just in case

                elif state == 1:
                    # Expecting the opening parenthesis
                    if cleaned_line == "(":
                        state = 2
                    elif cleaned_line:  # Don't reset on empty/comment lines
                        state = 0

                elif state == 2:
                    # Inside the surfaces block
                    if cleaned_line == ");":
                        state = 0
                        potential_name = None

                    elif cleaned_line == "{":

                        if potential_name:
                            surface_names.append(potential_name)
                            potential_name = None  # Reset after use
                    elif cleaned_line:

                        potential_name = cleaned_line.split()[0]

                        if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$',
                                        potential_name):
                            potential_name = None

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None
    except Exception as e:
        print(f"An error occurred while reading {file_path}: {e}")
        return None

    return surface_names


def read_dict(filepath):

    with open(filepath, 'r') as f:
        content = f.read()
    content = re.sub(r'//.*', '', content)
    content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
    write_interval = float(
        re.search(r'writeInterval\s+([\d.]+)', content).group(1))

    surface_names = extract_surface_names_from_foam_dict(filepath)
    fields_block = re.search(r'fields\s*\((.*?)\)', content, re.DOTALL)
    fields = re.findall(r'[\w\.]+',
                        fields_block.group(1)) if fields_block else []
    return write_interval, surface_names, fields


def read_topoGridDict(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    xVent = float(re.search(r'xVent\s+([\d\.\-eE]+)', content).group(1))
    yVent = float(re.search(r'yVent\s+([\d\.\-eE]+)', content).group(1))
    raster_match = re.search(r'rasterFile\s+(\S+);', content)
    DEMfile = raster_match.group(1) if raster_match else None
    return xVent, yVent, DEMfile


def vtk_to_grid(filename, xVent, yVent, grid_x, grid_y):

    mesh = pv.read(filename)
    points = mesh.points[:, :2] + np.array([xVent, yVent])
    data = {}

    for name in mesh.array_names:
        if 'alpha' in name:
            values = mesh[name]
            interp = griddata(points,
                              values, (grid_x, grid_y),
                              method='linear',
                              fill_value=0.0)
            data[name] = interp

    return data

def process_single_vtk(var_args, fixed_args):

    (xVent, yVent, grid_x, grid_y,output_root, extent_map, xi, yi, Zinit, extent, cell) = fixed_args
    (i, time_val, filename, surface) = var_args

    print('filename',filename)
    data = vtk_to_grid(filename, xVent, yVent, grid_x, grid_y)
    
    time_folder = os.path.join(output_root, f"{time_val:.6g}", surface)
    os.makedirs(time_folder, exist_ok=True)

    levels = np.linspace(-6.0, -0.5, 12)
    label_str = 'Log particles volume fraction'

    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    for var_name, grid in data.items():

        if 'alpha' in var_name:

            with np.errstate(divide='ignore', invalid='ignore'):
                val_plot = np.flipud(np.log10(grid))
                
            val_plot[val_plot < -6.0] = np.nan

            fig, ax = plt.subplots()

            ls = LightSource(azdeg=315, altdeg=45)

            ax.imshow(ls.hillshade(np.flipud(Zinit), vert_exag=1.0, dx=cell, dy=cell),
                      cmap='gray',
                      extent=extent)

            ax.set_aspect('equal', 'box')

            p1 = ax.imshow(val_plot,
                   cmap=cmap,
                   norm=norm,
                   interpolation='nearest',
                   extent=extent_map,
                   alpha=0.5)
            clb = plt.colorbar(p1)
            clb.set_label('Log10 of Volume fraction',
                  labelpad=5,
                  y=0.5,
                  rotation=90)

            png_path = os.path.join(time_folder, f"{var_name}.png")
            print('Saving',png_path)
            plt.savefig(png_path, dpi=200)
            plt.close(fig)

            asc_path = os.path.join(time_folder, f"{var_name}.asc")
            print('Saving',asc_path)
            save_ascii_grid(grid, xi, yi, asc_path)

    return i

def readASC(DEM_file):

    print('Reading DEM file: ', DEM_file)
    # Parse the topography header
    hdr = [getline(DEM_file, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    xs_DEM = lx + 0.5 * cell + np.linspace(0, (cols - 1) * cell, cols)
    ys_DEM = ly + 0.5 * cell + np.linspace(0, (rows - 1) * cell, rows)

    extent = lx, lx + cols * cell, ly, ly + rows * cell

    # Load the topography into a numpy array
    DEM = np.loadtxt(DEM_file, skiprows=6)
    DEM = np.flipud(DEM)
    DEM[DEM == nd] = 0.0

    xinit = xs_DEM
    yinit = ys_DEM

    xmin = np.amin(xinit)
    xmax = np.amax(xinit)

    print('xmin,xmax', xmin, xmax)

    ymin = np.amin(yinit)
    ymax = np.amax(yinit)

    print('ymin,ymax', ymin, ymax)

    Xinit, Yinit = np.meshgrid(xinit, yinit)
    Zinit = DEM

    return Xinit, Yinit, Zinit, cell, extent


def save_ascii_grid(data, xi, yi, filepath):
    print("Saving ", filepath)
    nodata = -9999
    with open(filepath, 'w') as f:
        f.write(f"ncols         {len(xi)}\n")
        f.write(f"nrows         {len(yi)}\n")
        f.write(f"xllcorner     {xi[0]:.6f}\n")
        f.write(f"yllcorner     {yi[0]:.6f}\n")
        f.write(f"cellsize      {xi[1] - xi[0]:.6f}\n")
        f.write(f"NODATA_value  {nodata}\n")
        for row in np.flipud(data):
            f.write(" ".join(f"{v:.6f}" if np.isfinite(v) else f"{nodata}"
                             for v in row) + "\n")

def main(args,mode):

    write_interval, surfaces, _ = read_dict("system/sampleSurfaceAlpha")
    xVent, yVent, DEMfile = read_topoGridDict("system/topoGridDict")
    Xinit, Yinit, Zinit, cell, extent = readASC(DEMfile)

    terrain_folder = "postProcessing/terrain/*"
    terrain_files = sorted(
        glob.glob(os.path.join(terrain_folder, "terrain.vtk")))

    timesteps = [float(file.split('/')[2]) for file in terrain_files]
    timevals = np.array(timesteps)
    vtk_files = [val for _, val in sorted(zip(timevals, terrain_files))]
    times = np.sort(timevals)

    mesh = pv.read(vtk_files[0])
    points = mesh.points[:, :2] + np.array([xVent, yVent])
        
    nalpha = 0
    alphalist = []
    for name in mesh.array_names:
        if 'alpha' in name:
            nalpha +=1
            alphalist.append(name)

    if mode == "full":
   
        print("Running with parameters:")
        print(f"lx={args.lx}, ly={args.ly}, res={args.res}, nrows={args.nrows}, ncols={args.ncols}, np={args.np}")

        lxOut = args.lx
        lyOut = args.ly
        resolutionOut = args.res
        nrowsOut = args.nrows
        ncolsOut = args.ncols
        nproc = args.np

        xmin = lxOut + 0.5*resolutionOut
        xmax = lxOut + (ncolsOut-0.5)*resolutionOut
        ymin = lyOut + 0.5*resolutionOut
        ymax = lyOut + (nrowsOut-0.5)*resolutionOut

        xi = np.linspace(xmin, xmax, ncolsOut)
        yi = np.linspace(ymin, ymax, nrowsOut)
        grid_x, grid_y = np.meshgrid(xi, yi)
        extent_map = xi[0]-0.5*resolutionOut, xi[-1]+0.5*resolutionOut, \
                yi[0]-0.5*resolutionOut, yi[-1]+0.5*resolutionOut

    else:

        if mode == "default":
            print("Using default grid parameters.")
            resolutionOut = 20
        elif mode == "res_only":
            resolutionOut = args.res


        xmin, ymin = points.min(axis=0)
        xmax, ymax = points.max(axis=0)
        xi = np.arange(xmin, xmax + resolutionOut, resolutionOut)
        yi = np.arange(ymin, ymax + resolutionOut, resolutionOut)
        ncolsOut = len(xi)
        nrowsOut = len(yi)
        grid_x, grid_y = np.meshgrid(xi, yi)
        extent_map = xi[0]-0.5*resolutionOut, xi[-1]+0.5*resolutionOut, \
            yi[0]-0.5*resolutionOut, yi[-1]+0.5*resolutionOut


    alphaMax = np.zeros((nrowsOut,ncolsOut))

    levels = np.linspace(-6.0, -0.5, 12)
    label_str = 'Log particles volume fraction'

    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    output_root = "postProcessing/raster_maps/flow/"
    os.makedirs(output_root, exist_ok=True)

    metadata = {}
    counter = 0

    print('surfaces', surfaces)

    fixed_args = (xVent, yVent, grid_x, grid_y,output_root, extent_map, xi, yi, Zinit, extent, cell)

    for surface in surfaces:

        """
        nc_path = os.path.join(output_root, f"{surface}.nc")
        ncfile = nc.Dataset(nc_path, "w", format="NETCDF4")
        ncfile.createDimension("time", None)

        metadata[surface] = {"times": times, "fields": []}
        data_collections = {}
        """

        var_args_list = [(i, times[i], filename,surface) for i, filename in enumerate(vtk_files)]

        with ProcessPoolExecutor(max_workers=args.np) as executor:
            results = list(executor.map(partial(process_single_vtk, fixed_args=fixed_args), var_args_list))


        for i, var_name in enumerate(alphalist):
                
            for j, filename in enumerate(vtk_files):
        
                time_val = times[j]
                time_folder = os.path.join(output_root, f"{time_val:.6g}", surface)
                asc_path = os.path.join(time_folder, f"{var_name}.asc")
                print('Reading',asc_path)
                _, _, alpha_val, _, _ = readASC(asc_path)
                alphaMax = np.maximum(alphaMax,alpha_val)     
 
            with np.errstate(divide='ignore', invalid='ignore'):
                val_plot = np.flipud(np.log10(alphaMax))
            val_plot[val_plot < -6.0] = np.nan

            fig, ax = plt.subplots()

            ls = LightSource(azdeg=315, altdeg=45)

            ax.imshow(ls.hillshade(np.flipud(Zinit), vert_exag=1.0, dx=cell, dy=cell),
                      cmap='gray',
                      extent=extent)

            ax.set_aspect('equal', 'box')

            p1 = ax.imshow(val_plot,
                   cmap=cmap,
                   norm=norm,
                   interpolation='nearest',
                   extent=extent_map,
                   alpha=0.5)
            clb = plt.colorbar(p1)
            clb.set_label('Log10 of Volume fraction',
                  labelpad=5,
                  y=0.5,
                  rotation=90)

            max_folder = os.path.join(output_root, surface)

            os.makedirs(max_folder, exist_ok=True)

            png_path = os.path.join(max_folder, f"{var_name}_max.png")
            print('Saving',png_path)
            plt.savefig(png_path, dpi=200)
            plt.close(fig)

            asc_path = os.path.join(max_folder, f"{var_name}_max.asc")
            print('Saving',asc_path)
            save_ascii_grid(alphaMax, xi, yi, asc_path)

            plt.show()        
        
        """   
        ncfile.createDimension("y", len(yi))
        ncfile.createDimension("x", len(xi))
        ncfile.createVariable("time", "f4", ("time", ))[:] = times
        ncfile.createVariable("x", "f4", ("x", ))[:] = xi
        ncfile.createVariable("y", "f4", ("y", ))[:] = yi

        # alphaMax = np.max(...

        for ialpha in range(nalpha):

            var_name = alphalist[ialpha]

            for i, filename in enumerate(vtk_files):

                print(filename)
                time_val = times[i]

                time_folder = os.path.join(output_root, f"{time_val:.6g}", surface)
                os.makedirs(time_folder, exist_ok=True)

                grid = gridAll[:,:,ialpha,i]
                print(var_name)

                val_plot = np.flipud(np.log10(grid))
                val_plot[val_plot < -6.0] = np.nan

                if counter == 0:

                    p1 = ax.imshow(val_plot,
                                   cmap=cmap,
                                   norm=norm,
                                   interpolation='nearest',
                                   extent=extent_map,
                                   alpha=0.5)
                    clb = plt.colorbar(p1)

                    label_str = 'Log10 of Volume fraction'
                    clb.set_label(label_str,
                                  labelpad=5,
                                  y=0.5,
                                  rotation=90)

                else:

                    p1.set_data(val_plot)

                png_path = os.path.join(time_folder, f"{var_name}.png")
                print(png_path)

                plt.savefig(png_path, dpi=200)
                counter += 1

                asc_path = os.path.join(time_folder, f"{var_name}.asc")
                save_ascii_grid(grid, xi, yi, asc_path)

        """

        # ncfile.close()

    #with open(os.path.join(output_root, "metadata.json"), "w") as f:
    #     json.dump(metadata, f, indent=2)


if __name__ == "__main__":

    args = parse_arguments()
    mode = validate_args(args)
    print(f"Mode selected: {mode}")
    main(args, mode)
