import os
import re
import glob
import csv
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from linecache import getline
import pandas as pd
from matplotlib.colors import LightSource
from matplotlib.colors import BoundaryNorm
import argparse # Added for command-line arguments
import sys       # Added for sys.exit with argparse errors
from matplotlib.ticker import FuncFormatter

# from sklearn.preprocessing import StandardScaler
# from sklearn.cluster import DBSCAN

# Global default values that can be overridden by CLI args
DEFAULT_RASTER_RESOLUTION = 50.0 # Was 'step_dens' globally
PARTICLE_IMPACT_VELOCITY_THRESHOLD = 1.0 # Was 'toll' globally


def log_to_percent_formatter(x, pos):
    """
    Converte un valore log10(percentuale) in una stringa di etichetta percentuale.
    x: il valore del tick (log10 della percentuale).
    pos: la posizione del tick (non usata qui).
    """
    if np.isneginf(x): # Gestisce log10(0)
        return "0%"
    
    percentage = 10**x
    
    if percentage == 0: # Potrebbe accadere per underflow vicino a zero
        return "0%"
    # Formattazione per leggibilità
    elif percentage < 0.001: # Percentuali molto piccole
        return f"{percentage:.1e}%" # Es. 1.2e-4%
    elif percentage < 0.1: # Percentuali piccole
        return f"{percentage:.3f}%" # Es. 0.012%
    elif percentage < 1:    # Percentuali < 1%
        return f"{percentage:.2f}%" # Es. 0.12%
    elif percentage < 10:   # Percentuali tra 1% e 10%
        return f"{percentage:.1f}%" # Es. 5.3%
    elif percentage <= 100: # Percentuali fino al 100%
        return f"{percentage:.0f}%" # Es. 75%
    else: # Non dovrebbe accadere se la percentuale massima è 100
        return f"{percentage:.0f}%"


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def read_csv_to_arrays(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        columns = {field: [] for field in reader.fieldnames}
        for row in reader:
            for key in reader.fieldnames:
                columns[key].append(float(row[key]))
    arrays = {key.strip(): np.array(values) for key, values in columns.items()}
    return arrays


def readASC(DEM_file_path): # Renamed variable for clarity
    print('Reading DEM file: ' + DEM_file_path)
    hdr = [getline(DEM_file_path, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell_size_dem, nd_val_dem = values # Renamed for clarity
    cols = int(cols)
    rows = int(rows)

    xs_DEM = lx + 0.5 * cell_size_dem + np.linspace(0, (cols - 1) * cell_size_dem, cols)
    ys_DEM = ly + 0.5 * cell_size_dem + np.linspace(0, (rows - 1) * cell_size_dem, rows)
    dem_extent = lx, lx + cols * cell_size_dem, ly, ly + rows * cell_size_dem

    DEM_data = pd.read_table(DEM_file_path, sep='\s+', header=None, skiprows=6).astype(float).values
    DEM_data = np.flipud(DEM_data)
    DEM_data[DEM_data == nd_val_dem] = 0.0 # Assuming 0 is a safe fill value for NODATA in DEM for hillshade

    x_coords_dem_centers = xs_DEM
    y_coords_dem_centers = ys_DEM
    xmin_dem = np.amin(x_coords_dem_centers)
    xmax_dem = np.amax(x_coords_dem_centers)
    print('xmin_dem_center, xmax_dem_center (DEM cell centers):', xmin_dem, xmax_dem)
    ymin_dem = np.amin(y_coords_dem_centers)
    ymax_dem = np.amax(y_coords_dem_centers)
    print('ymin_dem_center, ymax_dem_center (DEM cell centers):', ymin_dem, ymax_dem)

    X_dem_centers, Y_dem_centers = np.meshgrid(x_coords_dem_centers, y_coords_dem_centers)
    Z_dem_data = DEM_data
    return X_dem_centers, Y_dem_centers, Z_dem_data, cell_size_dem, dem_extent


def read_topoGridDict(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    xVent = float(re.search(r'xVent\s+([\d\.\-eE]+)', content).group(1))
    yVent = float(re.search(r'yVent\s+([\d\.\-eE]+)', content).group(1))
    raster_match = re.search(r'rasterFile\s+(\S+);', content)
    DEMfile_path_from_dict = raster_match.group(1) if raster_match else None # Renamed variable
    return xVent, yVent, DEMfile_path_from_dict

# --- Command Line Argument Parsing Functions ---
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process OpenFOAM Lagrangian particle data to create ballistic impact maps.")
    parser.add_argument("--lx", type=float, default=None, help="Lower-left x-coordinate of the output raster grid (optional).")
    parser.add_argument("--ly", type=float, default=None, help="Lower-left y-coordinate of the output raster grid (optional).")
    parser.add_argument("--res", type=float, default=None, help="Output raster grid resolution (cell size) (optional).")
    parser.add_argument("--nrows", type=int, default=None, help="Number of rows in the output raster grid (optional).")
    parser.add_argument("--ncols", type=int, default=None, help="Number of columns in the output raster grid (optional).")
    # --np argument removed as per request
    return parser.parse_args()

def validate_args(args):
    provided = {
        "lx": args.lx is not None,
        "ly": args.ly is not None,
        "res": args.res is not None,
        "nrows": args.nrows is not None,
        "ncols": args.ncols is not None
    }
    total_provided = sum(provided.values())

    if total_provided == 0:
        return "default"
    if total_provided == 1 and provided["res"]:
        return "res_only"
    if all(provided.values()):
        return "full"
    
    # If none of the above, it's an invalid combination
    error_message = (
        "Invalid combination of grid arguments (--lx, --ly, --res, --nrows, --ncols).\n"
        "Accepted usage patterns for grid definition:\n"
        "1. No grid arguments: Grid extents derived from DEM, resolution uses script default.\n"
        "2. Only --res: Grid extents derived from DEM, resolution from --res.\n"
        "3. All of --lx, --ly, --res, --nrows, --ncols: Fully user-defined grid.\n"
    )
    # Using sys.exit(error_message) or parser.error(error_message) might be cleaner
    # but for now, a raised ValueError is fine.
    # To use parser.error, parse_arguments would need the parser instance.
    # For now, let's print and raise.
    print("ERROR: " + error_message, file=sys.stderr)
    raise ValueError("Invalid command-line arguments for grid definition.")

# --- Main Function ---
def main():
    args = parse_arguments()
    grid_definition_mode = validate_args(args)
    print(f"Grid definition mode selected: {grid_definition_mode}")

    xVent, yVent, DEMfile_path = read_topoGridDict("system/topoGridDict")
    print('DEMfile_path from topoGridDict:', DEMfile_path)
    # X_dem_centers, Y_dem_centers are meshes of DEM cell centers
    # Z_dem_data is the elevation data
    # dem_cell_size is the resolution of the input DEM
    # dem_extent is (xmin_edge, xmax_edge, ymin_edge, ymax_edge) of the DEM
    X_dem_centers, Y_dem_centers, Z_dem_data, dem_cell_size, dem_extent = readASC(DEMfile_path)


    ls = LightSource(azdeg=45, altdeg=45)
    current_dir = os.getcwd()

    output_raster_root = os.path.join(current_dir, "postProcessing", "raster_maps", "ballistics")
    os.makedirs(output_raster_root, exist_ok=True)
    postprocessing_dir = os.path.join(current_dir, "postProcessing")
    os.makedirs(postprocessing_dir, exist_ok=True)

    search_pattern = os.path.join(current_dir, "postProcessing", "cloudInfo1", "*", "output.csv")
    cloud_files = glob.glob(search_pattern)

    if not cloud_files:
        print(f"No particle CSV files found with pattern: {search_pattern}")
        return

    try:
        timesteps = [float(os.path.basename(os.path.dirname(file))) for file in cloud_files]
    except ValueError as e:
        print(f"Error extracting timesteps from file paths. Example file: {cloud_files[0]}")
        print(f"Ensure files are in directories like 'postProcessing/cloudInfo1/NUMERICAL_TIMESTEP/output.csv'. Error: {e}")
        return

    timevals = np.array(timesteps)
    sorted_indices = np.argsort(timevals)
    sorted_files = [cloud_files[i] for i in sorted_indices]
    sorted_times = timevals[sorted_indices]
    n_times = len(sorted_files)

    if n_times == 0:
        print("No particle CSV files found. Exiting.")
        return
    print(f"Found {n_times} timesteps.")

    print("Starting scan for unique particles...")
    all_particle_ids_set = set()
    particle_properties_map = {} 

    for filename in sorted_files:
        data_scan = read_csv_to_arrays(filename)
        current_origProc = data_scan['origProc'].astype(int)
        current_origId = data_scan['origId'].astype(int)
        current_d = data_scan['d']
        current_rho = data_scan['rho']

        for j in range(len(current_origId)):
            particle_key = (current_origProc[j], current_origId[j])
            all_particle_ids_set.add(particle_key)
            if particle_key not in particle_properties_map:
                particle_properties_map[particle_key] = {
                    'd': current_d[j],
                    'rho': current_rho[j]
                }
    
    unique_global_particle_identifiers = sorted(list(all_particle_ids_set))
    particle_to_global_idx_map = {
        pid: i for i, pid in enumerate(unique_global_particle_identifiers)
    }
    
    nballistics_global = len(unique_global_particle_identifiers)
    print(f"Total unique particles found (nballistics_global): {nballistics_global}")

    if nballistics_global == 0:
        print("No particles found in any file. Exiting.")
        return

    d_global = np.array([particle_properties_map[pid]['d'] for pid in unique_global_particle_identifiers])
    rho_global = np.array([particle_properties_map[pid]['rho'] for pid in unique_global_particle_identifiers])
    
    plot_titles = []
    unique_diams_global = np.unique(d_global)
    num_unique_sizes = len(unique_diams_global)

    if num_unique_sizes < 5:
        diams_for_plot_categories = unique_diams_global
        for i in range(len(diams_for_plot_categories) + 1):
            if i < len(diams_for_plot_categories):
                plot_titles.append('Diameter = ' + str(diams_for_plot_categories[i]) + 'm')
            else:
                plot_titles.append('All sizes')
    else:
        dMin_global = np.amin(d_global)
        dMax_global = np.amax(d_global)
        print('Global dMin:', dMin_global, 'Global dMax:', dMax_global)
        diams_for_plot_categories = np.linspace(dMin_global, dMax_global, 4) 
        
        plot_titles.append('Diameter <= ' + f"{diams_for_plot_categories[0]:.2e}" + 'm')
        for i in range(len(diams_for_plot_categories) - 1):
             plot_titles.append(
                    f"{diams_for_plot_categories[i]:.2e}" + 'm < Diameter <= ' +
                    f"{diams_for_plot_categories[i+1]:.2e}" + 'm')
        plot_titles.append('Diameter > ' + f"{diams_for_plot_categories[-1]:.2e}" + 'm')
        plot_titles.append('All sizes')

    print('Diameters for plot categorization (diams_for_plot_categories):', diams_for_plot_categories)

    A_velocity = np.full((nballistics_global, 3, n_times), np.nan)
    B_position = np.full((nballistics_global, 3, n_times), np.nan)
    matr_data = np.full((n_times, 8, nballistics_global), np.nan)

    print("Starting population of arrays A_velocity, B_position, matr_data...")
    for i_time, filename in enumerate(sorted_files):
        data = read_csv_to_arrays(filename)
        current_origProc = data['origProc'].astype(int)
        current_origId = data['origId'].astype(int)
        
        x_curr, y_curr, z_curr = data['x'], data['y'], data['z']
        Ux_curr, Uy_curr, Uz_curr = data['Ux'], data['Uy'], data['Uz']

        for local_particle_idx in range(len(current_origId)):
            proc, oid = current_origProc[local_particle_idx], current_origId[local_particle_idx]
            global_idx = particle_to_global_idx_map[(proc, oid)]

            A_velocity[global_idx, :, i_time] = [Ux_curr[local_particle_idx], Uy_curr[local_particle_idx], Uz_curr[local_particle_idx]]
            B_position[global_idx, :, i_time] = [x_curr[local_particle_idx]+xVent, y_curr[local_particle_idx]+yVent, z_curr[local_particle_idx]]
            
            matr_data[i_time, 1:4, global_idx] = B_position[global_idx, :, i_time]
            matr_data[i_time, 4:7, global_idx] = A_velocity[global_idx, :, i_time]
            if not np.isnan(A_velocity[global_idx, 0, i_time]):
                 matr_data[i_time, 7, global_idx] = LA.norm(A_velocity[global_idx, :, i_time])

    print("Calculating impact times...")
    t_impact_idx = np.zeros(nballistics_global, dtype=int)
    time_impact_val = np.zeros(nballistics_global, dtype=float)

    for s_global_idx in range(nballistics_global):
        for l_time_idx in range(3, n_times):
            current_vel_norm = matr_data[l_time_idx, 7, s_global_idx]
            prev_Uz = matr_data[l_time_idx - 1, 6, s_global_idx]

            if not np.isnan(current_vel_norm) and not np.isnan(prev_Uz):
                if (current_vel_norm < PARTICLE_IMPACT_VELOCITY_THRESHOLD and prev_Uz < 0):
                    t_impact_idx[s_global_idx] = l_time_idx
                    time_impact_val[s_global_idx] = sorted_times[l_time_idx]
                    break 

    print("Calculating mean and max velocities...")
    velocities_summary = np.full((nballistics_global, 4), np.nan)
    for k_global_idx in range(nballistics_global):
        velocities_summary[k_global_idx, 0] = k_global_idx 
        velocities_summary[k_global_idx, 1] = d_global[k_global_idx]
        impact_timestep_index = t_impact_idx[k_global_idx]
        if impact_timestep_index > 0: 
            particle_vel_norms_trajectory = matr_data[:impact_timestep_index + 1, 7, k_global_idx]
            valid_vels = particle_vel_norms_trajectory[~np.isnan(particle_vel_norms_trajectory)]
            if len(valid_vels) > 0:
                velocities_summary[k_global_idx, 2] = np.mean(valid_vels)
                velocities_summary[k_global_idx, 3] = np.amax(valid_vels)

    C_vel_headers = ['global_index', 'diameter [m]', 'mean_velocity [m/s]', 'max_velocity [m/s]']
    df_velocities = pd.DataFrame(velocities_summary)
    df_velocities[0] = df_velocities[0].astype(int)
    velocities_csv_path = os.path.join(postprocessing_dir, "velocities.csv")
    df_velocities.to_csv(velocities_csv_path, header=C_vel_headers, index=False, na_rep='NaN')
    print(f"Velocities summary saved to: {velocities_csv_path}")

    print("Calculating impact properties...")
    r_global = d_global / 2.0
    V_global = (4.0 / 3.0) * np.pi * (r_global**3)
    m_global = rho_global * V_global
    impact_properties_matrix = np.full((nballistics_global, 9), np.nan)

    for s_global_idx in range(nballistics_global):
        impact_properties_matrix[s_global_idx, 0] = s_global_idx
        current_impact_timestep_idx = t_impact_idx[s_global_idx]
        if current_impact_timestep_idx > 0: 
            impact_properties_matrix[s_global_idx, 1:4] = [
                d_global[s_global_idx], rho_global[s_global_idx], time_impact_val[s_global_idx]
            ]
            impact_properties_matrix[s_global_idx, 4:7] = matr_data[current_impact_timestep_idx, 1:4, s_global_idx] # x, y, z
            if current_impact_timestep_idx - 1 >= 0:
                 impact_velocity_norm = matr_data[current_impact_timestep_idx - 1, 7, s_global_idx]
                 impact_properties_matrix[s_global_idx, 7] = impact_velocity_norm
                 if not np.isnan(impact_velocity_norm) and not np.isnan(m_global[s_global_idx]):
                     impact_properties_matrix[s_global_idx, 8] = 0.5 * m_global[s_global_idx] * (impact_velocity_norm**2)

    C_impact_headers = ['global_index', 'diameter [m]', 'density [kg/m3]', 'impact_time [s]',
                        'x_impact [m]', 'y_impact [m]', 'z_impact [m]', 'impact_velocity_norm [m/s]',
                        'landing_energy [J]']
    df_impacts = pd.DataFrame(impact_properties_matrix)
    df_impacts[0] = df_impacts[0].astype(int)
    impacts_csv_path = os.path.join(postprocessing_dir, "impacts.csv")
    df_impacts.to_csv(impacts_csv_path, header=C_impact_headers, index=False, na_rep='NaN')
    print(f"Impact properties saved to: {impacts_csv_path}")

    filtered_impact_data = impact_properties_matrix[
        ~np.isnan(impact_properties_matrix[:, 3]) & (impact_properties_matrix[:, 3] > 0)
    ]
    
    if filtered_impact_data.shape[0] == 0:
        print("No particles impacted according to criteria. Cannot generate raster maps.")
    else:
        # ... (la logica per x_impact_coords, y_impact_coords, diam_impacted e la definizione della griglia rimane invariata) ...
        # ... (grid definition logic: map_x_ll, map_y_ll, current_map_res, etc.) ...
        # ... (particle binning logic: count_ballistic_class) ...
        
        print(f"Number of impacted particles for mapping: {filtered_impact_data.shape[0]}")
        x_impact_coords, y_impact_coords, diam_impacted = filtered_impact_data[:,4], filtered_impact_data[:,5], filtered_impact_data[:,1]

        if grid_definition_mode == "full":
            # ... (codice per la modalità "full") ...
            print("Using user-defined grid parameters for raster maps.")
            map_x_ll = args.lx
            map_y_ll = args.ly
            current_map_res = args.res
            map_ncols = args.ncols
            map_nrows = args.nrows
            x_density_grid_centers = map_x_ll + 0.5 * current_map_res + np.arange(map_ncols) * current_map_res
            y_density_grid_centers = map_y_ll + 0.5 * current_map_res + np.arange(map_nrows) * current_map_res
        else: 
            data_extent_x_min, data_extent_x_max = dem_extent[0], dem_extent[1]
            data_extent_y_min, data_extent_y_max = dem_extent[2], dem_extent[3]
            if grid_definition_mode == "default":
                print(f"Using script default resolution ({DEFAULT_RASTER_RESOLUTION}) and DEM-derived extents for raster maps.")
                current_map_res = DEFAULT_RASTER_RESOLUTION
            else: # "res_only"
                print(f"Using user-defined resolution ({args.res}) and DEM-derived extents for raster maps.")
                current_map_res = args.res
            x_density_grid_centers = np.arange(data_extent_x_min + 0.5 * current_map_res, data_extent_x_max, current_map_res)
            y_density_grid_centers = np.arange(data_extent_y_min + 0.5 * current_map_res, data_extent_y_max, current_map_res)
            if len(x_density_grid_centers) == 0: x_density_grid_centers = np.array([data_extent_x_min + 0.5 * current_map_res])
            if len(y_density_grid_centers) == 0: y_density_grid_centers = np.array([data_extent_y_min + 0.5 * current_map_res])
            map_x_ll, map_y_ll = data_extent_x_min, data_extent_y_min
        
        xx_grid_centers, yy_grid_centers = np.meshgrid(x_density_grid_centers, y_density_grid_centers)
        map_nrows, map_ncols = xx_grid_centers.shape
        map_x_ur = map_x_ll + map_ncols * current_map_res
        map_y_ur = map_y_ll + map_nrows * current_map_res
        density_map_plot_extent = (map_x_ll, map_x_ur, map_y_ll, map_y_ur)

        print('Density map grid definition:')
        print(f'  LL corner: ({map_x_ll:.2f}, {map_y_ll:.2f})')
        print(f'  Cell size: {current_map_res:.2f}')
        print(f'  NCols: {map_ncols}, NRows: {map_nrows}')
        print(f'  UR corner: ({map_x_ur:.2f}, {map_y_ur:.2f})')
        print(f'  Plot extent: {density_map_plot_extent}')
        
        num_plot_classes = len(plot_titles)
        count_ballistic_class = np.zeros((num_plot_classes, map_nrows, map_ncols))

        for xi, yi, di in zip(x_impact_coords, y_impact_coords, diam_impacted):
            grid_ix_frac, grid_jy_frac = (xi - map_x_ll)/current_map_res, (yi - map_y_ll)/current_map_res
            i_col, j_row = int(np.clip(np.ceil(grid_ix_frac)-1,0,map_ncols-1)), int(np.clip(np.ceil(grid_jy_frac)-1,0,map_nrows-1))
            
            count_ballistic_class[-1, j_row, i_col] += 1 

            if num_unique_sizes < 5:
                for k_diam_idx in range(len(diams_for_plot_categories)):
                    if np.isclose(di, diams_for_plot_categories[k_diam_idx]):
                        count_ballistic_class[k_diam_idx, j_row, i_col] += 1; break
            else:
                assigned_to_specific_class = False
                if di <= diams_for_plot_categories[0]:
                    count_ballistic_class[0, j_row, i_col] += 1; assigned_to_specific_class = True
                else:
                    for k_bin_idx in range(len(diams_for_plot_categories) - 1):
                        if di <= diams_for_plot_categories[k_bin_idx + 1]:
                            count_ballistic_class[k_bin_idx + 1, j_row, i_col] += 1
                            assigned_to_specific_class = True; break
                if not assigned_to_specific_class:
                    count_ballistic_class[len(diams_for_plot_categories), j_row, i_col] += 1

        # Loop per il plotting
        for i_class in range(num_plot_classes):
            fig, ax = plt.subplots()
            ax.imshow(ls.hillshade(np.flipud(Z_dem_data), vert_exag=1.0, dx=dem_cell_size, dy=dem_cell_size),
                      cmap='gray', extent=dem_extent)
            ax.set_aspect('equal', 'box')

            zz_counts = np.squeeze(count_ballistic_class[i_class, :, :])
            sum_zz_counts = np.sum(zz_counts)
            
            if sum_zz_counts > 0:
                zz_percentage = zz_counts / sum_zz_counts * 100.0
                zz_log_percentage = np.log10(zz_percentage, out=np.full_like(zz_percentage, -np.inf), where=(zz_percentage > 0))
            else:
                zz_log_percentage = np.full_like(zz_counts, -np.inf)

            finite_vals = zz_log_percentage[np.isfinite(zz_log_percentage)]
            if finite_vals.size > 0:
                min_val_plot = np.min(finite_vals)
                max_val_plot = np.max(finite_vals)
                if np.isclose(max_val_plot, min_val_plot): 
                    max_val_plot = min_val_plot + 0.1 
            else: 
                min_val_plot = 0 # log10(1%)
                max_val_plot = 1 # log10(10%)
            
            # plot_levels sono ancora i valori logaritmici per la mappatura dei colori
            plot_levels = np.linspace(min_val_plot, max_val_plot, 11)
            cmap_terrain = plt.get_cmap('terrain_r') 
            norm_boundary = BoundaryNorm(plot_levels, ncolors=cmap_terrain.N, clip=True)
            
            im_overlay = ax.imshow(np.flipud(zz_log_percentage), cmap=cmap_terrain, norm=norm_boundary,
                                   interpolation='nearest', extent=density_map_plot_extent, alpha=0.65)

            ax.set_xlim(density_map_plot_extent[0], density_map_plot_extent[1])
            ax.set_ylim(density_map_plot_extent[2], density_map_plot_extent[3])

            # --- Modifica per la Colorbar ---
            # I ticks sono nelle posizioni dei livelli logaritmici, ma le etichette mostrano la percentuale
            clb = plt.colorbar(im_overlay, ticks=plot_levels, format=FuncFormatter(log_to_percent_formatter))
            clb.set_label('% Ballistics (Color Scale is Logarithmic)', labelpad=10) # Etichetta aggiornata
            # --- Fine Modifica Colorbar ---

            ax.set_title(plot_titles[i_class])
            
            safe_title_part = plot_titles[i_class].replace(" ", "_").replace("<=", "lte").replace(">", "gt").replace("=","eq").replace(".","p").replace(",","")
            png_file = os.path.join(output_raster_root, f"map_{safe_title_part}_ballistic.png")
            plt.savefig(png_file, dpi=200)
            plt.close(fig)

            nodata_val_asc = -9999.0
            asc_data_to_save = np.copy(zz_log_percentage) # Salviamo ancora i valori logaritmici nel file ASC
            asc_data_to_save[np.isneginf(asc_data_to_save)] = nodata_val_asc

            asc_file = os.path.join(output_raster_root, f"map_{safe_title_part}_ballistic.asc")
            header_asc = f"ncols     {map_ncols}\n"
            header_asc += f"nrows    {map_nrows}\n"
            header_asc += f"xllcorner {map_x_ll:.6f}\n"
            header_asc += f"yllcorner {map_y_ll:.6f}\n"
            header_asc += f"cellsize {current_map_res:.6f}\n"
            header_asc += f"NODATA_value {nodata_val_asc}\n"
            np.savetxt(asc_file, np.flipud(asc_data_to_save), header=header_asc, fmt='%.5f', comments='')

    print("Processing completed.")

# Assicurati che if __name__ == '__main__': sia alla fine
if __name__ == '__main__':
    main()
