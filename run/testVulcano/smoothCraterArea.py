import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import re

def read_openfoam_dict(file_path):
    """Legge un file OpenFOAM e restituisce un dizionario con i parametri."""
    data = {}
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File non trovato: {file_path}")
    
    with open(file_path, 'r') as file:
        for line in file:
            match = re.match(r'(\w+)\s+([\d\.\-eE]+);', line.strip())
            if match:
                key, value = match.groups()
                data[key] = float(value)
    return data


def read_asc_with_pandas(filepath):
    """Read a raster .asc file using pandas and return the data and metadata."""
    with open(filepath, 'r') as f:
        # Read the metadata (header)
        header = {}
        for _ in range(6):
            line = f.readline()
            key, value = line.strip().split()
            header[key] = float(value) if '.' in value else int(value)

    # Read the grid data as a pandas DataFrame
    data = pd.read_csv(filepath, skiprows=6, header=None, delim_whitespace=True)
    return data.values, header

def write_asc_with_pandas(filepath, data, header):
    """Write a numpy array to a raster .asc file with pandas."""
    with open(filepath, 'w') as f:
        # Write the header
        for key, value in header.items():
            f.write(f"{key} {value}\n")
        # Write the data
        np.savetxt(f, np.flipud(data), fmt='%.2f')

def replace_area_with_interpolation(data, x, y, radius, header):
    """Replace an area in the raster with interpolated values."""
    nrows, ncols = data.shape
    cellsize = header["cellsize"]
    xllcorner = header["xllcorner"]
    yllcorner = header["yllcorner"]

    # Create coordinate grids
    x_coords, y_coords = np.meshgrid(
        np.arange(ncols) * cellsize + xllcorner,
        np.arange(nrows) * cellsize + yllcorner
    )

    # Identify points within the circular area
    mask = ((x_coords - x)**2 + (y_coords - y)**2) <= (radius)**2
    border_mask = ((x_coords - x)**2 + (y_coords - y)**2 > radius**2) & (
        ((x_coords - x)**2 + (y_coords - y)**2 <= (radius+cellsize)**2)
    )

    # Extract border points for interpolation
    border_x = x_coords[border_mask]
    border_y = y_coords[border_mask]
    border_values = data[border_mask]

    xy = np.concatenate([border_x.reshape(-1, 1), border_y.reshape(-1, 1)], axis=1)
    from scipy.interpolate import RBFInterpolator
    rbf = RBFInterpolator(xy, border_values, kernel='linear', epsilon=2)

    masked_x = x_coords[mask]
    masked_y = y_coords[mask]

    z = []

    for x,y in zip(masked_x,masked_y):
    
        z.append(rbf([[x, y]])[0])
        
    # Apply interpolated values to the mask
    new_data = data.copy()
    new_data[mask] = z

    return new_data, mask, border_x, border_y, border_values

def plot_modified_points_3d(data, mask, x_coords, y_coords, border_x, border_y, border_values, output_plot="modified_points_3d.png"):
    """3D plot of the raster data highlighting modified points."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the full terrain
    # ax.plot_surface(x_coords, y_coords, data, cmap='terrain', alpha=0.8, rstride=1, cstride=1)

    # Highlight modified points
    ax.scatter(x_coords[mask], y_coords[mask], data[mask], color='red', s=10, label='Modified Points')

    ax.scatter(border_x, border_y, border_values, color='black', s=10)


    # Labels and title
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Elevation')
    ax.set_title('3D Plot of Elevation Data with Modified Points')
    ax.legend()

    # Save and show
    plt.savefig(output_plot)
    # plt.show()

def main():
    # Input and output file paths
    input_asc = "dsm_vulc_5m.asc"
    output_asc = "output_file.asc"
    output_plot = "modified_points.png"

    # Percorsi dei file
    geometry_file = os.path.join("constant", "geometryParameters")
    topo_grid_file = os.path.join("system", "topoGridDict")

    # Leggi i parametri dai file
    try:
        geometry_data = read_openfoam_dict(geometry_file)
        topo_grid_data = read_openfoam_dict(topo_grid_file)

        # Estrai i parametri desiderati
        r_crater_top = geometry_data.get("r_crater_top")
        x_vent = topo_grid_data.get("xVent")
        y_vent = topo_grid_data.get("yVent")

        print(f"r_crater_top: {r_crater_top}, Xvent: {x_vent}, Yvent: {y_vent}")
    except FileNotFoundError as e:
        print(e)


    # Point coordinates and radius for replacement
    point_x = x_vent  # Example x-coordinate (e.g., UTM or other units)
    point_y = y_vent  # Example y-coordinate
    radius = r_crater_top      # Radius in the same unit as coordinates

    # Read the input raster file
    data, header = read_asc_with_pandas(input_asc)

    data = np.flipud(data)
    
    # Get grid coordinates
    nrows, ncols = data.shape
    cellsize = header["cellsize"]
    xllcorner = header["xllcorner"]
    yllcorner = header["yllcorner"]
    x_coords, y_coords = np.meshgrid(
        np.arange(ncols) * cellsize + xllcorner,
        np.arange(nrows) * cellsize + yllcorner
    )

    # Replace area with interpolated values
    new_data, mask, border_x, border_y, border_values = replace_area_with_interpolation(data, point_x, point_y, radius, header)

    # Plot the modified points in 3D
    plot_modified_points_3d(new_data, mask, x_coords, y_coords, border_x, border_y, border_values, output_plot)

    # Write the modified data to a new raster file
    write_asc_with_pandas(output_asc, new_data, header)

    print(f"Modified raster saved to {output_asc}")

if __name__ == "__main__":
    main()




