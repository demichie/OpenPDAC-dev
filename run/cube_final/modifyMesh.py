import numpy as np
import sys
from linecache import getline
import pandas as pd
from scipy.interpolate import RegularGridInterpolator as RGI
from stl import mesh
from scipy.interpolate import RBFInterpolator
import argparse
import struct


def create_crater_stl(rbf,
                      x,
                      y,
                      dx_base,
                      dy_base,
                      z_base,
                      r_bottom,
                      r_top,
                      h,
                      n=100,
                      filename='crater.stl'):
    # Define the number of points around the circumference
    theta = np.linspace(0, 2 * np.pi, n, endpoint = False)

    # Generate points for the top and bottom circles
    x_circle_bottom = r_bottom * np.cos(theta) + dx_base
    y_circle_bottom = r_bottom * np.sin(theta) + dy_base
    x_circle_top = r_top * np.cos(theta)
    y_circle_top = r_top * np.sin(theta)

    z_circle_top = np.zeros_like(x_circle_bottom)
    z_circle_bottom = np.zeros_like(x_circle_bottom)

    for i, (xt, yt) in enumerate(zip(x_circle_top + x, y_circle_top + x)):

        z_circle_top[i] = rbf([[xt, yt]])

    for i, (xb, yb) in enumerate(zip(x_circle_bottom + x,
                                     y_circle_bottom + x)):

        z_circle_bottom[i] = rbf([[xb - dx_base, yb - dy_base]]) + z_base

    top_circle = np.column_stack((x_circle_top + x, y_circle_top + y,
                                  z_circle_top))
    bottom_circle = np.column_stack((x_circle_bottom + x, y_circle_bottom + y,
                                     z_circle_bottom))

    # Create faces for the sides of the crater
    faces = []
    for i in range(n):
        next_i = (i + 1) % n
        # Side faces
        faces.append(
            [bottom_circle[i], bottom_circle[next_i], top_circle[next_i]])
        faces.append([bottom_circle[i], top_circle[next_i], top_circle[i]])

    z_center_top = rbf([[x, y]])[0]
    z_center_bottom = rbf([[x, y]])[0]

    # Create faces for the top and bottom circles
    center_top = np.array([x, y, z_base + h + z_center_top])
    center_bottom = np.array(
        [x + dx_base, y + dy_base, z_base + z_center_bottom])

    for i in range(n):
        next_i = (i + 1) % n
        # Top face
        faces.append([center_top, top_circle[i], top_circle[next_i]])
        # Bottom face
        faces.append([center_bottom, bottom_circle[next_i], bottom_circle[i]])

    # Create the mesh
    faces_np = np.array(faces)
    crater_mesh = mesh.Mesh(np.zeros(faces_np.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces_np):
        for j in range(3):
            crater_mesh.vectors[i][j] = f[j]

    volume, cog, inertia = crater_mesh.get_mass_properties()
    print("Crater Volume = {0}".format(volume))
    print("Closed", crater_mesh.is_closed())

    # Save the mesh to file
    crater_mesh.save(filename)


def create_conduit_stl(rbf,
                       x,
                       y,
                       dx_base,
                       dy_base,
                       z_base,
                       r,
                       h,
                       n=100,
                       filename='conduit.stl'):
    # Define the number of points around the circumference
    theta = np.linspace(0, 2 * np.pi, n, endpoint = False)

    # Generate points for the top and bottom circles
    x_circle = r * np.cos(theta)
    y_circle = r * np.sin(theta)

    z_circle = np.zeros_like(x_circle)

    for i, (xt, yt) in enumerate(zip(x_circle, y_circle)):

        z_circle[i] = rbf([[xt, yt]])

    top_circle = np.column_stack(
        (x_circle + x, y_circle + y, np.full(n, z_base + h) + z_circle))
    bottom_circle = np.column_stack(
        (x_circle + x + dx_base, y_circle + y + dy_base,
         np.full(n, z_base) + z_circle))

    # Create faces for the sides of the conduit
    faces = []
    for i in range(n):
        next_i = (i + 1) % n
        # Side faces
        faces.append(
            [bottom_circle[i], bottom_circle[next_i], top_circle[next_i]])
        faces.append([bottom_circle[i], top_circle[next_i], top_circle[i]])

    z_center = rbf([[0, 0]])[0]

    # Create faces for the top and bottom circles
    center_top = np.array([x, y, z_base + h + z_center])
    center_bottom = np.array([x + dx_base, y + dy_base, z_base + z_center])

    for i in range(n):
        next_i = (i + 1) % n
        # Top face
        faces.append([center_top, top_circle[i], top_circle[next_i]])
        # Bottom face
        faces.append([center_bottom, bottom_circle[next_i], bottom_circle[i]])

    # Create the mesh
    faces_np = np.array(faces)
    conduit_mesh = mesh.Mesh(np.zeros(faces_np.shape[0],
                                      dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces_np):
        for j in range(3):
            conduit_mesh.vectors[i][j] = f[j]

    volume, cog, inertia = conduit_mesh.get_mass_properties()
    print("Conduit Volume = {0}".format(volume))
    print("Closed", conduit_mesh.is_closed())


    # Save the mesh to file
    conduit_mesh.save(filename)


def generate_circumference_points(n, r):
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = r * np.cos(angles)
    y = r * np.sin(angles)
    points = np.column_stack((x, y))
    return points


def format_points(points):
    return [f"({x} {y} {z})" for x, y, z in points]


def printProgressBar(iteration,
                     total,
                     prefix='',
                     suffix='',
                     decimals=1,
                     bar_length=100):
    """
    Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals
                                  in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = 'X' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' %
                     (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def is_integer(s):
    try:
        x = int(s)
        return True
    except ValueError:
        return False


def is_binary_format(content, maxline=20):
    """
    parse file header to judge the format is binary or not
    :param content: file content in line list
    :param maxline: maximum lines to parse
    :return: binary format or not
    """
    for lc in content[:maxline]:
        if b'format' in lc:
            if b'binary' in lc:
                return True
            return False
    return False


def parse_points_content(content, is_binary, skip=10):
    """
    parse points from content
    :param content: file contents
    :param is_binary: binary format or not
    :param skip: skip lines
    :return: points coordinates as numpy.array
    """
    n = skip
    while n < len(content):
        lc = content[n]
        if is_integer(lc):
            num = int(lc)
            if not is_binary:
                data = np.array(
                    [ln[1:-2].split() for ln in content[n + 2:n + 2 + num]],
                    dtype=float)
            else:
                buf = b''.join(content[n + 1:])
                disp = struct.calcsize('c')
                vv = np.array(
                    struct.unpack(
                        '{}d'.format(num * 3),
                        buf[disp:num * 3 * struct.calcsize('d') + disp]))
                data = vv.reshape((num, 3))
            return data, n
        n += 1
    return None


def main():

    msg = "Python script to modify the grid and create STL surfaces"

    # Initialize parser
    parser = argparse.ArgumentParser(description = msg)
    
    # Adding a boolean flag --verbose that will be set to True if provided
    parser.add_argument('-s','--STLonly', action='store_true', help='save only STL files')
    
    # Parse the arguments
    args = parser.parse_args()
        
    parser.parse_args()

    if args.STLonly:
    
        print ("Create STL files only.")
        onlySTL_flag = True
        
    else:   

        print ("Create STL files and modify points.")
        onlySTL_flag = False

    try:

        from modifyMeshDict import h_crater

    except ImportError:

        print('Missing parameter in dict: h_crater (float)')
        sys.exit(1)

    try:

        from modifyMeshDict import r_crater_top

    except ImportError:

        print('Missing parameter in dict: r_crater_top (float)')
        sys.exit(1)

    try:

        from modifyMeshDict import r_crater_bottom

    except ImportError:

        print('Missing parameter in dict: r_crater_bottom (float)')
        sys.exit(1)

    try:

        from modifyMeshDict import dx_crater_base

    except ImportError:

        print('Missing parameter in dict: dx_crater_base (float)')
        dx_crater_base = 0.0

    try:

        from modifyMeshDict import dy_crater_base

    except ImportError:

        print('Missing parameter in dict: dy_crater_base (float)')
        dy_crater_base = 0.0

    try:

        from modifyMeshDict import h_conduit

    except ImportError:

        print('Missing parameter in dict: h_conduit (float)')
        h_conduit = 0.0
        r_conduit = 0.0

    try:

        from modifyMeshDict import r_conduit

    except ImportError:

        print('Missing parameter in dict: r_conduit (float)')
        h_conduit = 0.0
        r_conduit = 0.0

    try:

        from modifyMeshDict import dx_conduit_base

    except ImportError:

        print('Missing parameter in dict: dx_conduit_base (float)')
        dx_crater_base = 0.0

    try:

        from modifyMeshDict import dy_conduit_base

    except ImportError:

        print('Missing parameter in dict: dy_conduit_base (float)')
        dy_crater_base = 0.0

    try:

        from modifyMeshDict import DEM_file

    except ImportError:

        print('Missing parameter in dict: DEM_file (str)')
        sys.exit(1)

    print('Reading DEM file: ' + DEM_file)
    # Parse the topography header
    hdr = [getline(DEM_file, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    # values are associated to the centers of the DEM pixels
    xs_DEM = lx + 0.5 * cell + np.linspace(0, (cols - 1) * cell, cols)
    ys_DEM = ly + 0.5 * cell + np.linspace(0, (rows - 1) * cell, rows)

    # Load the topography into a numpy array
    DEM = pd.read_table(DEM_file,
                        delim_whitespace=True,
                        header=None,
                        skiprows=6).astype(float).values
    DEM = np.flipud(DEM)
    DEM[DEM == nd] = 0.0

    try:

        from modifyMeshDict import xc

    except ImportError:

        print('Missing parameter in dict: xc (float)')
        sys.exit(1)

    try:

        from modifyMeshDict import yc

    except ImportError:

        print('Missing parameter in dict: yc (float)')
        sys.exit(1)

    xinit = xs_DEM - xc
    yinit = ys_DEM - yc
    Zinit = DEM

    try:

        from modifyMeshDict import dzVert

    except ImportError:

        print('Missing parameter in dict: dzVert (float)')
        dzVert = 0.0

    zVert = np.amax(Zinit) + dzVert

    # bilinear interpolation with values from original fine grid
    f = RGI((xinit, yinit), Zinit.T, method='linear', bounds_error=False)

    # create a list of points on the edge of the crater for a radial basis
    # interpolation to smooth the topography
    n_top_points = 100

    angles = np.linspace(0, 2 * np.pi, n_top_points, endpoint=False)
    xedge = r_crater_top * np.cos(angles)
    yedge = r_crater_top * np.sin(angles)

    zedge = np.zeros_like(xedge)

    for i in range(np.shape(zedge)[0]):

        zedge[i] = f([[xedge[i], yedge[i]]])

    xedge = xedge.reshape(-1, 1)
    yedge = yedge.reshape(-1, 1)

    xy = np.concatenate([xedge, yedge], axis=1)
    rbf = RBFInterpolator(xy, zedge, epsilon=2)

    # number of points for the STL surfaces
    npoints = 100

    create_crater_stl(rbf,
                      x=0,
                      y=0,
                      dx_base=dx_crater_base,
                      dy_base=dy_crater_base,
                      z_base=-(h_crater),
                      r_bottom=r_crater_bottom,
                      r_top=r_crater_top,
                      h=h_crater,
                      n=npoints,
                      filename='./constant/triSurface/crater.stl')

    create_conduit_stl(rbf,
                       x=dx_crater_base,
                       y=dy_crater_base,
                       dx_base=dx_conduit_base,
                       dy_base=dy_conduit_base,
                       z_base=-(h_conduit + h_crater),
                       r=r_conduit,
                       h=h_conduit,
                       n=npoints,
                       filename='./constant/triSurface/conduit.stl')

    # list of points generated by blockMesh
    fn = "./constant/polyMesh/points"

    try:
        with open(fn, "rb") as pointfile:
            content = pointfile.readlines()
            points, nheader = parse_points_content(content,
                                                   is_binary_format(content),
                                                   skip=10)

    except FileNotFoundError:
        print('file not found: %s' % fn)
        sys.exit(1)

    xB = points[:, 0]
    yB = points[:, 1]
    zB = points[:, 2]

    xBmin = np.min(xB)-0.1
    xBmax = np.max(xB)+0.1

    yBmin = np.min(yB)-0.1
    yBmax = np.max(yB)+0.1

    zBmin = np.min(zB)
    zBmax = np.max(zB)

    zmin = -100.0
    zmax = 10000.0
    xmin = xBmin
    xmax = xBmax
    ymin = yBmin
    ymax = yBmax

    if onlySTL_flag:

        try:

            from modifyMeshDict import zRefine
            zmax = zRefine
            zmin = zBmin-1.0

        except ImportError:

            print('Missing parameter in dict: zRefine (float)')

    # Define the 8 vertices of the cube
    vertices = np.array([[xmin, ymin, zmin], [xmax, ymin, zmin],
                         [xmax, ymax, zmin], [xmin, ymax, zmin],
                         [xmin, ymin, zmax], [xmax, ymin, zmax],
                         [xmax, ymax, zmax], [xmin, ymax, zmax]])
    # Define the 12 triangles composing the cube
    faces = np.array([[0, 3, 1], [1, 3, 2], [0, 4, 7], [0, 7, 3], [4, 5, 6],
                      [4, 6, 7], [5, 1, 2], [5, 2, 6], [2, 3, 6], [3, 7, 6],
                      [0, 1, 5], [0, 5, 4]])

    # Create the mesh
    cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, face in enumerate(faces):
        for j in range(3):
            cube.vectors[i][j] = vertices[face[j]]

    if onlySTL_flag:

        # Write the mesh to file "cube.stl"
        cube.save('./constant/triSurface/refine.stl')

    else:
    
        # Write the mesh to file "cube.stl"
        cube.save('./constant/triSurface/cube.stl')
    

    if onlySTL_flag:
    
        return

    pointsNew = np.copy(points)
    npoints = points.shape[0]

    try:

        from modifyMeshDict import exp_factor

    except ImportError:

        print('Missing parameter in dict: exp_factor (float)')
        exp_factor = 0.0

    try:

        from modifyMeshDict import exp_shape

    except ImportError:

        print('Missing parameter in dict: exp_shape (float)')
        exp_factor = 1.0


    for i, point in enumerate(points):

        x = point[0]
        y = point[1]
        z = point[2]

        if (x**2 + y**2 > r_crater_top**2):

            zInterp = f([[x, y]])

        else:

            zInterp = rbf([[x, y]])

        zRel = (zBmax - z) / zBmax
        zRel = np.minimum(1.0,np.maximum(0.0, zRel))

        zNew = z + zRel * zInterp
        pointsNew[i, 2] = zNew

        if (z >= 0.0):

            # the points above the topography are enlarged
            # horizontally with increasing z

            if (dzVert > 0):

                # enlarge from a fixed height above the maximum
                # topography and the top, thus from an horizontal
                # plane to the top
                z2Rel = np.maximum(0, (zNew - zVert) / (zBmax - zVert))

            elif (dzVert < 0):

                # enlarge from a fixed height above the topography
                # and the top, thus from a terrain following surface
                # to the top
                zVert = zInterp - dzVert

                z2Rel = np.maximum(0, (zNew - zVert) / (zBmax - zVert))

            else:

                # enlarge from the topography to the top
                z2Rel = (zNew - zInterp) / (zBmax - zInterp)

            z2Rel = z2Rel**exp_shape

            xNew = x * (1.0 + z2Rel * (exp_factor - 1.0))
            yNew = y * (1.0 + z2Rel * (exp_factor - 1.0))

        else:

            # the points in the crater and the conduit are
            # translated horizontally

            xNew = x
            yNew = y

            # z_rel_crater is 0 at top of crater and 1 at the base
            z_rel_crater = (z + h_crater) / h_crater
            z_rel_crater = np.maximum(0.0, np.minimum(1.0, z_rel_crater))
            # z_rel_conduit is 0 at top of conduit and 1 at the base
            z_rel_conduit = (z + h_crater + h_conduit) / h_conduit
            z_rel_conduit = np.maximum(0.0, np.minimum(1.0, z_rel_conduit))

            xNew = x + z_rel_crater * dx_crater_base + z_rel_conduit * dx_conduit_base
            yNew = y + z_rel_crater * dy_crater_base + z_rel_conduit * dy_conduit_base

        pointsNew[i, 0] = xNew
        pointsNew[i, 1] = yNew

        printProgressBar(np.asarray(i) + 1,
                         npoints,
                         prefix='Progress:',
                         suffix='Complete',
                         decimals=1,
                         bar_length=50)

    # Find the line with the number of points and insert the formatted points
    num_points = len(pointsNew)
    formatted_points = format_points(pointsNew)

    # Read the content of the file
    with open(fn, 'r') as file:
        lines = file.readlines()

    # Create the new file
    new_file_name = './constant/polyMesh/points.new'
    with open(new_file_name, 'w') as new_file:
        # Write the first n lines
        for i in range(nheader + 2):
            new_file.write(lines[i])

        # Write the formatted points
        for point in formatted_points:
            new_file.write(point + '\n')

        # Write the next 4 lines from the original file
        for i in range(nheader + 2 + num_points, nheader + 2 + num_points + 4):
            new_file.write(lines[i])

    print(f"The new file '{new_file_name}' has been created successfully.")


if __name__ == "__main__":

    main()
