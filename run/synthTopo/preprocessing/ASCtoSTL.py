from matplotlib import pyplot
from shapely.geometry import LineString, Point
from shapely.affinity import translate
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from stl import mesh
from scipy import interpolate
from linecache import getline
import sys
import os.path
import time
import pandas as pd

from ASCtoSTLdict import *



def saveDicts(xmin,xmax,ymin,ymax,Zinit,delta_mesh,path):

    # Define the 3D block and compute the number of cells in each direction
    nx = int(np.ceil((xmax-xmin)/delta_mesh))
    ny = int(np.ceil((ymax-ymin)/delta_mesh))
    zmin = np.min(np.min(Zinit)) - offset_mesh
    zmax = np.max(np.max(Zinit)) + z_atm
    nz = int(np.ceil((zmax-zmin)/delta_mesh))
    
    # In order to write blockMeshDict dictionary, write the vertices of the
    # block and the discretization
    fid1 = open('blockMesh.mod','w')
    fid1.write('vertices\n')
    fid1.write('(\n')
    fid1.write('(%f %f %f) \n' % (xmin,ymin,zmin))
    fid1.write('(%f %f %f) \n' % (xmax,ymin,zmin))
    fid1.write('(%f %f %f) \n' % (xmax,ymax,zmin))
    fid1.write('(%f %f %f) \n' % (xmin,ymax,zmin))
    fid1.write('(%f %f %f) \n' % (xmin,ymin,zmax))
    fid1.write('(%f %f %f) \n' % (xmax,ymin,zmax))
    fid1.write('(%f %f %f) \n' % (xmax,ymax,zmax))
    fid1.write('(%f %f %f) \n' % (xmin,ymax,zmax))
    fid1.write(');\n\n')
    fid1.write('blocks\n')
    fid1.write('(\n')
    fid1.write('    hex (0 1 2 3 4 5 6 7) (%d %d %d) simpleGrading (1 1 1) \n' % (nx,ny,nz))
    fid1.write(');\n\n')
    fid1.close()

    # Write blockMeshDict dictionary by concatenating the header with the other
    # parts. Save it as path/system/blockMeshDict
    path_system = path + 'system/'
    command = "cat templates/blockMesh.start blockMesh.mod templates/blockMesh.end > " + path_system + "blockMeshDict"


    # Call the operating system to execute the specified commands
    os.system(command)
    os.system('rm blockMesh.mod')
    
    fout = open((path_system + 'snappyHexMeshDict'),'w')
    with open('./templates/snappyHexMeshDict.template') as f:
	    content = f.readlines()
	    textdata = [x.strip() for x in content]
    for i in range(len(textdata)): 
        textdata[i] = textdata[i].replace('xxxxxx', str(zmax-10.0))
        fout.write('%s\n' % textdata[i])
    fout.close()
    
# Print iterations progress
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
        decimals    - Optional  : positive number of decimals in percent complete (Int)
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


line = LineString(points)

print('Reading DEM file: ' + DEM_file)
# Parse the topography header
hdr = [getline(DEM_file, i) for i in range(1, 7)]
values = [float(h.split()[-1].strip()) for h in hdr]
cols, rows, lx, ly, cell, nd = values
cols = int(cols)
rows = int(rows)

xs_DEM = lx + 0.5 * cell + np.linspace(0, (cols - 1) * cell, cols)
ys_DEM = ly + 0.5 * cell + np.linspace(0, (rows - 1) * cell, rows)

# Load the topography into a numpy array
DEM = pd.read_table(DEM_file, delim_whitespace=True, header=None,
                    skiprows=6).astype(float).values
DEM = np.flipud(DEM)
DEM[DEM == nd] = 0.0

xinit = np.linspace(0, (cols - 1) * cell, cols) - xc
yinit = np.linspace(0, (rows - 1) * cell, rows) - yc

xinit = xs_DEM - xc
yinit = ys_DEM - yc
 
xmin = np.amin(xinit)+offset_mesh
xmax = np.amax(xinit)-offset_mesh

print('xmin,xmax',xmin,xmax)

ymin = np.amin(yinit)+offset_mesh
ymax = np.amax(yinit)-offset_mesh

print('ymin,ymax',ymin,ymax)

Xinit, Yinit = np.meshgrid(xinit, yinit)
Zinit = DEM

if 'domain_size_x' in locals(): 

    xmin = np.maximum(xmin,-0.5 * domain_size_x)
    xmax = np.minimum(xmax,0.5 * domain_size_x)

if 'domain_size_y' in locals(): 

    ymin = np.maximum(ymin,-0.5 * domain_size_y)
    ymax = np.minimum(ymax,0.5 * domain_size_y)

if saveDicts_flag:

    saveDicts(xmin,xmax,ymin,ymax,Zinit,delta_mesh,path)

# translate the linestring (relative reference system with (lx,ly)=(0,0))
line = translate(line, -xc, -yc)

bb = line.bounds

print('bb',bb)

# interpolate with values from original fine grid
f = interpolate.interp2d(xinit, yinit, Zinit, kind='linear')

# coarsening of the original grid
xinit = xinit[::resample]
yinit = yinit[::resample]

Xinit = Xinit[::resample, ::resample]
Yinit = Yinit[::resample, ::resample]
Zinit = Zinit[::resample, ::resample]

# list of points of DEM for the refined nested grid
x_check = []
y_check = []
z_check = []

# list of displacements
dx = []
dy = []
dz = []

# boundaing box of first refinement level
xmin_bb = bb[0] - nlevels * dist0
ymin_bb = bb[1] - nlevels * dist0
xmax_bb = bb[2] + nlevels * dist0
ymax_bb = bb[3] + nlevels * dist0

# x bounding box indexes
idx_min = np.searchsorted(xinit, xmin_bb, side="left") - 1
idx_max = np.searchsorted(xinit, xmax_bb, side="left")

x0 = xinit[idx_min:idx_max + 1]

# y bounding box indexes
idy_min = np.searchsorted(yinit, ymin_bb, side="left") - 1
idy_max = np.searchsorted(yinit, ymax_bb, side="left")

y0 = yinit[idy_min:idy_max + 1]

# add original points on the west of first bounding box
x_check.extend(Xinit[:, 0:idx_min].ravel())
y_check.extend(Yinit[:, 0:idx_min].ravel())
z_check.extend(Zinit[:, 0:idx_min].ravel())

dx.extend(0.0 * Xinit[:, 0:idx_min].ravel())
dy.extend(0.0 * Xinit[:, 0:idx_min].ravel())
dz.extend(0.0 * Xinit[:, 0:idx_min].ravel())

# add original points on the north of first bounding box
x_check.extend(Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
y_check.extend(Yinit[0:idy_min, idx_min:idx_max + 1].ravel())
z_check.extend(Zinit[0:idy_min, idx_min:idx_max + 1].ravel())

dx.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
dy.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
dz.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())

# add original points on the south of first bounding box
x_check.extend(Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
y_check.extend(Yinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
z_check.extend(Zinit[idy_max + 1:, idx_min:idx_max + 1].ravel())

dx.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
dy.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
dz.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())

# add original points on the east of first bounding box
x_check.extend(Xinit[:, idx_max + 1:].ravel())
y_check.extend(Yinit[:, idx_max + 1:].ravel())
z_check.extend(Zinit[:, idx_max + 1:].ravel())

dx.extend(0.0 * Xinit[:, idx_max + 1:].ravel())
dy.extend(0.0 * Xinit[:, idx_max + 1:].ravel())
dz.extend(0.0 * Xinit[:, idx_max + 1:].ravel())

dist_lev0 = 1.e10

# loop over revinement levels
for ilevel in range(nlevels):

    print('')
    print('Refinement level:', ilevel)

    X0, Y0 = np.meshgrid(x0, y0)

    X0_1d = X0.ravel()
    Y0_1d = Y0.ravel()

    nxy = len(X0_1d)

    dist_lev = dist0 * (nlevels - ilevel)
    print('Distance range:', dist_lev0, dist_lev)
    print('Points to check:', nxy)

    i = 0

    for j, (x, y) in enumerate(zip(X0_1d, Y0_1d)):

        printProgressBar(np.asarray(i) + 1,
                         nxy,
                         prefix='Progress:',
                         suffix='Complete',
                         decimals=1,
                         bar_length=50)

        dist = line.distance(Point(x, y))

        if (dist > dist_lev) and (dist <= dist_lev0):

            dx.append(0.0)
            dy.append(0.0)
            dz.append(0.0)

            x_check.append(x)
            y_check.append(y)

            z = f(x, y)

            z_check.append(float(z))

        i += 1
        
    # bounding box of new refinement level
    xmin_bb = bb[0] - dist_lev
    ymin_bb = bb[1] - dist_lev
    xmax_bb = bb[2] + dist_lev
    ymax_bb = bb[3] + dist_lev

    idx_min = np.searchsorted(x0, xmin_bb, side="left") - 1
    idx_max = np.searchsorted(x0, xmax_bb, side="left")
    nx = 2 * (idx_max - idx_min) + 1

    x1 = np.linspace(x0[idx_min], x0[idx_max], nx)

    idx_min = np.searchsorted(y0, ymin_bb, side="left") - 1
    idx_max = np.searchsorted(y0, ymax_bb, side="left")
    ny = 2 * (idx_max - idx_min) + 1
    y1 = np.linspace(y0[idx_min], y0[idx_max], ny)

    x0 = x1
    y0 = y1
    dist_lev0 = dist_lev

print('')
print('Refinement level:', nlevels)

X1, Y1 = np.meshgrid(x1, y1)

X1_1d = X1.ravel()
Y1_1d = Y1.ravel()

nxy = len(X1_1d)
print('Distance range:', dist_lev0, 0)
print('Points to check:', nxy)

n_out = len(x_check)

i = 0

for j, (x, y) in enumerate(zip(X1_1d, Y1_1d)):

    printProgressBar(np.asarray(i) + 1,
                     nxy,
                     prefix='Progress:',
                     suffix='Complete',
                     decimals=1,
                     bar_length=50)
    dist = line.distance(Point(x, y))

    if (dist <= dist_lev):

        dz_rel = -(1.0 - (dist / dist0)**enne)**(1.0 / enne)

        dx.append(dz_rel * xb)
        dy.append(dz_rel * yb)
        dz.append(depth * dz_rel)

        z = f(x, y)

        x_check.append(x)
        y_check.append(y)
        z_check.append(float(z))

    i += 1

x_old = np.array(x_check)
y_old = np.array(y_check)
z_old = np.array(z_check)

x_new = np.array(x_check) + np.array(dx)
y_new = np.array(y_check) + np.array(dy)
z_new = np.array(z_check) + np.array(dz)

points = []

for i, (x, y) in enumerate(zip(x_check, y_check)):

    points.append([x, y])

points = np.asarray(points)

print('')
print('Building Delaunay triangulation')

tri = Delaunay(points)

print('Delaunay triangulation completed')

# split the triangulation
inner_tri_list = []

tri_simpl = tri.simplices.copy()

vert = range(n_out, len(x_check) + 1)

inner_tri_list = np.in1d(tri_simpl[:, 0], vert)
inner_tri_list = np.logical_or(inner_tri_list, np.in1d(tri_simpl[:, 1], vert))
inner_tri_list = np.logical_or(inner_tri_list, np.in1d(tri_simpl[:, 2], vert))

inner_tri_list = np.arange(tri_simpl.shape[0])[inner_tri_list]

tri_in = tri_simpl[inner_tri_list, :]

outer_tri_list = list(set(range(tri_simpl.shape[0])) - set(inner_tri_list))

tri_out = tri_simpl[outer_tri_list, :]

# original elevation 3D point refined grid
vertices_old = np.column_stack((x_old, y_old, z_old))

# modified elevation 3D point refined grid
vertices = np.column_stack((x_new, y_new, z_new))

# Create the full mesh
print('')
print('Saving full stl')
faces = tri.simplices

surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        surface.vectors[i][j] = vertices[f[j], :]

output_dir = '../constant/triSurface'
# Check whether the specified output path exists or not
isExist = os.path.exists(output_dir)

if not isExist:

    # Create a new directory because it does not exist
    os.makedirs(output_dir)
    print('The new directory ' + output_dir + ' is created!')

surface.save('../constant/triSurface/surface.stl')

# Create the inside mesh
print('Saving in stl')

faces = tri_in

surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        surface.vectors[i][j] = vertices[f[j], :]

surface.save('../constant/triSurface/surface_in.stl')

# Create the outside mesh
print('Saving out stl')
faces = tri_out

surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        surface.vectors[i][j] = vertices[f[j], :]

surface.save('../constant/triSurface/surface_out.stl')

# Create the inside mesh closed on top
print('Saving closed in stl')

faces = tri_in

surface = mesh.Mesh(np.zeros(2 * faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        surface.vectors[i][j] = vertices[f[j], :]
        surface.vectors[i + faces.shape[0]][j] = vertices_old[f[j], :]

surface.save('../constant/triSurface/surface_in_closed.stl')

# plt.show()
