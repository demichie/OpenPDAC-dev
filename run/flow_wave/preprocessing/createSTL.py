import numpy as np

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
import stl
    
from createSTLDict import *

def combined_stl(paths, save_path="./combined.stl"):

    with open(save_path,'w') as outfile:
        for fname in paths:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line.replace('surface.stl','groundAndWall'))

def saveAsc(x,y,h,DEM_file):

    
    dx = x[1]-x[0]
    xmin = x[0] 
    
    dy = y[1]-y[0]
    ymin = y[0] 

    # Save initial topography on ascii raster file
    header = "ncols     %s\n" % h.shape[1]
    header += "nrows    %s\n" % h.shape[0]
    header += "xllcenter " + str(xmin) + "\n"
    header += "yllcenter " + str(ymin) + "\n"
    header += "cellsize " + str(dx) + "\n"
    header += "NODATA_value -9999\n"

    resampled_DEM = DEM_file.replace('.asc','_resampled.asc')

    print('')
    print('Saving resample DEM:',resampled_DEM)

    np.savetxt(resampled_DEM,
               np.flipud(h),
               header=header,
               fmt='%1.5f',
               comments='')

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
    
    
def main():

    cellsize = 10.0

    xmin = 0.0
    xmax = 1100.0
    
    x = np.arange(xmin,xmax,cellsize)
    
    ymin = 0.0
    ymax = 300.0
    
    y = np.arange(ymin,ymax,cellsize)
    
    X,Y = np.meshgrid(x,y)
        
    W = wave_amplitude * np.sin( np.maximum(0.0,X-wave_start) / wave_period )
    Z = np.tan(np.deg2rad(slope))*X + W
    
    X_1d = X.ravel()
    Y_1d = Y.ravel()
    Z_1d = Z.ravel()

    nxy = len(X_1d)
    
    points = []

    for i, (x, y) in enumerate(zip(X_1d, Y_1d)):

        points.append([x, y])

    points = np.asarray(points)
    
    print('')
    print('Building Delaunay triangulation')

    tri = Delaunay(points)

    print('Delaunay triangulation completed')

    # Create the full mesh
    print('')
    print('Saving full stl')
    faces = tri.simplices

    # modified elevation 3D point refined grid
    # vertices = np.column_stack((X_1d, Y_1d, Z_1d))
    vertices = np.column_stack((X_1d, Z_1d, Y_1d))

    groundAndWall = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            groundAndWall.vectors[i][j] = vertices[f[j], :]

    groundAndWall.save('../constant/triSurface/partial/surface.stl', mode=stl.Mode.ASCII)

    direc = "../constant/triSurface/partial/"
    paths = [os.path.join(direc, i) for i in os.listdir(direc)]
    combined_stl(paths,save_path="../constant/triSurface/combined.stl")

if __name__ == '__main__':
    main()       

