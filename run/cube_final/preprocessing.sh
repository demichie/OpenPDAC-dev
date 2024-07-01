foamCleanCase

blockMesh

# NEXT 5 LINES ONLY FOR REFINEMESH
# python modifyMesh.py -s
# topoSet -dict topoSetDict.refine
# refineMesh
# cp -r 0.005/polyMesh constant/
# rm -rf 0.005

checkMesh -allGeometry -allTopology

python modifyMesh.py

cp constant/polyMesh/points constant/polyMesh/points.old
cp constant/polyMesh/points.new constant/polyMesh/points

# CHECK FOR CLOSED SURFACES
#surfaceCheck constant/triSurface/crater.stl
#surfaceCheck constant/triSurface/conduit.stl

checkMesh -allGeometry -allTopology

surfaceFeatures
snappyHexMesh -overwrite
checkMesh -allGeometry -allTopology

#foamToVTK -cellSet underdeterminedCells -fields '(none)'
#foamToVTK -cellSet concaveCells -fields '(none)'
#foamToVTK -faceSet concaveFaces -fields '(none)'
#foamToVTK -faceSet skewFaces -fields '(none)'
#foamToVTK -faceSet lowQualityTetFaces -fields '(none)'

topoSet

# set fvSolution and cloudProperties for initialization of atm. profile and ballistics
cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

# remove previous 0 folder and copy not-initialized fields
rm -rf 0
cp -r org.0 0

# (first run for initialization of the solution: fileds and ballistics) 


#FOR PARALLEL RUN:
#sbatch MPIJob_init.script
#squeue

#FOR SCALAR RUN:
foamRun

# set different values in the crater/conduit zone
setFields



# set the run parameters for the actual simulation
cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR PARALLEL RUN ON PC:
decomposePar
mpirun -np xx foamRun -parallel

#FOR SCALAR RUN:
foamRun


