# module load openFOAM-10
# source $FOAM_BASHRC

rm -rf constant/triSurface/*
foamCleanCase

cd preprocessing
unzip synthDEM.zip
python3 ASCtoSTL.py
# python createSphere.py
cd ..
surfaceCheck constant/triSurface/surface_crater_closed.stl
surfaceCheck constant/triSurface/surface_conduit_closed.stl
surfaceCheck constant/triSurface/surface_total_closed.stl

cp ./system/controlDict.init ./system/controlDict
blockMesh 
checkMesh -allTopology -allGeometry

snappyHexMesh -overwrite
topoSet -dict topoSetDict-conduit
# topoSet -dict topoSetDict-sphere

cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

rm -rf 0
cp -r org.0 0

#FOR PARALLEL RUN:
#sbatch MPIJob_init.script
#squeue

#FOR SCALAR RUN:
foamRun
setFields

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR SCALAR RUN:
foamRun

python plotBallistics.py

rm -rf VTK
foamToVTK -fields "(alpha.particles)"

