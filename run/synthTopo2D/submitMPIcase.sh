rm -rf constant/triSurface/*.stl
foamCleanCase

cd preprocessing
python3 ASCtoSTL.py

cd ..
surfaceCheck constant/triSurface/surface_crater_closed.stl
surfaceCheck constant/triSurface/surface_conduit_closed.stl
surfaceCheck constant/triSurface/surface_total_closed.stl

cp ./system/controlDict.init ./system/controlDict
blockMesh > log.blockMesh
checkMesh -allTopology -allGeometry > log.checkMesh0

snappyHexMesh -overwrite > log.snappyHexMesh
extrudeMesh > log.extrudeMesh
changeDictionary > log.changeDictionary

checkMesh -allTopology -allGeometry > log.checkMesh1

topoSet -dict topoSetDict-conduit > log.topoSet

cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

rm -rf 0
cp -r org.0 0

foamRun > log.foamRun0
setFields > log.setFields

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

foamRun > log.foamRun1 &


