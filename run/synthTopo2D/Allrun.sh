rm -rf constant/triSurface/*.stl
rm -rf 0
foamCleanCase
echo "foamCleanCase completed"

cd preprocessing
python3 ASCtoSTL.py
echo "ASCtoSTL completed"

cd ..
surfaceCheck constant/triSurface/surface_crater_closed.stl > log.surfaceCheck0
surfaceCheck constant/triSurface/surface_conduit_closed.stl > log.surfaceCheck1
surfaceCheck constant/triSurface/surface_total_closed.stl > log.surfaceCheck2

echo "surfaceCheck completed"

cp ./system/controlDict.init ./system/controlDict
blockMesh > log.blockMesh
echo "blockMesh completed"

checkMesh -allTopology -allGeometry > log.checkMesh0
echo "checkMesh completed"

snappyHexMesh -overwrite > log.snappyHexMesh
echo "snappyHexMesh completed"

extrudeMesh > log.extrudeMesh
echo "extrudeMesh completed"

changeDictionary > log.changeDictionary
echo "changeDictionary completed"

checkMesh -allTopology -allGeometry > log.checkMesh1
echo "checkMesh completed"

topoSet -dict topoSetDict-conduit > log.topoSet
echo "topoSet completed"

cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

rm -rf 0
cp -r org.0 0

foamRun > log.foamRun0
echo "foamRun0 completed"

setFields > log.setFields
echo "setFields completed"

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

foamRun > log.foamRun1
echo "foamRun1 completed"


