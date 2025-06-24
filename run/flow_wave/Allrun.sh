foamCleanCase

cd preprocessing
python createSTL.py > log.createSTL
echo "createSTL completed"
cd ..

cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

blockMesh > log.blockMesh
echo "blockMesh completed"

checkMesh -allTopology -allGeometry > log.checkMesh0
echo "checkMesh0 completed"

snappyHexMesh -overwrite > log.snappyHexMesh
echo "snappyHexMesh completed"

extrudeMesh > log.extrudeMesh
echo "snappyHexMesh completed"

checkMesh -allTopology -allGeometry > log.checkMesh1
echo "checkMesh1 completed"


changeDictionary > log.changeDictionary
echo "changeDictionary completed"

rm -rf 0
cp -r org.0 0
mv 0/alpha.air.init 0/alpha.air
mv 0/alpha.particles.init 0/alpha.particles
mv 0/T.air.init 0/T.air
mv 0/T.particles.init 0/T.particles
mv 0/U.air.init 0/U.air
mv 0/U.particles.init 0/U.particles

foamRun > log.foamRun0
echo "foamRun0 completed"

mv 0/alpha.air.run 0/alpha.air
mv 0/alpha.particles.run 0/alpha.particles
mv 0/T.air.run 0/T.air
mv 0/T.particles.run 0/T.particles
mv 0/U.air.run 0/U.air
mv 0/U.particles.run 0/U.particles

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution

foamRun > log.foamRun1
echo "foamRun1 completed"



