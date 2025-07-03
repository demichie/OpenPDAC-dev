foamCleanCase

python smoothCraterArea.py

blockMesh > log.blockMesh
echo blockMesh completed


# remove previous 0 folder and copy not-initialized fields
rm -rf 0
cp -r org.0 0

decomposePar > log.decomposePar
echo decomposePar completed

mpirun -np 8 checkMesh -parallel -allGeometry -allTopology -writeSets > log.checkMesh0
echo checkMesh0 completed

mpirun -np 8 topoSet -parallel > log.topoSet
echo topoSet completed

# Grid deformation
mpirun -np 8 topoGrid -parallel > log.topoGrid
echo topoGrid completed

mpirun -np 8 checkMesh -parallel -allGeometry -allTopology -writeSets > log.checkMesh1
echo checkMesh1 completed

