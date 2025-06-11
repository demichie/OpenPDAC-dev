# set fvSolution and cloudProperties for initialization of atm. profile and ballistics
cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

mpirun -np 8 foamRun -parallel > log.foamRun0
echo fisrt foamRun completed

mpirun -np 8 setFields -parallel > log.setFields
echo setFields completed


