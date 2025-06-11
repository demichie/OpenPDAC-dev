cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

mpirun -np 8 foamRun -parallel > log.foamRun &

