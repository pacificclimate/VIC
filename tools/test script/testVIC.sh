#!/bin/bash
echo starting test

codeDir="VIC_4.1.2_cpp_trunk"
programName="vicNl"
globalOptionsFile="glb_prb_base_BASIN_SCENARIO_19502006_VIC4.1.2_netcdf_auto.txt"
export curDir=$(pwd)
TIMESTAMP=$(date +"%Y_%m_%d__%H_%M_%S")
export outputName="automated_4.1.2_netcdf_$TIMESTAMP"
echo "Output for this test will be in out/$outputName"
echo ""

fail()
{
    echo "automatic test failed"
    exit 1 #exit shell script
}

#Create directory for output
mkdir out/$outputName

#Replace output line in global options file
perl -pi.bak -e 's/RESULT_DIR.*$/RESULT_DIR\t$ENV{curDir}\/out\/$ENV{outputName}/g' $globalOptionsFile

#Make the program
pushd $codeDir > /dev/null
echo "building the program"
make clean
rm $programName
make
#check that it built correctly
if [ -f "$programName" ]
then
    echo "Program built correctly"
else
    echo "Program not built - check output for compile errors"
    fail
fi
popd > /dev/null

#Run the code
echo "running the program"
$codeDir/$programName -g $globalOptionsFile

#Compare output files with pristine version hashes
pushd out/$outputName > /dev/null
echo ""
echo $(md5sum -c ../4.1.2_pristine_forcings_v2_1950-2006.md5 2>/dev/null | grep -cE 'OK$') out of $(cat ../4.1.2_pristine_forcings_v2_1950-2006.md5 | wc -l) files compare OK

popd > /dev/null

echo "Finished."