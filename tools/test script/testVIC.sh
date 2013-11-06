#!/bin/bash

netCDF=false;
stateOutputASCII=false
stateOutputBinary=false
stateOutputNetCDF=false

for var in "$@"
do
    if [ "$var" == "-netCDF" ]; then
        echo "netCDF option enabled";
        netCDF=true;
    fi
    if [ "$var" == "-stateOutputASCII" ]; then
        echo "state output ASCII enabled";
        stateOutputASCII=true;
    fi
    if [ "$var" == "-stateOutputBinary" ]; then
        echo "state output Binary enabled";
        stateOutputBinary=true;
    fi
done

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

#Replace output line in global options file
perl -pi.bak -e 's/RESULT_DIR.*$/RESULT_DIR\t$ENV{curDir}\/out\/$ENV{outputName}/g' $globalOptionsFile

if $netCDF ; then
    perl -pi.bak -e 's/OUTPUT_FORMAT.*$/OUTPUT_FORMAT\tNETCDF/g' $globalOptionsFile
else
    perl -pi.bak -e 's/OUTPUT_FORMAT.*$/OUTPUT_FORMAT\tASCII/g' $globalOptionsFile
fi

if $stateOutputASCII ; then
    perl -pi.bak -e 's/.*STATENAME.*$/STATENAME\tfrs.state/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_OUTPUT_FORMAT.*$/STATE_OUTPUT_FORMAT\tASCII/g' $globalOptionsFile
elif $stateOutputBinary ; then
    perl -pi.bak -e 's/.*STATENAME.*$/STATENAME\tfrs.state/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_OUTPUT_FORMAT.*$/STATE_OUTPUT_FORMAT\tBINARY/g' $globalOptionsFile
else 
    perl -pi.bak -e 's/.*STATENAME.*$/#STATENAME\tfrs.state/g' $globalOptionsFile
fi

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

#Create directory for output
mkdir out/$outputName

#Run the code
echo "running the program"
$codeDir/$programName -g $globalOptionsFile

#Compare output files with pristine version hashes
echo ""
echo "validating output results"
if $netCDF ; then
    Rscript vic_output_netcdf_compare.r out/$outputName/results.nc vicNetCDFCorrectOutputs.nc
else
    pushd out/$outputName > /dev/null
    echo $(md5sum -c ../4.1.2_pristine_forcings_v2_1950-2006.md5 2>/dev/null | grep -cE 'OK$') out of $(cat ../4.1.2_pristine_forcings_v2_1950-2006.md5 | wc -l) files compare OK
    popd > /dev/null
fi

if $stateOutputASCII ; then
    md5sum frs.state_19951231 stateASCIIOutput
elif $stateOutputBinary ; then
    md5sum frs.state_19951231 stateBinaryOutput
fi

echo "Finished."
