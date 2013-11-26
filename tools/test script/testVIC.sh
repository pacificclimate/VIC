#!/bin/bash

netCDF=false;
stateOutputASCII=false
stateOutputBinary=false
stateOutputNetCDF=false

stateInputASCII=false
stateInputBinary=false
stateInputNetCDF=false

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
    if [ "$var" == "-stateOutputNetCDF" ]; then
        echo "state output NetCDF enabled";
        stateOutputNetCDF=true;
    fi
    if [ "$var" == "-stateInputBinary" ]; then
        echo "state input Binary enabled";
        stateInputBinary=true;
    fi
    if [ "$var" == "-stateInputASCII" ]; then
        echo "state input ASCII enabled";
        stateInputASCII=true;
    fi
    if [ "$var" == "-stateInputNetCDF" ]; then
        echo "state input NetCDF eneabled";
        stateInputNetCDF=true;
    fi
done

echo starting test

codeDir="VIC_4.1.2_cpp_trunk"
correctResultsDir="correctResults"
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
    perl -pi.bak -e 's/.*STATE_FORMAT.*$/STATE_FORMAT\tASCII/g' $globalOptionsFile
elif $stateOutputBinary ; then
    perl -pi.bak -e 's/.*STATENAME.*$/STATENAME\tfrs.state/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_FORMAT.*$/STATE_FORMAT\tBINARY/g' $globalOptionsFile
elif $stateOutputNetCDF ; then
    perl -pi.bak -e 's/.*STATENAME.*$/STATENAME\tfrs.state/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_FORMAT.*$/STATE_FORMAT\tNETCDF/g' $globalOptionsFile
else 
    perl -pi.bak -e 's/.*STATENAME.*$/#STATENAME\tfrs.state/g' $globalOptionsFile
fi

if $stateInputASCII ; then
    perl -pi.bak -e 's/.*INIT_STATE.*$/INIT_STATE\tstateASCIIOutput/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_FORMAT.*$/STATE_FORMAT\tASCII/g' $globalOptionsFile
elif $stateInputBinary ; then
    perl -pi.bak -e 's/.*INIT_STATE.*$/INIT_STATE\tstateBinaryOutput/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_FORMAT.*$/STATE_FORMAT\tBINARY/g' $globalOptionsFile
elif $stateInputNetCDF ; then
    perl -pi.bak -e 's/.*INIT_STATE.*$/INIT_STATE\tstateNetCDFOutput/g' $globalOptionsFile
    perl -pi.bak -e 's/.*STATE_FORMAT.*$/STATE_FORMAT\tNETCDF/g' $globalOptionsFile
else 
    perl -pi.bak -e 's/.*INIT_STATE.*$/#INIT_STATE\tinputState/g' $globalOptionsFile
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
    Rscript vic_output_netcdf_compare.r out/$outputName/results.nc $correctResultsDir/vicNetCDFCorrectOutputs.nc
else
    pushd out/$outputName > /dev/null
    echo $(md5sum -c ../../$correctResultsDir/timeSeriesChecksums.md5 2>/dev/null | grep -cE 'OK$') out of $(cat ../../$correctResultsDir/timeSeriesChecksums.md5 | wc -l) files compare OK
    popd > /dev/null
fi

if $stateOutputASCII ; then
    md5sum frs.state_19951231 stateASCIIOutput
elif $stateOutputBinary ; then
    md5sum frs.state_19951231 stateBinaryOutput
fi

if $stateInputASCII ; then
    pushd out/$outputName > /dev/null
    md5sum -c ../../correctResults/stateInputFluxesOutputASCII/checksum
    popd > /dev/null
elif $stateInputBinary ; then
    pushd out/$outputName > /dev/null
    md5sum -c ../../correctResults/stateInputFluxesOutputBinary/checksum
    popd > /dev/null
elif $stateInputNetCDF ; then
    pushd out/$outputName > /dev/null
    md5sum -c ../../correctResults/stateInputFluxesOutputBinary/checksum
    popd > /dev/null
fi

echo "Finished."

