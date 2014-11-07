#!/bin/bash

netCDF=false;
stateOutputASCII=false
stateOutputBinary=false
stateOutputNetCDF=false

stateInputASCII=false
stateInputBinary=false
stateInputNetCDF=false

correctResultsDir="/home/james/code/hg/VIC/correctResults"
globalOptionsFile="glb_prb_base_BASIN_SCENARIO_19502006_VIC4.1.2_netcdf_auto.txt"
outputDir="out"

usage()
{
    echo "Unknown option $key"
    echo "Usage: testVIC.sh [--netCDF] [--state_output ASCII|binary] [--state_input ASCII|binary|NetCDF] [--global <global_options_file> [--output-dir <path>] [--output-comparison-dir <path>]";
    exit 1
}

echo $# $@

while [[ $# > 0 ]]
do
key="$1"
shift

echo $key $1

case $key in
    --netCDF)
        echo "netCDF option enabled"
        netCDF=true
    ;;
    --state_output)
        if [ "$1" == "ASCII"]; then
            echo "state output ASCII enabled"
            stateOutputASCII=true
        elif [ "$1" == "binary"]; then
            echo "state output Binary enabled"
            stateOutputBinary=true
        else
            echo "Incorrect --state_output value\nUsage: --state_output [ASCII|binary]"
            exit 1
        fi
        shift
    ;;
    --state_input)
        if [ "$1" == "ASCII"]; then
            echo "state input ASCII enabled"
            stateInputASCII=true
        elif [ "$1" == "binary"]; then
            echo "state input Binary enabled"
            stateInputBinary=true
        elif [ "$1" == "NetCDF"]; then
            echo "state input NetCDF enabled"
            stateInputNetCDF=true
        else
            echo "Incorrect --state_input value\nUsage: [--state_input ASCII|binary|NetCDF]"
            exit 1
        fi
        shift
    ;;
    --global)
        globalOptionsFile="$1"
        shift
    ;;
    --output-dir)
        outputDir="$1"
        shift
    ;;
    --output-comparison-dir)
        correctResultsDir="$1"
        shift
    ;;
    *)
        usage
    ;;
esac
done

echo starting test

codeDir="VIC_4.1.2_cpp_trunk"
programName="vicNl"
export curDir=$(pwd)
TIMESTAMP=$(date +"%Y_%m_%d__%H_%M_%S")
export outputName="automated_4.1.2_netcdf_$TIMESTAMP"
echo "Output for this test will be in $outputDir/$outputName"
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
mkdir $outputDir/$outputName

#Save the runtime config with the output
cp $globalOptionsFile $outputDir/$outputName

#Run the code
echo "running the program"
$codeDir/$programName -g $globalOptionsFile

#Compare output files with pristine version hashes
echo ""
if $netCDF ; then
    echo "validating NetCDF output results"
    Rscript vic_output_netcdf_compare.r $outputDir/$outputName/results.nc $correctResultsDir/vicNetCDFCorrectOutputs.nc
else
    echo "validating plain ASCII output results"
    pushd $outputDir/$outputName > /dev/null
    echo $(md5sum -c $correctResultsDir/timeSeriesChecksums.md5 2>/dev/null | grep -cE 'OK$') out of $(cat $correctResultsDir/timeSeriesChecksums.md5 | wc -l) files compare OK
    popd > /dev/null
fi

if $stateOutputASCII ; then
    md5sum frs.state_19951231 $correctResultsDir/stateASCIIOutput
elif $stateOutputBinary ; then
    md5sum frs.state_19951231 $correctResultsDir/stateBinaryOutput
fi

if $stateInputASCII ; then
    pushd $outputDir/$outputName > /dev/null
    md5sum -c $correctResultsDir/stateInputFluxesOutputASCII/checksum
    popd > /dev/null
elif $stateInputBinary ; then
    pushd $outputDir/$outputName > /dev/null
    md5sum -c $correctResultsDir/stateInputFluxesOutputBinary/checksum
    popd > /dev/null
elif $stateInputNetCDF ; then
    pushd $outputDir/$outputName > /dev/null
    md5sum -c $correctResultsDir/stateInputFluxesOutputBinary/checksum
    popd > /dev/null
fi

echo "Finished."

