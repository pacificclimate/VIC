#!/usr/bin/python

import sys

def printUsage():
        print "Usage: inputFile outputFile"
	print "This converts a soil file (inputFile) into the new format (outputFile) by adding default values for the new parameters"

def getArguments():
        if (len(sys.argv) != 3):
                printUsage()
                quit()
        return sys.argv[1:]     # remove program name

# Order of new variables in the soil file is as follows (as defined by read_soilparam.c)
#NEW_SNOW_ALB = 0.85;    
#SNOW_ALB_ACCUM_A = 0.94;
#SNOW_ALB_ACCUM_B = 0.58;
#SNOW_ALB_THAW_A = 0.82; 
#SNOW_ALB_THAW_B = 0.46; 
#MIN_RAIN_TEMP = 0.0;   
#MAX_SNOW_TEMP = 6.0;    
#PADJ = 1.0;             
#T_LAPSE = 6.5;          
#PGRAD = 1.0;            
#AREA = 0;               
#GLAC_SURF_THICK = 100.0;
#GLAC_SURF_WE = 91.7;    
#GLAC_KMIN = 0.05;       
#GLAC_DK = 0.75;         
#GLAC_A = 0.01;          
#GLAC_ALBEDO = 0.3;      
#GLAC_ROUGH = 0.002;     

valuesToAdd = ["0.85", "0.94", "0.58", "0.82", "0.46", "0.0", "6.0", "1.0", "6.5", "1.0", "0.0", "100.0", "91.7", "0.05", "0.75", "0.01", "0.3", "0.002"]

def addAllValues():
	inputFilename, outputFilename = getArguments()
	if (inputFilename == outputFilename):
		print "Error: input and output files cannot be the same!"
		printUsage()
		quit()
	inputFile = open(inputFilename, 'r')
	outputFile = open(outputFilename, 'w')
	
	stringToAppend = ""
	for i in range(len(valuesToAdd)):
		stringToAppend = stringToAppend + " " + valuesToAdd[i];

	print "appending %s to all lines..."%stringToAppend

	for line in inputFile:
		line = line.replace("\n", "")
		line = line + stringToAppend + "\n"
		outputFile.write(line)


if __name__ == "__main__":
	addAllValues()

