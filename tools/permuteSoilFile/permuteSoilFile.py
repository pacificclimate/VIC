#!/usr/bin/python

import re
import sys

def printUsage():
	print "Usage: inputFile outputFile options"
	print "options is either append [-append=number] or replace [-replace=index,number]"
	print "for example:"
	print "    input.txt output.txt -append=0.25"
	print "will read input.txt and append 0.25 to every line and write the result to output.txt"
	print "    input.txt output.txt -replace=0,10.57"
	print "will read input.txt and replace the first column (index 0) with the value 10.57 and write the result to output.txt"

def getArguments():
	if (len(sys.argv) != 4):
		printUsage()
		quit()
	return sys.argv[1:]	# remove program name


def convert():
	inputFilename, outputFilename, operation = getArguments()
	if (inputFilename == outputFilename):
		print "Error: input and output files cannot be the same!"
		printUsage()
		quit()
	inputFile = open(inputFilename, 'r')
	outputFile = open(outputFilename, 'w')
	
	if (operation.find("replace=") >=0):
		values = operation.split("=")[1].split(",")
		if (len(values) != 2):
			printUsage()
			quit()
		replace(inputFile, outputFile, values[0], values[1])
	elif (operation.find("append=")):
		value = operation.split("=")[1]
		append(inputFile, outputFile, value)

def replace(inputFile, outputFile, index, value):
	indexValue = int(index)
	for line in inputFile:
		numbers = line.split(" ")
		if (len(numbers) <= indexValue or indexValue < 0):
			print "Error: index (%s) must be in the range of the current number of items (%d) on the following line:" % (indexValue, len(numbers))
			print line
			print numbers
			quit()
		numbers[indexValue] = value
		modifiedLine = ' '.join(numbers)
		outputFile.write(modifiedLine)
		


def append(inputFile, outputFile, value):
	
	for line in inputFile:
		line = line.replace("\n", "")
		line = line + " " + str(value) + "\n"
		outputFile.write(line)
		


if __name__ == "__main__":
	convert()

