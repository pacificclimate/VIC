#! /usr/bin/python
import re
import sys

regexCellLine = re.compile(r'^([0-9]+)[ \t]+([0-9]+)')
regexVegLine = re.compile(r'^[ \t]+([0-9]+)[ \t]+([.0-9]+)(.*)')
insideVegLine = re.compile(r'[ \t]+([.0-9]+)')
snowBandGeneral = re.compile(r'^([0-9]+)(.*)')
snowBandData = re.compile(r'[ \t]+([.0-9]+)')


def loadSnowBandData(inputSnowBandFileName):
	with open(inputSnowBandFileName) as bandsFile:
		bands = dict()
		for line in bandsFile:
			tuples = snowBandGeneral.findall(line)
			if len(tuples) != 0:
				bands[tuples[0][0]] = tuples[0][1]

		if len(bands) <= 0:
			print "Error while parsing snow bands file (%s)" % inputSnowBandFileName
			quit(1)

		sampleBandData = bands.values()[0]
		bandData = snowBandData.findall(sampleBandData)
		numBands = len(bandData) / 3
		print "Number of bands calculated from snowband file is %d" % numBands
		return (bands, numBands)

def convertVegParam(inputVegFileName, inputSnowBandFileName, outputVegParamFileName):

	(bandData, numBands) = loadSnowBandData(inputSnowBandFileName)

	vegOutputFile = open(outputVegParamFileName, "w")

	with open(inputVegFileName) as vegInputFile:
		curCellNumber = -1
		curVegCount = 0
		maxVegInCell = 0
		curVegParamString = ""
		for line in vegInputFile:	# read line by line rather than putting the whole file in memory
			tuples = regexCellLine.findall(line)
			if len(tuples) != 0:
				# this is a cell number line
				curVegCount = 0
				curCellNumber = int(tuples[0][0])
				maxVegInCell = int(tuples[0][1])
				#print line + "found tuples" + str(tuples)
				cellVegString = "%d %d\n" % (curCellNumber, maxVegInCell * numBands)
				vegOutputFile.write(cellVegString)
			else:
				# this is a veg line under a cell number
				veggies = regexVegLine.findall(line)
				if len(veggies) != 0:
					#insides = insideVegLine.findall(str(veggies[0][2]))
					#print line + "found veg line" + str(veggies) + " insides " + str(insides)
					
					curBandValues = snowBandData.findall(bandData[str(curCellNumber)])
					for x in range(0, numBands):
						hruCV = float(veggies[0][1]) * float(curBandValues[x])
						# write out:	vegIndex Cv vegParams snowBandIndex
						outString = "\t%s\t%f\t%s\t%d\n" % (veggies[0][0], hruCV, veggies[0][2], x)
						vegOutputFile.write(outString)
					
					curVegCount = curVegCount + 1
			
			if curVegCount > maxVegInCell:
				print "Unexpected line! Expecting %d vegetation parameters but this is number %d! on line %s" % (maxVegInCell, curVegCount, line)
				print "quitting..."
				quit(1);

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Invalid arguments, usage:"
		print "./convertVegParam inputVegParamFile inputSnowBandFile outputVegParamFile"
		quit(1)
	convertVegParam(sys.argv[1], sys.argv[2], sys.argv[3])

