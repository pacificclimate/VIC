#! /usr/bin/python
import re
import sys
import hashlib

f = open("metaData.txt", 'r')
metaDataContents = f.readlines()

dataFile = open("outputDefinitions.txt", 'r')
dataContents = dataFile.read()


for line in metaDataContents:
	tuples = re.findall(r'mapping\[\"(.*)\"\].*\"\".*\"\".*\"\".*\"\".*', line)
	if not tuples:
		sys.stdout.write(line) # no extra newline
	else:
		for item in tuples:
			tempName = "FIXME_"
			m = hashlib.md5()
			m.update(item)
			tempName += str(int(m.hexdigest(), 16))
			tempName = item		# ignore FIXME_md5hash, just put in the VIC internal variable name
			searchString = ".*" + item + "\s.*/\*\s*(.*?)\s*\[(.*?)\].*\s*\*/"
			definition = re.findall(searchString, dataContents)
			unitString = definition[0][1]
			if unitString == "fraction":
				unitString = "%"
			outString = "  mapping[\"" + item + "\"] = \t\tVariableMetaData(\"" + unitString + "\", \"" + tempName + "\", \"\", \"" + definition[0][0] + "\", \"\");\n"
			sys.stdout.write(outString)
