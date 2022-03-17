#!/usr/bin/env python
# Author: Jie Lu & Haktan Suren
# Date: 11/07/2016
# Modified: 05/29/2019
# Purpose: reformat TVC output file (.xls) to a matrix of marker by sample (barcode) where markers are the rows and the barcodes are the columns.

# Allowed File Formats
# HOTSPOT file: should be bed format (chr, start, end, marker name) , first row is the header. We only use marker name from hotspot bed file but it MUST be the 4th column.
# VariantCaller XLS file: We use the following column headers "Ref, Variant, Allele Call, Barcode, Allele Source, Allele Name, Sample Name". Column headers should match exactly to the names mentioned before. Order of the columns may change in the XLS file. Samples/Barcodes MUST be ordered per sample within the XLS file. 

# Output
# Genotype Calls are made like following: Homozygous: Var/Var, Absent: Ref/Ref, Heterozygous: Ref/Var, No Call: ./..

# Assumptions: 1) Only works for biallelic markers. 2) Multiallelic markers should be indicated as multiple rows in HOTSPOT file (assigning different name: e.g. M1_1 and M1_2). 3) Long indel markers are not reported in the output.  

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import common
import csv
from glob import glob
from itertools import izip
import json
import re

def convertToGenotype(ref, var, call):
	"""
	Converts genotype calls to allelic representation as x/x.
	This is the final format used in the output file
	"""
	if call == "No Call":
		return "./."
	elif call == "Homozygous":
		return "%s/%s" % (var, var)
	elif call == "Absent":
		return "%s/%s" % (ref, ref)
	elif call == "Heterozygous":
		return "%s/%s" % (ref, var)

path = sys.argv[1]
vid = sys.argv[2]

outputStyle = '2'

common.printtime('Path: '+path)
common.printtime('Output Style: '+outputStyle)


os.chdir(path)

# Getting the variantCaller ID, if it's 0, pick the latest run. IDs are incrementally created so the highest ID belongs to the latest run.
if vid == '0' :
	vids = glob('%s/../variantCaller_out.*' % path)
	vidsIDs = []
	for vidR in vids:
		vidsIDs.append(re.sub(".+\.(\d+)$","\\1",vidR))
	vidsIDs.sort()
	vid = vidsIDs.pop()
else :
	vid = sys.argv[2]

# Parse startplugin.json to get some constants about the run e.g. Run name, location of hotspot file etc.
varDir = os.path.abspath( '%s/../variantCaller_out.%s' % ( path, vid ))
startjson = varDir + "/startplugin.json"

with open(startjson) as jsonFile:
  tvc_jsonParams = json.load(jsonFile)

runName = tvc_jsonParams['expmeta'].get('run_name','')
xlsfile = '%s/%s.xls' % (varDir, runName)

common.printtime('Run Name: '+runName)
common.printtime('XLS File: '+xlsfile)

hotspotsFiles = []

hotspotFile = tvc_jsonParams['pluginconfig']['meta'].get('targetloci','')
if hotspotFile != '' and hotspotFile not in hotspotsFiles:
	common.printtime('HOTSPOT File From targetloci: '+hotspotFile)
	hotspotsFiles.append(hotspotFile)

if len(hotspotsFiles) == 0:
	for i in tvc_jsonParams['plan']['barcodedSamples']:
		barcodes = tvc_jsonParams['plan']['barcodedSamples'][i]['barcodes']
		for barcode in barcodes:
			hotspotFile = tvc_jsonParams['plan']['barcodedSamples'][i]['barcodeSampleInfo'][barcode]['hotSpotRegionBedFile']
			if hotspotFile not in hotspotsFiles:
				common.printtime('HOTSPOT File %s: %s' % (len(hotspotsFiles)+1, hotspotFile))
				hotspotsFiles.append(hotspotFile)

# hotspotfile = "/results/uploads/BED/6/ThreeAmigosCanineReference/unmerged/detail/AgriSeq_Canine_CP_Hotspot_mod1.bed" #fake incase we need to test the script

outputFile = '%s-%s-Matrix.xls' % (vid,runName)

common.printtime('Output File: '+outputFile)

"""
Create an array of marker names from hotspot file.
["Marker1","Marker2",...]
"""
hotspots = []
for hotspotfile in hotspotsFiles:
	with open(hotspotfile, 'r') as HOTSPOT:
		HOTSPOT.readline()  # assuming a header line in bed file
		for line in HOTSPOT:
			try:
				hotspot = line.strip().split("\t")[3] # assuming the hotspot name is the 4th column
				hotspots.append(hotspot)
			except IndexError: 
				raise Exception("Please check the hotspot.bed format")

TMPFILE = open("tmpout.tmp", 'w')
with open(xlsfile, 'r') as IN:
	line = IN.readline()
	header = line.strip().split("\t")
	# Get the position of the column using a column name match. So even the order of the column changes in future, the script can run without any interuptions. 
	for i in range(len(header)):
		if header[i] == "Ref":  # find the index of the GENOTYPE column
			ref_idx = i
		elif header[i] == "Variant":
			var_idx = i
		elif header[i] == "Allele Call":
			call_idx = i
		elif header[i] == "Barcode":
			barcode_idx = i
		elif header[i] == "Allele Source":
			allele_source_idx = i
		elif header[i] == "Allele Name": # find the index of the HOTSPOT column
			hotspot_idx = i
		elif header[i] == "Sample Name":
			sample_idx = i
			
	current_barcode = None
	tmpdict = {}
	
	# initialize the hotspot dictionary
	for i in hotspots:
		TMPFILE.write("\t" + i)  # write the hotspot name line
		tmpdict[i] = "NA"
	TMPFILE.write("\n")
	
	#This part could be written using dict (hash) to increase the readibility.
	#One of the biggest assumption with this part is XLS file should be sorted by Samples
	for line in IN:
		vect = line.strip().split("\t")
		if vect[allele_source_idx] == "Hotspot":
			barcode = vect[barcode_idx]
			ref = vect[ref_idx]
			var = vect[var_idx]
			call = vect[call_idx]
			hotspot = vect[hotspot_idx]
			hotspotlist = list(set(hotspot.split(",")))  #some hotspots are composite site separated by ",". This is due to duplicated markers in the hotspot file (same coordinates for multiple different names).
			
			if barcode == current_barcode:
				sample = vect[sample_idx]
				for hotspot in hotspotlist:
					if hotspot not in tmpdict:
						common.printtime("Unknown hotspot (not included in final report):" + hotspot)
					else:
						tmpdict[hotspot] = convertToGenotype(ref, var, call)
				
			else:                           # Assuming the barcodes are sorted, if reaching the next barcode
				if current_barcode != None:
					TMPFILE.write(sample+";"+current_barcode)
					for i in hotspots:
						TMPFILE.write("\t"+tmpdict[i])
					TMPFILE.write("\n")
					#print len(tmpdict)  for debug
				
				current_barcode = barcode
				tmpdict = {}
				for i in hotspots:
					tmpdict[i] = "NA"
					
				for hotspot in hotspotlist:
					if hotspot not in tmpdict:
						common.printtime("Unknown hotspot(not included in final report):" + hotspot)
					else:
						tmpdict[hotspot] = convertToGenotype(ref, var, call)

# write the last barcode
TMPFILE.write(sample+";"+current_barcode)
for i in hotspots:
	TMPFILE.write("\t" + tmpdict[i])
TMPFILE.write("\n")
			
HOTSPOT.close()
TMPFILE.close()

# Transpose the matrix
a = izip(*csv.reader(open("tmpout.tmp", "rb"), delimiter="\t"))
csv.writer(open("tmpout2.tmp", "wb"),delimiter="\t").writerows(a)

if outputStyle == "2":
	os.system("mv -f tmpout2.tmp %s" % outputFile)
	
# Clean up all temporary files
os.system("rm *.tmp")			
