import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import common
import json
import re
import time

# Doesn't work for INDELs
# Doesn't work for multi-allelic markers
# all the absent calls are converted to "homozygous reference"
# -/- means absent for insertion or deletion
# ./. means no call 
# ATGC calls are shown for all non-SNP markers

def forwardToTopbottom(ABallele_list, genotypes):
	'''convert forward strand format to topbottom format given a list of genotypes'''
	ABallele_dict = {}
	ABallele_dict[ABallele_list[0]] = "A"
	ABallele_dict[ABallele_list[1]] = "B"
	topbottoms = []
	for genotype in genotypes:
		if genotype == 'NA':
			topbottoms.append("NA")
		else:
			allele1, allele2 = genotype.split("/")
			if allele1 == ".":
				topbottoms.append("./.")
			else:
				tb1 = ABallele_dict[allele1]
				tb2 = ABallele_dict[allele2]
				topbottom = tb1+tb2
				if topbottom == "BA":
					topbottom = "AB"
				topbottoms.append(topbottom)
	return topbottoms

def forwardToTop(REFallele_list, ABallele_list, genotypes):
	'''convert forward strand format to topbottom format given a list of genotypes'''
	ABallele_dict = {}
	ABallele_dict[REFallele_list[0]] = REFallele_list[0]
	ABallele_dict[REFallele_list[1]] = REFallele_list[1]
	if ABallele_list[0] == "BOT":
		ABallele_dict[REFallele_list[0]] = ABallele_list[1]
		ABallele_dict[REFallele_list[1]] = ABallele_list[2]
		
	topbottoms = []
	for genotype in genotypes:
		if genotype == 'NA':
			topbottoms.append("NA")
		else:
			allele1, allele2 = genotype.split("/")
			if allele1 == ".":
				topbottoms.append("./.")
			else:
				tb1 = ABallele_dict[allele1]
				tb2 = ABallele_dict[allele2]
				topbottom = tb1+tb2
				if topbottom == "GA":
					topbottom = "AG"
				if topbottom == "CA":
					topbottom = "AC"
				if topbottom == "TA":
					topbottom = "AT"
				if topbottom == "GC":
					topbottom = "CG"
				topbottoms.append(topbottom)
	return topbottoms

path = sys.argv[1]
vid = sys.argv[2]

os.chdir(path)

if vid == '0' :
	vids = glob('%s/../variantCaller_out.*' % path)
	vidsIDs = []
	for vidR in vids:
		vidsIDs.append(re.sub(".+\.(\d+)$","\\1",vidR))
	vidsIDs.sort()
	vid = vidsIDs.pop()
else :
	vid = sys.argv[2]

varDir = os.path.abspath( '%s/../variantCaller_out.%s' % ( path, vid ))
startjson = varDir + "/startplugin.json"

with open(startjson) as jsonFile:
  tvc_jsonParams = json.load(jsonFile)

runName = tvc_jsonParams['expmeta'].get('run_name','')
common.printtime('Run Name: '+runName)

matrixFile = '%s-%s-Matrix.xls' % (vid,runName)
common.printtime('matrixFile: '+matrixFile)

topBottomFile = '%s-%s.topbottom.txt' % (vid,runName)
common.printtime('topBottomFile: '+topBottomFile)

topFile = '%s-%s.top.txt' % (vid,runName)
common.printtime('topFile: '+topFile)

TOPBOTTOM = open(topBottomFile, 'r')
TOP = open(topFile, 'r')
MATRIX = open(matrixFile, 'r')
outbase = re.sub(r".xls", "", matrixFile)
OUT1 = open(outbase+".topbottom.xls", "w")
OUT2 = open(outbase+".top.xls", "w")
#timestamp = MATRIX.readline()
header = MATRIX.readline()
#OUT.write(timestamp)
OUT1.write(header)
OUT2.write(header)

ABalleles = {}
for line in TOPBOTTOM:
	line = line.rstrip().split("\t")
	id, alleles = line[3], line[9]
	ABallele = alleles.split(";")[1:]
	if id not in ABalleles:
		ABalleles[id] = ABallele

TOPBOTTOM.close()	

TOPalleles = {}
REFalleles = {}
for line in TOP:
	line = line.rstrip().split("\t")
	id, alleles = line[3], line[9]
	ref, alt = line[6].split(";")[0:2]
	ref = re.sub(r"REF=", "", ref)
	alt = re.sub(r"OBS=", "", alt)
	if id not in REFalleles:
		REFalleles[id] = [ref, alt]
	TOPallele = alleles.split(";")[0:]
	if id not in TOPalleles:
		TOPalleles[id] = TOPallele
TOP.close()	

for line in MATRIX:
	line = line.rstrip().split("\t")
	id = line[0]
	genotypes = line[1:]
	topgenotypes = line[1:]
	if id in ABalleles:
		ABallele_list = ABalleles[id]
		genotypes = forwardToTopbottom(ABallele_list, genotypes)
		OUT1.write(line[0]+"\t"+"\t".join(genotypes)+"\n")
	
	if id in TOPalleles:
		TOPallele_list = TOPalleles[id]
		REFallele_list = REFalleles[id]
		topgenotypes = forwardToTop(REFallele_list, TOPallele_list, topgenotypes)
		OUT2.write(line[0]+"\t"+"\t".join(topgenotypes)+"\n")

MATRIX.close()
OUT1.close()
OUT2.close()

