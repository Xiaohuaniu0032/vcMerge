import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import common
import json
import re
import time

# Step 1: Read the hotspot file and get the position information of the markers
# Step 2: Extract the marker flanking sequences from the reference genome - creates an intermediate file (.seq)
# Step 3: Read intermediate file to convert the ref / alt alleles TOP/BOTTOM alleles

def snp_to_top_bot(snp):
    strand_determination = {'A/T': '?', 'T/A': '?', 'G/C': '?', 'C/G': '?', 'A/C': 'TOP', 'A/G': 'TOP', 'C/A': 'TOP',
                            'G/A': 'TOP', 'T/C': 'BOT', 'T/G': 'BOT', 'C/T': 'BOT', 'G/T': 'BOT'}
    pseudo_strand_determination = {'A/T': '?', 'T/A': '?', 'G/C': '?', 'C/G': '?', 'A/C': 'TOP', 'A/G': 'TOP',
                               'C/A': 'TOP', 'G/A': 'TOP', 'T/C': 'BOT', 'T/G': 'BOT', 'C/T': 'BOT', 'G/T': 'BOT',
                               'A/A': '?', 'T/T': '?', 'C/C': '?', 'G/G': '?'}
    normal_allele_designation = {'A/C': ['A', 'C'], 'A/G': ['A', 'G'], 'C/A': ['A', 'C'], 'G/A': ['A', 'G'],
                             'T/C': ['T', 'C'], 'T/G': ['T', 'G'], 'C/T': ['T', 'C'], 'G/T': ['T', 'G']}
    special_allele_top_designation = {'A/T': ['A', 'T'], 'T/A': ['A', 'T'], 'C/G': ['C', 'G'], 'G/C': ['C', 'G']}
    special_allele_bot_designation = {'A/T': ['T', 'A'], 'T/A': ['T', 'A'], 'C/G': ['G', 'C'], 'G/C': ['G', 'C']}

    regex = re.compile('([ATGCN]+)\[([ATGC]/[ATGC])\]([ATGCN]+)')
    re_match = re.match(regex, snp)
    if re_match:
        left_flank = re_match.group(1)[::-1]
        genotype = re_match.group(2)
        right_flank = re_match.group(3)
        strand = strand_determination[genotype]
        if strand == 'TOP' or strand == 'BOT':
            a_allele = normal_allele_designation[genotype][0]
            b_allele = normal_allele_designation[genotype][1]
            return {'strand': strand, 'A allele': a_allele, 'B allele': b_allele}
        else:
            for l_base, r_base in zip(left_flank, right_flank):
                pseudo_genotype = l_base + '/' + r_base
                pseudo_strand = pseudo_strand_determination[pseudo_genotype]
                if pseudo_strand != '?':
                    if l_base == 'A' or l_base == 'T':
                        strand = 'TOP'
                        a_allele = special_allele_top_designation[genotype][0]
                        b_allele = special_allele_top_designation[genotype][1]
                        return {'strand': strand, 'A allele': a_allele, 'B allele': b_allele}
                    else:
                        strand = 'BOT'
                        a_allele = special_allele_bot_designation[genotype][0]
                        b_allele = special_allele_bot_designation[genotype][1]
                        return {'strand': strand, 'A allele': a_allele, 'B allele': b_allele}
                    break
            else:
                common.printtime('ERROR: Strand walking procedure terminated with no possible strand determination!')
    else:
        common.printtime('ERROR: {} is not of the format \'([ATGC]+)\[([ATGC]/[ATGC])\]([ATGC]+)\''.format(snp))

    return


def snp_to_top(snp):
    strand_determination = {'A/T': '?', 'T/A': '?', 'G/C': '?', 'C/G': '?', 'A/C': 'TOP', 'A/G': 'TOP', 'C/A': 'TOP',
                            'G/A': 'TOP', 'T/C': 'BOT', 'T/G': 'BOT', 'C/T': 'BOT', 'G/T': 'BOT'}
    pseudo_strand_determination = {'A/T': '?', 'T/A': '?', 'G/C': '?', 'C/G': '?', 'A/C': 'TOP', 'A/G': 'TOP',
                               'C/A': 'TOP', 'G/A': 'TOP', 'T/C': 'BOT', 'T/G': 'BOT', 'C/T': 'BOT', 'G/T': 'BOT',
                               'A/A': '?', 'T/T': '?', 'C/C': '?', 'G/G': '?'}
    normal_allele_designation = {'A/C': ['A', 'C'], 'A/G': ['A', 'G'], 'C/A': ['C', 'A'], 'G/A': ['G', 'A'],
                             'T/C': ['A', 'G'], 'T/G': ['A', 'C'], 'C/T': ['G', 'A'], 'G/T': ['C', 'A']}
    special_allele_top_designation = {'A/T': ['A', 'T'], 'T/A': ['T', 'A'], 'C/G': ['C', 'G'], 'G/C': ['G', 'C']}
    special_allele_bot_designation = {'A/T': ['T', 'A'], 'T/A': ['A', 'T'], 'C/G': ['G', 'C'], 'G/C': ['C', 'G']}

    regex = re.compile('([NATGC]+)\[([ATGC]/[ATGC])\]([NATGC]+)')
    re_match = re.match(regex, snp)
    if re_match:
        left_flank = re_match.group(1)[::-1]
        genotype = re_match.group(2) # Genotype alleles without brackets from hotspot (example: A/G)
        right_flank = re_match.group(3)
        strand = strand_determination[genotype] # Strand determination
        if strand == 'TOP' or strand == 'BOT': # Use unambiguos alleles only (A/G) 
            a_allele = normal_allele_designation[genotype][0]
            b_allele = normal_allele_designation[genotype][1]
            return {'strand': strand, 'A allele': a_allele, 'B allele': b_allele}
        else: # For unambiguos (A/T or G/C) walk on both sides
            for l_base, r_base in zip(left_flank, right_flank):
                pseudo_genotype = l_base + '/' + r_base
                pseudo_strand = pseudo_strand_determination[pseudo_genotype]
                if pseudo_strand != '?': # Ignoring A/A T/T C/C G/G while walking
                    if l_base == 'A' or l_base == 'T':
                        strand = 'TOP'
                        a_allele = special_allele_top_designation[genotype][0]
                        b_allele = special_allele_top_designation[genotype][1]
                        return {'strand': strand, 'A allele': a_allele, 'B allele': b_allele}
                    else:
                        strand = 'BOT'
                        a_allele = special_allele_bot_designation[genotype][0]
                        b_allele = special_allele_bot_designation[genotype][1]
                        return {'strand': strand, 'A allele': a_allele, 'B allele': b_allele}
                    break
            else:
                common.printtime('ERROR: Strand walking procedure terminated with no possible strand determination!')
    else:
        common.printtime('ERROR: {} is not of the format \'([NATGC]+)\[([ATGC]/[ATGC])\]([NATGC]+)\''.format(snp))

    return


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

hotspotsFiles = {}

hotspotFile = tvc_jsonParams['pluginconfig']['meta'].get('targetloci','')
referenceFile = tvc_jsonParams['pluginconfig']['meta'].get('reference','')
if hotspotFile != '' and referenceFile != '' and referenceFile not in hotspotsFiles:
	common.printtime('Hotspot file From targetloci: '+hotspotFile)
	common.printtime('Reference File From meta: %s' % referenceFile)
	hotspotsFiles[referenceFile] = hotspotFile

if len(hotspotsFiles) == 0:
	for i in tvc_jsonParams['plan']['barcodedSamples']:
		barcodes = tvc_jsonParams['plan']['barcodedSamples'][i]['barcodes']
		for barcode in barcodes:
			hotspotFile = tvc_jsonParams['plan']['barcodedSamples'][i]['barcodeSampleInfo'][barcode]['hotSpotRegionBedFile']
			referenceFile = tvc_jsonParams['plan']['barcodedSamples'][i]['barcodeSampleInfo'][barcode]['reference']
			if referenceFile not in hotspotsFiles:
				common.printtime('HOTSPOT File %s: %s' % (len(hotspotsFiles)+1, hotspotFile))
				common.printtime('Reference File %s: %s' % (len(hotspotsFiles)+1, referenceFile))
				hotspotsFiles[referenceFile] = hotspotFile

runName = tvc_jsonParams['expmeta'].get('run_name','')
common.printtime('Run Name: '+runName)

outputFile = '%s-%s' % (vid,runName)
common.printtime('outputFile: '+outputFile)

flanking = 100
OUT = open(outputFile+".seq", 'w')

for referenceFile in hotspotsFiles:
	
	hotspotfile = hotspotsFiles[referenceFile]
	referencePath = "/results/referenceLibrary/tmap-f3/%s/%s.fasta" % (referenceFile,referenceFile)
	referenceIndexPath = "%s.fai" % referencePath

	chrSize = {}
	RIN = open(referenceIndexPath, 'r') # hotspot file
	for line in RIN:
		line = line.rstrip().split("\t")
		chrSize[str(line[0])] = int(line[1])
	RIN.close()

	IN = open(hotspotfile, 'r') # hotspot file
	IN.readline()
	#CHR1    3249056 3249057 ARS-USMARC-Parent-DQ381153-rs29012842   0       +       REF=G;OBS=T;ANCHOR=A    SP_78.2684
	for line in IN:
		all_line = line
		line = line.rstrip().split("\t")
	
		chr, start, end, id = line[0:4]
		ref, alt = line[6].split(";")[0:2]
		ref = re.sub(r"REF=", "", ref)
		alt = re.sub(r"OBS=", "", alt)
		var_len = int(end)-int(start)

		if var_len == 1 and len(ref) == 1 and len(alt) == 1: #works for SNPs only
			if os.path.exists('tmp.bed'):
				os.system("rm tmp.bed")
			BED = open("tmp.bed", 'w')	
			bed_start_5 = max(int(start) - flanking + 1, 0)
			bed_end_5 = int(start)
			
			bed_start_3 = int(end)  
			bed_end_3 = int(end) + flanking - 1
			
			bed_start_5 = max(1, bed_start_5)
			bed_end_5 = min(chrSize[chr] if chr in chrSize else bed_end_5, bed_end_5)
			
			bed_start_3 = max(1, bed_start_3)
			bed_end_3 = min(chrSize[chr] if chr in chrSize else bed_end_3, bed_end_3)
			
			BED.write(chr+"\t"+str(bed_start_5) + "\t" + str(bed_end_5)+"\n")
			BED.write(chr+"\t"+str(bed_start_3) + "\t" + str(bed_end_3)+"\n")
			BED.close()
						
			try:
				os.system("bedtools getfasta -fi %s -bed tmp.bed -fo tmp.bed.fasta -tab" % referencePath)
			except:
				common.printtime("Error: no sequence can be retrieved for %s (%s:%s-%s)" % (id, chr, start, end))
				continue
			FASTA = open("tmp.bed.fasta", 'r')
			seq = []
			for fasta in FASTA:
				fasta = fasta.rstrip().split("\t")
				fasta[1].replace("BDHMRSUVWXYK","N") # Replacing the IUPAC characters with N
				seq.append(fasta[1])
			FASTA.close()
			#system("rm tmp.bed.fasta")
			final_seq = "%s[%s/%s]%s" % (seq[0], ref, alt, seq[1])
			OUT.write(("\t").join(line[0:9])+"\t"+final_seq+"\n")	
		else:
			common.printtime("Warning: %s (%s:%s-%s) is not a SNP" % (id, chr, start, end))
	IN.close()
OUT.close()
	
IN = open(outputFile+".seq", 'r')
OUT1 = open(outputFile+".topbottom.txt", 'w')
OUT2 = open(outputFile+".top.txt", 'w')
for line in IN:
	all_line = line.rstrip()
	line = line.rstrip().split("\t")
	seq = line[8]
	topbottom = snp_to_top_bot(seq)
	top = snp_to_top(seq)
	if topbottom != "None":
		OUT1.write(all_line+"\t"+topbottom["strand"]+";"+topbottom["A allele"] + ";" + topbottom["B allele"] + "\n")
		OUT2.write(all_line+"\t"+top["strand"]+";"+top["A allele"] + ";" + top["B allele"] + "\n")
	else:
		OUT1.write(all_line + "\tNA\n")
		OUT2.write(all_line + "\tNA\n")
IN.close()
OUT1.close()
OUT2.close()