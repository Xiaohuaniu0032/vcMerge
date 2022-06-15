from Bio import SeqIO
from Bio.Seq import Seq,translate
import sys
from csv import reader

def cds_region(feature_location):
    cdsregion = []
    for part in feature_location.parts:
        start = part.nofuzzy_start
        end = part.nofuzzy_end
        strand=part.strand
        cdsregion.append([start,end])
    return cdsregion

def pname(position,alt,start,gene_dict,gene,seq):
    pos = (position-start)/3
    AAcodon = (position-start)%3
    print AAcodon,'This is the codon postion'
    if AAcodon == 0:
        namepos = pos
        AApos = pos-1
        newseq= seq[position-3]+seq[position-2]+alt
    elif AAcodon == 1:
        namepos = pos+1
        AApos = pos
        newseq = alt+seq[position]+seq[position+1]
    elif AAcodon == 2:
        namepos = pos+1
        AApos = pos
        newseq = seq[position-2]+alt+seq[position]
    newAA = translate(newseq)
    print 'Amino Acid position: '+ str(AApos)
    if (position>=13468) and (position <21555): # orf1b gene, the AAposition should add orf1a GENE LENGTH +4401
        AApos += 4401
    oriAA = gene_dict[gene]['AA'][AApos]
    if newAA==oriAA:
        pname = ''.join([oriAA,str(namepos),'='])
    else:
        pname = ''.join([oriAA,str(namepos),str(newAA)])
    return pname

def delname(position,ref,start,gene_dict,gene,seq):
    pos_start = (position-start)/3
    AAcodon_start = (position-start)%3
    print AAcodon_start,'This is the codon postion'
    if AAcodon_start == 1: # del nuc is the first codon
        namepos = pos_start+1
        if (position>=13468) and (position <21555): # orf1b gene, the AAposition should add orf1a GENE LENGTH +4401
            AApos = pos_start + 4401
        else:
            AApos = pos_start
        delAA_list = []
        for i in range(len(ref)/3):
            delname = gene_dict[gene]['AA'][AApos+i]+str(namepos+i)+'-'
            delAA_list.append(delname)
    else: # del nuc break the AA codon
        if AAcodon_start == 2: # codon2
            newseq = seq[position-2]+seq[position+len(ref)-1]+seq[position+len(ref)]
            namepos = pos_start+1
            if (position>=13468) and (position <21555): # orf1b gene, the AAposition should add orf1a GENE LENGTH +4401
                AApos = pos_start + 4401
            else:
                AApos = pos_start
        elif AAcodon_start == 0: # codon3
            newseq= seq[position-3]+seq[position-2]+seq[position+len(ref)-1]
            namepos = pos_start
            if (position>=13468) and (position <21555): # orf1b gene, the AAposition should add orf1a GENE LENGTH +4401
                AApos = pos_start -1 + 4401
            else:
                AApos = pos_start -1
        oriAA = gene_dict[gene]['AA'][AApos]
        newAA = translate(newseq)
        delAA_list = []
        if oriAA != newAA:
            for i in range(len(ref)/3):
                delname = gene_dict[gene]['AA'][AApos+i]+str(namepos+i)+'-'
                delAA_list.append(delname)
        else:
            for i in range(len(ref)/3):
                delname = gene_dict[gene]['AA'][AApos+i+1]+str(namepos+i+1)+'-'
                delAA_list.append(delname)
    print 'Amino Acid position: '+ str(AApos)
    return ';'.join(delAA_list)

def run(gbfile,infile):
    gb_record = SeqIO.read(gbfile,'genbank')
    nucid = gb_record.id
    nucseq = gb_record.seq
    features = gb_record.features
    cds_features = filter(lambda i:i.type=='CDS',features)
    gene_dict = {}
    position_list = []
    for cds_feature in cds_features:
        qualifier = cds_feature.qualifiers
        location = cds_feature.location
        cds_position = cds_region(location)
        gene = qualifier['gene'][0]
        AA =  qualifier['translation'][0]
        gene_dict[gene]={}
        gene_dict[gene].update(AA=AA,location=cds_position)
    
    row = reader(open(infile), delimiter = '\t')
    header = next(row)
    newheader = header[:4]+['Gene','Amino Acid','Defining mutations']+header[4:]
    writeline = []
    writeline.append('\t'.join(newheader))
    for r in row:
        mutation = '\t'.join(r[1:4])
        print mutation
        ref = r[2]
        alt = r[3]
        position_start = int(r[1])
        position_end = position_start+len(alt)
        if (position_end<265) or (position_start>29674):
            newline = '\t'.join(r[:4]+['UTR','NA','NA']+r[4:])
            writeline.append(newline)
        else:
            Gene = None
            for gene in gene_dict:
                for exon in gene_dict[gene]['location']:
                    cds_start,cds_end = exon
                    if position_start>int(cds_start) and position_start <=int(cds_end):
                        Gene = gene
                        GeneName = gene.replace('orf','ORF')
                        start = cds_start
                        break
            if not Gene:
                newline = '\t'.join(r[:4]+['NA','NA','NA']+r[4:])
                writeline.append(newline)
            else:
                if ref == '-': # insertion, no example
                    codon_start = position_start-int(start)-1
                    codon_end = position_start+len(alt)-int(start)-1
                    AA_start = codon_start/3
                    AA_end = codon_end/3
                    AminoAcid = 'ins'+gene_dict[Gene]['AA'][AA_start]+str(AA_start+1)+'_'+gene_dict[Gene]['AA'][AA_end]+str(AA_end+1)
                    #AminoAcid = 'ins'+str(AA_start+1)+'_'+str(AA_end+1)
                    newline = '\t'.join(r[:4]+[GeneName,AminoAcid,'unknown']+r[4:])
                    writeline.append(newline)
                elif alt == '-': # deletion, no codon1,2,3
                    if len(ref) %3 != 0: # frameshit deletion
                        #AminoAcid = 'del'+str(AA_start+1)+'_'+str(AA_end+1)
                        DefineMutations = 'frameshift'
                    else:
                        DefineMutations = 'Nonsynonymous'
                    AminoAcid = delname(position_start,ref,int(start),gene_dict,Gene,nucseq)
                    newline = '\t'.join(r[:4]+[GeneName,AminoAcid,DefineMutations]+r[4:])
                    writeline.append(newline)
                else: # SNP
                    if len(ref) != len(alt): #frameshift
                        newline = '\t'.join(r[:4]+[GeneName,'unknown','unknown']+r[4:])
                        writeline.append(newline)
                    else:
                        AAcodon = (position_start-int(start))%3
                        if (AAcodon== 1) and len(alt) == 3: # both mutations in a 3-codon nuc
                            AApos = (position_start-int(start))/3
                            if (position_start>=13468) and (position_start <21555): # orf1b gene, the AAposition should add orf1a GENE LENGTH +4401
                                AApos += 4401
                            oriAA = gene_dict[Gene]['AA'][AApos]
                            newAA = translate(alt)
                            namepos = AApos+1
                            AminoAcid = ''.join([oriAA,str(namepos),str(newAA)])
                            DefineMutations = 'Nonsynonymous'
                            newline = '\t'.join(r[:4]+[GeneName,AminoAcid,DefineMutations]+r[4:])
                            writeline.append(newline)
                        else:
                            for i in range(len(alt)):
                                position = position_start+i
                                ref_allele = ref[i]
                                alt_allele = alt[i]
                                if ref_allele == alt_allele: # no mutation
                                    pass
                                else:
                                    print position,alt_allele,start,Gene,GeneName
                                    AminoAcid = pname(position,alt_allele,int(start),gene_dict,Gene,nucseq)
                                    if '=' in AminoAcid: # Synonymous mutation
                                        DefineMutations = 'Synonymous'
                                    else:
                                        DefineMutations = 'Nonsynonymous'
                                    newline = '\t'.join(['2019-nCoV',str(position),ref_allele,alt_allele]+[GeneName,AminoAcid,DefineMutations]+r[4:])
                                    print "New Line Info: ",position,ref_allele,alt_allele,GeneName,AminoAcid
                                    writeline.append(newline)
    outname = infile[:-3]+'annt.xls'
    with open(outname,'wb') as handle:
        handle.write('\n'.join(writeline))


def main():
    usage = "python ./nuc2pro.py 2019-NCov.gb SARS_CoV_2_variantCaller_out.17.TSVC_variants.merged.vcf.xls"
    pypath = sys.argv[0] # path of this script
    gbname = sys.argv[1] # first parameter: genbank file in the same path with this script "2019-NCov.gb"
    infile = sys.argv[2] # second parameter: mutation file
    run(gbname,infile)

if __name__=="__main__":
    main()
