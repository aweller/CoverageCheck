# Copyright (C) 2014 Andreas M. Weller <andreas.m.weller@gmail.com>
#
# reads a file that maps GeneID <> Transcript ID <> RefSeq ID of the format
#
#Ensembl Gene ID	Ensembl Transcript ID	RefSeq mRNA [e.g. NM_001195597]
#ENSG00000261657	ENST00000566782	NM_001277303
#ENSG00000261657	ENST00000562780	
#
#
#and a file with exon information of the format
#
#Ensembl Gene ID	Ensembl Transcript ID	Exon Chr Start (bp)	Exon Chr End (bp)	Constitutive Exon	Exon Rank in Transcript	Chromosome Name	Associated Gene Name Strand
#ENSG00000252952	ENST00000517143	23791571	23791673	1	1	13	RNU6-58P
#ENSG00000207298	ENST00000384568	99677488	99677599	1	1	13	RNU6-83P
#
# finds the transcipt with the lowest RefSeq number per gene (i.e. the reference one used in diagnostics)
# returns output in BED format

import sys

id_mapping_file = sys.argv[1]
exon_info_file = sys.argv[2]

#########################################################################

gid_to_refseq = {}
accepted_transids = {}

# parse all transcipts for each gene
with open(id_mapping_file) as handle:
    for row in handle:
        if "Ensembl" in row or not "NM" in row: continue
        f = row.strip().split("\t")
        
        gid, transid, refseq = f
        refseq_int = int(refseq.split("_")[1])
        
        if not gid_to_refseq.get(gid):
            gid_to_refseq[gid] = []
            
        gid_to_refseq[gid].append([refseq_int, refseq, transid])

# decide on the oldest transcript per gene and save it 
for gid, refseqs in gid_to_refseq.iteritems():
    refseqs.sort(key=lambda x: x[0])
    accepted_transids[gid] = refseqs[0][2]
    
#########################################################################

strand_names = {"-1":"-", "1":"+"}

# parse the exon info and translate correct transcripts to bed
with open(exon_info_file) as handle:
    for row in handle:
        if "Ensembl" in row: continue
        f = row.strip().split("\t")
        
        gid, transid, start, stop, _, exon_count, chrom, gene, strand = f
        
        if accepted_transids.get(gid):
            if accepted_transids.get(gid) == transid:
                #print gene, transid
        
                output = ["chr"+chrom, start, stop, gene, exon_count, strand_names[strand]]
                print "\t".join([str(x) for x in output])
    