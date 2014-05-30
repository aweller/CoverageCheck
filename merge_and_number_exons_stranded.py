# reads a file of the format
#
# chrom start           stop            strand  gene 
# 22	51227178	51227226	1	RPL23AP82
#
# and merges overlapping exons
# also adds a column with the no of the exon within the gene

class Region():
    
    """
    Each instance corresponds to one contiguos run of bases 
    """
    
    def __init__(self, row):
        f = row.strip().split("\t")
        self.chrom, self.start, self.stop, self.strand, self.gene = f
        
        self.first_start = int(self.start)
        
        self.last_start = int(self.start)
        self.last_stop = int(self.stop)

         
    def add_row(self,row):
        f = row.strip().split("\t")
        self.chrom, self.start, self.stop, self.strand, self.gene = f
        
        self.last_start = int(self.start)
        self.last_stop = int( self.stop)
        
    def print_output(self, exon_count):
        
        self.region_start = self.first_start 
        self.region_stop = self.stop
        
        output = [self.chrom, self.region_start, self.region_stop, self.gene, self.strand, exon_count]
        
        return output

####################################################################################################################

import sys

infile = sys.argv[1]

exon_count = 1
exon = None
last_stop = 0
last_gene = "NA"        

output_rows = []
strand_dict = {}
exon_count_per_gene = {"NA":1}

with open(infile) as handle:
    for row in handle:
        if "Exon" in row: continue
        
        f = row.strip().split("\t")
        chrom, start, stop, strand, gene = f

        start, stop, gene = int(f[1]), int(f[2]), f[4]
        
        if not strand_dict.get(gene):
            strand_dict[gene] = strand 
        
        if not exon:
            exon = Region(row)
        
        ##########################################################
        
        if gene != last_gene: # new gene, needs reset
            result = exon.print_output(exon_count_per_gene[last_gene])
            output_rows.append(result)
            exon_count_per_gene[gene] = 1
            exon = Region(row)
                
        elif start > last_stop: # exon has ended, but gene continues 
            result = exon.print_output(exon_count_per_gene[gene])
            output_rows.append(result)
            exon_count_per_gene[gene] += 1
            exon = Region(row)
            
        else: # exon continues
            exon.add_row(row)
        
        last_start = start
        last_stop = stop
        last_gene = gene        

####################################################################################################
# go through the output_rows and reverse the exon counts for the negative strands

strand_names = {"-1":"-", "1":"+"}

for row in output_rows:
    chrom, region_start, region_stop, gene, strand, exon_count = row
    count_range = range(1, exon_count_per_gene[gene]+1)
    
    if strand_dict[gene] == "1":    # plus strand, exon count is ok
        pass
    
    else: # minus strand, exon count needs to be reversed
        reverse_exon_count = count_range[-exon_count]
        row[-1] = reverse_exon_count 
    
    row[0] = "chr" + row[0]
    row[-2] = strand_names[strand]
    
    print "\t".join([str(x) for x in row])

        



