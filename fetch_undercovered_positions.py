# read an output file from
# bedtools coverage -abam Q2PL2_H01_N.bam -b TSB_148_gene_panel_HP_amplicons.bed -d > test_coverage.tsv
# and find bases with coverage or strand_ratio below the threshold

import pandas as pd
import sys
import subprocess
import os
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

#TODO:
# deal with failed bedtool runs which leave empty files behind

#####################################################################################################################

class UndercoveredRegion():
    
    """
    Each instance corresponds to one contiguos run of bases 
    """
    
    def __init__(self, row):
        f = row.strip().split("\t")
        self.chrom, self.start, self.stop, self.amplicon, na, self.strand,  self.pos, self.value = f
        self.gene = self.amplicon.split("_")[0]
        
        self.first_start = int(self.start)
        self.first_pos = int(self.pos)
        
        self.last_start = int(self.start)
        self.last_stop = int(self.stop)
        self.last_pos = int(self.pos)
        
        self.rows = []
        self.values = []
        self.rows.append(row)
        self.values.append(self.value)
         
    def add_row(self,row):
        f = row.strip().split("\t")
        self.chrom, self.start, self.stop, self.amplicon, na, self.strand,  self.pos, self.value = f
        
        self.last_start = int(self.start)
        self.last_stop = int( self.stop)
        self.last_pos = int( self.pos)
        
        self.rows.append(row)
        self.values.append(self.value)
        
    def print_output(self):
        
        self.region_start = self.first_start + self.first_pos
        self.region_stop = self.last_start + self.last_pos
        
        mean_value = round(np.mean([float(x) for x in self.values]), 2)
        region_length = len(self.values)
        
        output = [self.gene, self.chrom, self.region_start, self.region_stop, self.strand, mean_value, region_length]
        
        return "\t".join([str(x) for x in output])

#####################################################################################################################

def find_undercovered_positions(coverage_file, dp_cutoff = None, whitelist_filename = None):
    
    region = None
    last_gene = None
    last_pos = 0
    last_gene = None
    
    total_bases = 0
    undercovered_bases = 0
        
    sample = coverage_file[:-13]
    undercovered_output_name = sample + "_undercovered_regions.txt"
    undercovered_output = open(undercovered_output_name, "w")
    
    output_header = ["gene", "chrom", "start", "stop", "strand", "mean_coverage", "size"]
    undercovered_output.write("\t".join(output_header) + "\n")
    
    if whitelist_filename:
        whitelist = {}
        with open(whitelist_filename) as handle:
            for row in handle:
                f = row.strip().split("\t")
                chrom, pos = f[:2]
                chrompos = chrom + "\t" + pos
                whitelist[chrompos] = row
    
    whitelist_output_name = sample + "_undercovered_regions.txt"
    whitelist_output = open(whitelist_output_name, "w")
    

    
    with open(coverage_file) as handle:
        for row in handle:
            f = row.strip().split("\t")
            chrom, start, stop, amplicon, na, strand,  pos, dp = f
            chrompos = chrom + "\t" + str(pos+start)
            
            start = int(start)
            stop = int(stop)
            pos = int(pos)
            dp = int(dp)
            gene = amplicon.split("_")[0]
            
            ##########################################################
            
            total_bases += 1
            if dp < dp_cutoff:
                undercovered_bases += 1
                
                if whitelist_filename:
                    if chrompos in whitelist:
                        whitelist_output.write( whitelist[chrompos] )
            
            ##########################################################
            
            if gene != last_gene or last_pos != (pos - 1):
                if region:
                    result = region.print_output()
                    undercovered_output.write(result + "\n")
                    region = None

            elif dp < dp_cutoff:
                if region:
                    region.add_row(row)
                else:
                    region = UndercoveredRegion(row)
            else:
                if region:
                    result = region.print_output()
                    undercovered_output.write(result + "\n")
                    region = None
            
            #print row.strip() + "\t" + str(region) 
            
            last_chrom = chrom
            last_pos = pos
            last_gene = gene
    
    undercovered_output.close()
    whitelist_output.close()
    
    ######################################################################
    
    covered_bases = total_bases - undercovered_bases
    pass_percent = 100 * round( covered_bases / float(total_bases), 3)
    
    print "%s percent (%s/%s) of positions in %s have at least the minimum coverage of %sX." % (pass_percent, covered_bases, total_bases, sample, dp_cutoff)

#####################################################################################################################

def find_strandbiased_positions(coverage_file, max_strand_ratio = None, whitelist_filename = None):
    
    filename = coverage_file
    header = ["chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]
    
    rawdf = pd.read_csv(filename, sep="\t", names=header)
    rawdf["pos"] = rawdf.start +rawdf.amplicon_pos
    rawdf["chrompos"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["pos"]) ]), axis = 1 ) 
    rawdf["gene"] = rawdf.apply(lambda x : x["amplicon"].split("_")[0], axis = 1 ) 
    
    # df: per base
    
    df = dict(chrom = rawdf.chr.groupby(rawdf.chrompos).min(), 
              pos = rawdf.pos.groupby(rawdf.chrompos).min(),
              gene = rawdf.gene.groupby(rawdf.chrompos).min(),
            start = rawdf.start.groupby(rawdf.chrompos).min(),
            stop = rawdf.stop.groupby(rawdf.chrompos).min(),
              minus_dp = rawdf[rawdf.strand == "-"].dp.groupby(rawdf.chrompos).max(),    
              plus_dp = rawdf[rawdf.strand == "+"].dp.groupby(rawdf.chrompos).max(),  
              dp = rawdf.dp.groupby(rawdf.chrompos).sum(),    
    )    
    df = pd.DataFrame(df).reset_index()
    
    def get_strand_ratios(row):
        minus = row["minus_dp"]
        plus = row["plus_dp"]
        
        if plus + minus == 0:
            return 0
        if plus == 0:
            return minus
        elif minus == 0:
            return plus
        else:
            larger = max([plus, minus])
            smaller = min([plus, minus])
            
            return abs(larger/float(smaller))
        
    df["strand_ratio"] = df.apply(get_strand_ratios, axis =1 )
    
    #################################################    
    
    region = None
    last_gene = None
    last_pos = 0
    last_gene = None
    
    total_bases = 0
    biased_bases = 0
    
    sample = coverage_file[:-13]
    strandbiased_output_name = sample + "_strandbiased_regions.txt"
    strandbiased_output = open(strandbiased_output_name, "w")
    
    output_header = ["gene", "chrom", "start", "stop", "strand", "mean_strand_bias", "size"]
    strandbiased_output.write("\t".join(output_header) + "\n")
    
    if whitelist_filename:
        whitelist = {}
        with open(whitelist_filename) as handle:
            for row in handle:
                f = row.strip().split("\t")
                chrom, pos = f[:2]
                chrompos = chrom + "\t" + pos
                whitelist[chrompos] = row
    
    whitelist_output_name = sample + "_strandbiased_regions.txt"
    whitelist_output = open(whitelist_output_name, "w")
    
    for index, pandas_dict in df.iterrows():
        chrompos, chrom, dp, gene, minus_dp, plus_dp, pos, start, stop, strand_ratio = pandas_dict
        row = [chrom, start, stop, gene, ".", "NA", (pos-start), strand_ratio]
        row = "\t".join([str(x) for x in row])
        
        start = int(start)
        stop = int(stop)
        pos = int(pos)
        
        ##########################################################
        
        total_bases += 1
        if strand_ratio > max_strand_ratio:
            biased_bases += 1
            
            if whitelist_filename:
                if chrompos in whitelist:
                    whitelist_output.write( whitelist[chrompos] )
        
        ##########################################################
        
        if gene != last_gene or last_pos != (pos - 1):
            if region:
                result = region.print_output()
                strandbiased_output.write(result + "\n")
                region = None

        elif strand_ratio > max_strand_ratio:
            if region:
                region.add_row(row)
            else:
                region = UndercoveredRegion(row)
        else:
            if region:
                result = region.print_output()
                strandbiased_output.write(result + "\n")
                region = None
                
        last_chrom = chrom
        last_pos = pos
        last_gene = gene
    
    strandbiased_output.close()
    whitelist_output.close()
    
    ######################################################################
    
    even_bases = total_bases - biased_bases
    pass_percent = 100 * round( even_bases / float(total_bases), 3)
    
    print "%s percent (%s/%s) of positions in %s have a strand bias below the threshold (%s:1)." % (pass_percent, even_bases, total_bases, sample, max_strand_ratio)
        
#####################################################################################################################

def run_bedtools(bam, output, bed = None):
    bedtools_cmd = "bedtools coverage -s -d -abam %s -b %s > %s" % (bam, bed, output)
    #print bedtools_cmd
    subprocess.call(bedtools_cmd, shell=True)

def check_bed_filetype(filename):
    """
    Check if the bedfile is actually a bed and convert it to one if its an Illumina Manifest file
    """
    
    rows = open(filename).readlines()
    
    if not "Header" in rows[0]: 
        
        if all([len(row.split("\t")) == 6 for row in rows]):
            print "Correct bed file detected:", filename
            return filename
        
        else:
            print "Error in input bed file:", filename
            print "Not all rows contain 6 fields as is expected."
            print "Sorry, aborting..."
            sys.exit()
    
    else: # this is an Illumina manifest file
        
        print "Input is not a bed file, but an Illumina manifest file:", filename
        
        expected_bed = filename[:-4] + "_plusminus.bed"
                
        if os.path.exists(expected_bed):
            print "Switching to existing bed file:", expected_bed
            return expected_bed
        
        else:
            print "Converting manifest to bed."
            
            import manifest2bed as m2b
            all_out, plusminus = m2b.convert_manifest(filename)
            
            print "Switching to newly created bed file:", plusminus
            return plusminus

#####################################################################################################################
#####################################################################################################################

def main():
    
    bed = sys.argv[1]
    bed = check_bed_filetype(bed)
    
    print "-" * 100
    
    dp_cutoff = 50
    max_strand_ratio = 5
    
    print "Minimum accepted coverage per base: %sX" % (dp_cutoff)
    print "Maximum accepted coverage ratio between strands: %s:1" % (max_strand_ratio)
    
    print "-" * 100
    print "-" * 100

    target_file = None
    
    if len(sys.argv) > 2:
        target_file = sys.argv[2]
    
    if target_file:
        find_undercovered_regions(target_file, graphs = True)
    else:
        bams = [x for x in os.listdir("./") if x.endswith(".bam")]    
        
        for bam in bams:
            bedtools_output = bam.split(".")[0] + "_coverage.tsv"
            
            if not os.path.exists(bedtools_output):
                print "Running bedtools coverage for:", bam
                run_bedtools(bam, output = bedtools_output, bed=bed)
            else:
                print "Coverage file %s already present, skipping calculation." % bedtools_output
            
            find_undercovered_positions(bedtools_output, dp_cutoff = dp_cutoff)
            find_strandbiased_positions(bedtools_output, max_strand_ratio = max_strand_ratio)
            
            print "Finished", bam
            print "-" * 100
    


##################################################################################

if __name__ == "__main__":
    main()

    
    