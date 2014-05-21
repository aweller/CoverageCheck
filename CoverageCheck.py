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
import plot_exon_coverage as plotting
import plot_exon_coverage_all_samples as all_sample_plotting
import unite_coverage_files as UniteCoverage


#TODO:
# deal with failed bedtool runs which leave empty files behind

#####################################################################################################################

class BadRegion():
    
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
        
        output = [self.gene, self.chrom, self.region_start, self.region_stop, mean_value, region_length]
        
        return "\t".join([str(x) for x in output])

#####################################################################################################################

def parse_coverage_file_into_dataframe(coverage_file):
    """
    Parses the coverage file, transforms into a per-base matrix and returns a pandas DF
    """
    
    header = ["chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]    
    rawdf = pd.read_csv(coverage_file, sep="\t", names=header)
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
    df = df.sort(columns = ["chrom", "pos"])
    
    return df#####################################################################################################################

def find_bad_positions(coverage_matrix, target_folder = None, trait = None, samplename = None, trait_cutoff = None, whitelist_filename = None):    
    
    region = None
    last_gene = None
    last_pos = 0
    last_gene = None
    
    total_bases = 0
    bad_bases = 0
    
    sample = samplename
    bad_output_name = target_folder + sample + "_%s_regions.tsv" % (trait)
    bad_output = open(bad_output_name, "w")
    
    output_header = ["gene", "chrom", "start", "stop", "mean_strand_bias", "size"]
    bad_output.write("\t".join(output_header) + "\n")
    
    if whitelist_filename:
        whitelist = {}
        with open(whitelist_filename) as handle:
            for row in handle:
                if row[0] == "#": continue
                f = row.strip().split("\t")
                chrom, pos = f[:2]
                chrompos = chrom + "\t" + pos
                whitelist[chrompos] = row
    
        whitelist_output_name = target_folder + sample + "_variants_failed_%s.txt" % (trait)
        whitelist_output = open(whitelist_output_name, "w")
    
    for index, pandas_dict in coverage_matrix.iterrows():
        chrompos, chrom, dp, gene, minus_dp, plus_dp, pos, start, stop, strand_ratio = pandas_dict
        row = [chrom, start, stop, gene, ".", "NA", (pos-start), strand_ratio]
        row = "\t".join([str(x) for x in row])
        
        start = int(start)
        stop = int(stop)
        pos = int(pos)
        
        ##########################################################
        
        pass_check = False
        
        if trait == "coverage":
            if dp > trait_cutoff:
                pass_check = True
        
        elif trait == "strandbias":
            if strand_ratio < trait_cutoff:
                pass_check = True
        
        
        ##########################################################
        
        total_bases += 1
        if not pass_check:
            bad_bases += 1
            
            if whitelist_filename:
                if chrompos in whitelist:
                    whitelist_output.write( whitelist[chrompos] )
        
        ##########################################################
        
        if gene != last_gene or last_pos != (pos - 1):
            if region:
                result = region.print_output()
                bad_output.write(result + "\n")
                region = None

        elif not pass_check:
            if region:
                region.add_row(row)
            else:
                region = BadRegion(row)
        else:
            if region:
                result = region.print_output()
                bad_output.write(result + "\n")
                region = None
                
        last_chrom = chrom
        last_pos = pos
        last_gene = gene
    
    bad_output.close()
    
    if whitelist_filename:
        whitelist_output.close()
    
    ######################################################################
    
    good_bases = total_bases - bad_bases
    pass_percent = 100 * round( good_bases / float(total_bases), 3)
    
    if trait == "strandbias":
        print "%s percent (%s/%s) of positions in %s have a strand bias below the threshold (%s:1)." % (pass_percent, good_bases, total_bases, sample, trait_cutoff)
    elif trait == "coverage":
        print "%s percent (%s/%s) of positions in %s have at least the minimum coverage of %sX." % (pass_percent, good_bases, total_bases, sample, trait_cutoff)
    
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

def run(bed, target_folder, min_dp, max_strand_ratio, whitelist=None):
        
    print "-" * 100
    
    print "Minimum accepted coverage per base: %sX" % (min_dp)
    print "Maximum accepted coverage ratio between strands: %s:1" % (max_strand_ratio)
    
    target_folder = target_folder + "/" # just to be on the safe side    
    bams = [x for x in os.listdir(target_folder) if x.endswith(".bam")]    
    
    print "Target folder: %s" % target_folder
    if whitelist:
        print "Expected variants to check: %s" % whitelist
    else:
        print "No expected variants specified."
    
    print "Found %s bams to process." % len(bams)
    print "-" * 100
    
    all_sample_filename = "all_samples.tsv"
    UniteCoverage.unite(all_sample_filename, target_folder = target_folder)
    all_sample_plotting.plot_exon_coverage(all_sample_filename, target_folder = target_folder)
    
    for bam in bams:
        samplename = bam.split(".")[0]
        bedtools_output = target_folder + samplename + "_coverage.tsv"
        
        if not os.path.exists(bedtools_output):
            print "Running bedtools coverage for:", bam
            run_bedtools(target_folder + bam, output = bedtools_output, bed=bed)
        #else:
        #    print "Coverage file %s already present, skipping coverage calculation from bam." % bedtools_output
        
        coverage_matrix = parse_coverage_file_into_dataframe(bedtools_output)
        
        ##############################################################################################
        # run
        
        find_bad_positions(coverage_matrix, target_folder = target_folder, trait = "coverage",
                           samplename = samplename, trait_cutoff = min_dp, whitelist_filename=whitelist)
        find_bad_positions(coverage_matrix, target_folder = target_folder, trait = "strandbias",
                           samplename = samplename, trait_cutoff = max_strand_ratio, whitelist_filename=whitelist)
    
        ##############################################################################################
        # plot
        
        plotting.plot_exon_coverage(bedtools_output, target_folder = target_folder)
        
        print "-" * 100
    
##################################################################################

if __name__ == "__main__":
    bed = sys.argv[1]
    whitelist = sys.argv[2]
    
    if "/" in bed:
        target_folder = "/".join(bed.split("/")[:-1])
    else:
        target_folder = "."
    
    min_dp = 50
    max_strand_ratio = 5
    
    run(bed, target_folder, min_dp, max_strand_ratio, whitelist=whitelist)

    
    