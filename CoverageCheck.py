# read an output file from
# bedtools coverage -abam Q2PL2_H01_N.bam -b TSB_148_gene_panel_HP_amplicons.bed -d > test_coverage.tsv
# and find bases with coverage or strand_ratio below the threshold

# general modules
import pandas as pd
import sys
import subprocess
import os
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import pprint

# personal modules from the same folder
import plot_exon_coverage as plotting
import plot_exon_coverage_all_samples as all_sample_plotting
import unite_coverage_files as UniteCoverage
import CoverageCheckClasses

#####################################################################################################################

def parse_coverage_file_into_dataframe(coverage_file):
    """
    Parses the coverage file, transforms into a per-base matrix and returns a pandas DF
    """
    
    header = ["chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]    
    rawdf = pd.read_csv(coverage_file, sep="\t", names=header)
    rawdf["pos"] = rawdf.start +rawdf.amplicon_pos
    rawdf["chrompos"] = rawdf.apply(lambda x : "\t".join([str(x["chr"]), str(x["pos"]) ]), axis = 1 ) 
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
    
    return df
    
#####################################################################################################################

def find_bad_positions(coverage_matrix, target_folder = None, trait = None, samplename = None,
                       trait_cutoff = None, whitelist = None):    
    
    region = None
    last_gene = None
    last_pos = 0
    last_gene = None
    
    total_bases = 0
    bad_bases = 0
    
    sample = samplename
    bad_output_name = target_folder + sample + "_failed_regions_%s.tsv" % (trait)
    bad_output = open(bad_output_name, "w")
    
    output_header = ["gene", "chrom", "start", "stop", "mean_strand_bias", "size"]
    bad_output.write("\t".join(output_header) + "\n")

    for index, pandas_dict in coverage_matrix.iterrows():
        chrompos, chrom, dp, gene, minus_dp, plus_dp, pos, start, stop, strand_ratio = pandas_dict
        row = [chrom, start, stop, gene, ".", "NA", (pos-start), strand_ratio]
        row = "\t".join([str(x) for x in row])
        
        start = int(start)
        stop = int(stop)
        pos = int(pos)
        
        if whitelist:
            if whitelist.dict.get(chrompos):
                whitelist.add_coverage(chrompos, dp)   
                whitelist.add_strand_ratio(chrompos, strand_ratio)   
        
        ##########################################################
        
        pass_check = False
        
        if trait == "coverage":
            if dp > trait_cutoff:
                pass_check = True
        
        elif trait == "strandbias":
            if strand_ratio < trait_cutoff:
                pass_check = True
        
        total_bases += 1
        if not pass_check:
            bad_bases += 1
        
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
                region = CoverageCheckClasses.BadRegion(row)
        else:
            if region:
                result = region.print_output()
                bad_output.write(result + "\n")
                region = None
                
        last_chrom = chrom
        last_pos = pos
        last_gene = gene
    
    bad_output.close()
    
    ######################################################################
    
    good_bases = total_bases - bad_bases
    pass_percent = 100 * round( good_bases / float(total_bases), 3)
    
    if trait == "strandbias":
        print "%s percent (%s/%s) of positions in %s have a strand bias below the threshold (%s:1)." % (pass_percent, good_bases, total_bases, sample, trait_cutoff)
    elif trait == "coverage":
        print "%s percent (%s/%s) of positions in %s have at least the minimum coverage of %sX." % (pass_percent, good_bases, total_bases, sample, trait_cutoff)
    
    return whitelist
    
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

def remove_empty_files_from_folder(folder):
    for filename in os.listdir(folder):
        if os.path.getsize(folder + filename) == 0:
            os.remove(folder + filename)

#####################################################################################################################
#####################################################################################################################

def run(bed, target_folder, min_dp, max_strand_ratio, whitelist_filename=None):
    
    remove_empty_files_from_folder(target_folder) # remove empty files that might have been left over from previous runs
    bed = check_bed_filetype(bed)
        
    print "-" * 100
    
    print "Minimum accepted coverage per base: %sX" % (min_dp)
    print "Maximum accepted coverage ratio between strands: %s:1" % (max_strand_ratio)
    
    target_folder = target_folder + "/" # just to be on the safe side    
    bams = [x for x in os.listdir(target_folder) if x.endswith(".bam")]    
    
    print "Target folder: %s" % target_folder
    if whitelist_filename:
        variant_no = len([x for x in open(whitelist_filename).readlines() if x[0] != "#"])
        print "Found %s expected variants to check in %s" % (variant_no, whitelist_filename)
        
        if variant_no == 0:
            print "ERROR: expected variant file doesn't contain variants."
            sys.exit()
        
    else:
        print "No expected variants specified."
    
    print "Found %s bams to process." % len(bams)
    print "-" * 100

    ##############################################################################################
    # process each samples separately
    
    #for bam in bams:
    #    
    #    samplename = bam.split(".")[0]
    #    sample_output_folder = "%s/%s_results/" % (target_folder, samplename)
    #    
    #    if not os.path.exists(sample_output_folder):
    #        os.makedirs(sample_output_folder)
    #    
    #    ##############################################################################################
    #    # create the coverage file (if necessary) and parse it into a pandas DF
    #    
    #    bedtools_output = target_folder + samplename + "_coverage.tsv"
    #    
    #    if not os.path.exists(bedtools_output):
    #        print "Running bedtools coverage for:", bam
    #        run_bedtools(target_folder + bam, output = bedtools_output, bed=bed)
    #    
    #    coverage_matrix = parse_coverage_file_into_dataframe(bedtools_output)
    #    
    #    ##############################################################################################
    #    # Initialize the class instance thats collects information on each expected variant
    #    
    #    whitelist = None
    #    if whitelist_filename:
    #        whitelist = CoverageCheckClasses.ExpectedVariants(whitelist_filename, samplename=samplename, folder=sample_output_folder,
    #                                    dp_cutoff = min_dp, strandbias = max_strand_ratio)
    #    
    #    ##############################################################################################
    #    # run
    #    
    #    whitelist = find_bad_positions(coverage_matrix, target_folder = sample_output_folder, trait = "coverage",
    #                       samplename = samplename, trait_cutoff = min_dp, whitelist=whitelist)
    #    whitelist = find_bad_positions(coverage_matrix, target_folder = sample_output_folder, trait = "strandbias",
    #                       samplename = samplename, trait_cutoff = max_strand_ratio, whitelist=whitelist)
    #    
    #    if whitelist_filename:
    #        whitelist.print_output()
    #    
    #    ##############################################################################################
    #    # plot
    #    
    #    plotting.plot_exon_coverage(bedtools_output, target_folder = sample_output_folder)
    #    print "-" * 100
    #    
    ##############################################################################################
    # create summaries across all samples
    # this has to happen after the individual bam treatment because it assumes that coverage files for each bam are present
    
    all_sample_filename = "all_samples.tsv"
    UniteCoverage.unite(all_sample_filename, target_folder = target_folder)
    
    byte_size = os.path.getsize(target_folder+all_sample_filename)
    mb_size = byte_size/1.049e+6
    
    if mb_size > 1000:
        print "The united coverage size of all samples is %sMB. Skipping summary plots." % (round(mb_size,2))    
    else:
        all_sample_plotting.plot_exon_coverage(all_sample_filename, target_folder = target_folder)
    
    all_sample_plotting.plot_exon_coverage(all_sample_filename, target_folder = target_folder)
    
##################################################################################

if __name__ == "__main__":
    bed = sys.argv[1]
    whitelist_filename = sys.argv[2]
    
    if whitelist_filename == "None":
        whitelist_filename = None
    
    if "/" in bed:
        target_folder = "/".join(bed.split("/")[:-1]) + "/"
    else:
        target_folder = "./"
    
    min_dp = 50
    max_strand_ratio = 5
    
    run(bed, target_folder, min_dp, max_strand_ratio, whitelist_filename=whitelist_filename)

    
    