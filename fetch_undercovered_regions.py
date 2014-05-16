# read an output file from
# bedtools coverage -abam Q2PL2_H01_N.bam -b TSB_148_gene_panel_HP_amplicons.bed -d > test_coverage.tsv
# and find regions with a set % of bases below the threshold

import pandas as pd
import sys
import subprocess
import os
import matplotlib.pylab as plt
import seaborn as sns

#TODO:
# deal with failed bedtool runs which leave empty files behind

def find_undercovered_regions(bam, dp_cutoff = None, pass_cutoff = None, strand_diff_cutoff = None, output = None, graphs = False):
    
    if not dp_cutoff:
        dp_cutoff = 50
    if not pass_cutoff:
        pass_cutoff = 0.5
    if not strand_diff_cutoff:
        strand_diff_cutoff = 0.5
    
    ##################################################################
    # read in the bedtools output as rawdf
    header = ["chr", "start", "stop", "amplicon", "na", "strand",  "pos", "dp"]
    rawdf = pd.read_csv(bam, sep="\t", names=header)
    rawdf["above_cutoff"] = rawdf.dp > dp_cutoff
    rawdf["name"] = rawdf.amplicon 
    rawdf["amplicon"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["start"]), str(x["stop"]) ]), axis = 1 )
    
    # check for empty bams
    covered_amplicons = len( rawdf[rawdf.dp > 0].amplicon.groupby(rawdf.amplicon) )
    total_amplicons = len( rawdf.amplicon.groupby(rawdf.amplicon) )
    covered_ratio = round( (covered_amplicons/float(total_amplicons)) * 100, 2)
    
    #print "Only %s percent (%s/%s) of amplicons in %s are covered at all." % (covered_ratio, covered_amplicons, total_amplicons, bam)
     
    if covered_ratio < 10:
        print "Only %s percent (%s/%s) of amplicons in %s have any reads mapping to them." % (covered_ratio, covered_amplicons, int(total_amplicons), bam)
        print "This looks like a bad bam or the wrong manifest file."
        print "Skipping %s..." % bam 
        return None
    else:
        print "%s percent (%s/%s) of amplicons in %s have reads mapping to them. Proceeding..." % (covered_ratio, covered_amplicons, int(total_amplicons), bam)

    
    ##################################################################
    # transform to per-amplicon/per-strand df: adf
    adf_dict= dict(name = rawdf.name.groupby([rawdf.amplicon, rawdf.strand]).max(),
                size = rawdf.pos.groupby([rawdf.amplicon, rawdf.strand]).max(),
               total_pass = rawdf[rawdf.above_cutoff].dp.groupby([rawdf.amplicon, rawdf.strand]).count(),
               median_dp = rawdf.dp.groupby([rawdf.amplicon, rawdf.strand]).median(),
               )
    adf = pd.DataFrame(adf_dict) 
    adf = adf.reset_index()
    adf = adf.fillna(0)    
    adf["pass_ratio"] = adf.total_pass / adf.size.astype(float)
    
    ##################################################################
    # summarize in the final per-amplicon df 

    df_dict= dict(name = adf.name.groupby(adf.amplicon).max(),
                    size = adf.size.groupby(adf.amplicon).max(),
                    total_pass = adf.total_pass.groupby(adf.amplicon).sum(),
                   minus_pass_ratio = adf[adf.strand == "-"].pass_ratio.groupby(adf.amplicon).max(),
                    plus_pass_ratio = adf[adf.strand == "+"].pass_ratio.groupby(adf.amplicon).max(),
                    minus_median_dp = adf[adf.strand == "-"].median_dp.groupby(adf.amplicon).max(),
                    plus_median_dp = adf[adf.strand == "+"].median_dp.groupby(adf.amplicon).max(),
                   )
    
    df = pd.DataFrame(df_dict) 
    df = df.reset_index()
    df = df.fillna(0)    

    df["total_pass_ratio"] =  df.total_pass / (df.size.astype(float) * 2.0)
    df["strand_dp_ratio"] =  abs(1- ( df.minus_median_dp / df.plus_median_dp.astype(float) ))
    
    df["pass_ratio_ok"] = df.total_pass_ratio > pass_cutoff
    df["strand_ratio_ok"] = df.strand_dp_ratio < strand_diff_cutoff
    
    df = df.fillna(0)    
    
    ##################################################################
    # print flagged amplicons
    
    stem = bam.split(".")[0].strip("_coverage")
    coverage_output =  stem + "_undercovered_amplicons.tsv"
    strand_output = stem + "_strandbiased_amplicons.tsv"

    df[~df.pass_ratio_ok].to_csv(coverage_output, sep = "\t") 
    df[~df.strand_ratio_ok].to_csv(strand_output, sep = "\t")
    
    total_amplicons = float( df.name.nunique() )
    
    amplicons_passing_coverage_cutoff = len(df[~df.pass_ratio_ok])
    coverage_pass_ratio = round( (amplicons_passing_coverage_cutoff/total_amplicons) * 100, 2)
    print "%s percent (%s/%s) of amplicons in %s pass the coverage cutoff." % (coverage_pass_ratio, amplicons_passing_coverage_cutoff, total_amplicons, bam)
    
    amplicons_passing_strand_cutoff = len(df[~df.strand_ratio_ok])
    strand_pass_ratio = round( (amplicons_passing_strand_cutoff/total_amplicons) * 100, 2)
    print "%s percent (%s/%s) of amplicons in %s pass the strand cutoff." % (strand_pass_ratio, amplicons_passing_strand_cutoff, total_amplicons, bam)

    ##################################################################
    # generate graphs
    
    if graphs:
        plt.figure()
        plt.hist(df[df.strand_dp_ratio < 1].strand_dp_ratio, bins = 40)
        plt.xlabel("Difference between median coverage on left vs right strand ")
        plt.ylabel("Amplicon count")
        graphname = stem + "_strand_comparison_median_dp_histogram.png"
        plt.savefig(graphname, dpi = 300)
        
        plt.figure()
        plt.hist(df.total_pass_ratio, bins = 20)
        plt.xlabel("Ratio of bases above coverage cutoff")
        plt.ylabel("Amplicon count")
        graphname = stem + "_bases_above_cutoff_histogram.png"
        plt.savefig(graphname, dpi = 300)
        
        plt.figure()
        plt.scatter(df.plus_median_dp,df.minus_median_dp) 
        plt.xlabel("Median coverage of bases on PLUS strand")
        plt.ylabel("Median coverage of bases on MINUS strand")
        graphname = stem + "_strand_comparison_median_dp_scatterplot.png"
        plt.savefig(graphname, dpi = 300)
        
        plt.figure()
        plt.scatter(df.plus_pass_ratio, df.minus_pass_ratio) 
        plt.xlabel("Ratio of bases above cutoff on PLUS strand")
        plt.ylabel("Ratio of bases above cutoff on MINUS strand")
        graphname = stem + "_strand_comparison_bases_above_cutoff_scatterplot.png"
        plt.savefig(graphname, dpi = 300)
    
def run_bedtools(bam, output, bed = None):
    bedtools_cmd = "bedtools coverage -s -d -abam %s -b %s > %s" % (bam, bed, output)
    #
    print bedtools_cmd
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

def main():
    
    
    bed = sys.argv[1]
    bed = check_bed_filetype(bed)
    
    print "-" * 50
    
    dp_cutoff = 100
    pass_cutoff = 0.8
    strand_diff_cutoff = 0.5
    
    print "Minimum coverage per base: %sX \nMinimum ratio of accepted bases per amplicon: %s" % (dp_cutoff, pass_cutoff)
    print "Maximum coverage difference between strands: %s" % (strand_diff_cutoff)
    
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
                print "Already present:", bedtools_output
            
            find_undercovered_regions(bedtools_output, graphs = False, dp_cutoff = dp_cutoff, pass_cutoff = pass_cutoff, strand_diff_cutoff = strand_diff_cutoff)
            print "Finished", bam
            print "---------------------------------------------------"
    


##################################################################################

if __name__ == "__main__":
    main()

    
    