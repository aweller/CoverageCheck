# Copyright (C) 2014 Andreas M. Weller <andreas.m.weller@gmail.com>
#
# read a bedtools output file from
#
# bedtools coverage -abam Q2PL2_H01_N.bam -b TSB_148_gene_panel_HP_amplicons.bed -d > test_coverage.csv
#
# and find bases with coverage or strand_ratio below the threshold

# general modules
import pandas as pd
import sys
import subprocess
import os
import matplotlib.pylab as plt
import numpy as np
import pprint
import logging
import re

import seaborn as sns
sns.set(font="serif") # prevents error message about fonts

# personal modules from the same folder
import plot_exon_coverage_v2 as plotting
import plot_exon_coverage_all_samples as all_sample_plotting
import unite_coverage_files as UniteCoverage
import CoverageCheckClasses
from CoverageCheckConfig import *

script_folder = os.path.dirname(os.path.realpath(__file__)) + "/"

DISTNAME = 'coveragecheck'
MAINTAINER = 'Andreas Weller'
MAINTAINER_EMAIL = 'andreas.m.weller@gmail.com'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/aweller/CoverageCheck/'
VERSION = '0.1'

#####################################################################################################################

def parse_exons_into_dataframe_and_dict(exon_filename):
    """
    Parse a list of all human exons into 2 datastructures:
    1. exons (pandas df)
    2. exons_per_gene (dictionary of gene:[start, stop, exon_no])
    """
    
    #exon_filename = "/home/andreas/bioinfo/core/general/data/HumanExons_Ensembl_v65_merged.csv"
    #exon_filename = "/home/andreas/bioinfo/core/general/data/HumanExons_Ensembl_v75_all_genes_merged.csv"

    header = ["chrom", "exon_start", "exon_stop", "gene", "exon_no", "strand"]
    exons = pd.read_csv(exon_filename, sep="\t", names=header)
    exons["gene_upper"] = exons.gene.str.upper()
    exons = exons.sort(columns = ["gene", "exon_start", "exon_stop"])
        
    global exons_per_gene
    exons_per_gene = {}
    for _, row in exons.iterrows():
        gene = row["gene"].upper()
        start, stop = int(row["exon_start"]), int(row["exon_stop"])
        exon_no = row["exon_no"]
        
        if not exons_per_gene.get(gene):
            exons_per_gene[gene] = []
            
        exons_per_gene[gene].append((start, stop, exon_no))
    
    return exons, exons_per_gene

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

def parse_amplicons_into_df(bed, exons):
    """
    Parse the amplicons in the bedfile into a DF
    Get rid of amplicons that are too far away from any exon (as defined in the exon file)
    """
    
    header = ["chrom", "start", "stop", "amplicon", "na", "strand"]
    
    raw_amplicons = pd.read_csv(bed, sep="\t", names=header)
    raw_amplicons["gene"] = raw_amplicons.apply(lambda x: x["amplicon"].split("_")[0].upper(), axis=1)
    raw_amplicons["altgene"] = raw_amplicons.apply(lambda x: x["amplicon"].split("_")[1].upper(), axis=1)

    def crosscheck_gene_chromosome(row):
        """ Check if this row is close enough to an exon """
        offset = 5000
        
        gene = row["gene"]
        altgene = row["altgene"]
        chrom = row["chrom"]
        start = row["start"]
        stop = row["stop"]
        
        gexons = exons[exons.gene == gene]
        
        if len(gexons) > 0:
            altgexons = exons[exons.gene == altgene]
            if len(gexons) < len(altgexons):
                gexons = altgexons
            
            chrom_ok = chrom == gexons.chrom.unique()[0]
            start_ok = start > gexons.exon_start.min() - offset
            stop_ok = stop < gexons.exon_start.max() + offset
            
            return all([chrom_ok, start_ok, stop_ok])
        else:
            #print "gexons empty for", gene
            return False
    
    raw_amplicons["correct_chrom"] = raw_amplicons.apply(crosscheck_gene_chromosome, axis=1)
    amplicons = raw_amplicons[raw_amplicons.correct_chrom]
    amplicons = amplicons.drop_duplicates(cols=["start", "stop"]) # as in the  'separate strand' beds each exon is present twice
    
    return amplicons

#####################################################################################################################

def find_bad_positions(coverage_matrix, target_folder = None, trait = None, samplename = None,
                       trait_cutoff = None, whitelist = None):    
    
    """
    Walk through all bases and find contigous regions of bases that fail the coverage/strandbias cutoff
    Create an output file of these regions
    If a base an expected variant (recorded in a ExpectedVariants instance ('whitelist'), record the trait to the instance)
    Return the updated whitelist
    """
    
    region = None
    last_gene = None
    last_pos = 0
    last_gene = None
    
    total_bases = 0
    bad_bases = 0
    
    sample = samplename
    bad_output_name = target_folder + sample + "_failed_regions_%s_cutoff_%s.csv" % (trait, trait_cutoff)
    bad_output = open(bad_output_name, "w")
    
    if trait == "strandbias":
        output_header = ["gene", "chrom", "start", "stop", "mean_strand_bias", "size"]
    elif trait == "coverage":
        output_header = ["gene", "chrom", "start", "stop", "mean_coverage", "size"]
    
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
        sampleinfo.add_strandbias(sample, [pass_percent, good_bases, total_bases])
        logging.info( "%s percent (%s/%s) of positions in %s have a strand bias below the threshold (%s:1)."
                     % (pass_percent, good_bases, total_bases, sample, trait_cutoff) )
        
    elif trait == "coverage":
        sampleinfo.add_coverage(sample, [pass_percent, good_bases, total_bases])
        logging.info( "%s percent (%s/%s) of positions in %s have at least the minimum coverage of %sX."
                     % (pass_percent, good_bases, total_bases, sample, trait_cutoff) )
    
    return whitelist
    
#####################################################################################################################

def run_bedtools_coverage(bam, output, bed = None):
    """ Run bedtools coverage to get the per-base coverage in the target bam. """
    
    bedtools_cmd = "bedtools coverage -s -d -abam %s -b %s > %s" % (bam, bed, output)
    logging.debug( bedtools_cmd )
    output_code = subprocess.call(bedtools_cmd, shell=True)
    
    if output_code != 0:
        logging.critical( "Bedtools coverage run error." )
        logging.critical( bedtools_cmd )
        logging.critical( "Sorry, aborting..." )
        sys.exit()

def run_bedtools_intersect(bed):
    """ Run bedtools intersect to reduce the exon file to the exons with an amplicon overlap. """
    
    output = bed.replace(".bed", "_covered_exon_locations.bed") 
    
    bedtools_cmd = "intersectBed -u -a %s/input/%s -b %s > %s" % (script_folder, human_exon_file, bed, output)
    logging.debug( bedtools_cmd )
    output_code = subprocess.call(bedtools_cmd, shell=True)
    
    if output_code != 0:
        logging.critical( "Bedtools intersect run error." )
        logging.critical( bedtools_cmd )
        logging.critical( "Sorry, aborting..." )
        sys.exit()
        
    return output

def check_bed_filetype(filename):
    """
    Check if the bedfile is actually a bed and convert it to one if its an Illumina Manifest file
    """
    
    rows = open(filename).readlines()
    
    if not "Header" in rows[0]: 
        
        if all([len(row.split("\t")) == 6 for row in rows]):
            logging.info( "Correct bed file detected: "+ filename )
            return filename
        
        else:
            logging.critical( "Error in input bed file:", filename )
            logging.critical( "Not all rows contain 6 fields as is expected." )
            logging.critical( "Sorry, aborting..." )
            sys.exit()
    
    else: # this is an Illumina manifest file
        
        logging.error( "Input is not a bed file, but an Illumina manifest file: " + filename )
        
        expected_bed = filename[:-4] + "_plusminus.bed"
                
        if os.path.exists(expected_bed):
            logging.info( "Switching to existing bed file: " + expected_bed )
            return expected_bed
        
        else:
            logging.info( "Converting manifest to bed." )
            
            import manifest2bed as m2b
            all_out, plusminus = m2b.convert_manifest(filename)
            
            logging.info( "Switching to newly created bed file: " + plusminus )
            return plusminus

def fix_gene_names_in_bedfile(bed, gene_alias_filename):
    """
    Creates a new bedfile that has all genenames changed according to the gene_aliases file
    Return the name of the new bedfile
    """
    
    fixed_bed_name = bed.replace(".bed", "_fixed_genenames.bed")
    
    logging.info( "Gene Name Aliases detected: " + gene_alias_filename )
        
    if os.path.exists(fixed_bed_name):
        logging.info( "Switching to existing fixed bed file: " + fixed_bed_name )
        return fixed_bed_name
    
    else:
        logging.info( "Fixing gene names in bedfile." )
    
        gene_aliases = {}
        with open(gene_alias_filename) as handle:
            for row in handle:
                row = row.strip()
                old, new = re.split("[\t ]", row) # split the row regardsless of tab or space
                gene_aliases[old] = new
        
        out = open(fixed_bed_name, "w")
        with open(bed) as handle:
            for row in handle:
                f = row.strip().split("\t")
                amplicon = f[3]
                gene = amplicon.split("_")[0]
                new_gene = gene_aliases.get(gene, gene)
                
                new_amplicon = amplicon.replace(gene, new_gene)
                
                f[3] = new_amplicon
                result = "\t".join(f)
                
                out.write(result + "\n")
                
                if amplicon != new_amplicon:
                    logging.debug("Replaced %s with %s" % (amplicon, new_amplicon))
        out.close()
        return fixed_bed_name

def remove_empty_files_from_folder(folder):
    for filename in os.listdir(folder):
        try:
            if os.path.getsize(folder +"/"+ filename) == 0:
                os.remove(folder +"/"+ filename)
        except:
            pass
        
#####################################################################################################################
#####################################################################################################################

def run(bed, target_folder, min_dp, max_strand_ratio, whitelist_filename=None, gene_alias_filename=None, target_bams_filename=None):
    
    if whitelist_filename == "None":
        whitelist_filename = None
    if gene_alias_filename == "None":
        gene_alias_filename = None
    if target_bams_filename == "None":
        target_bams_filename = None

    remove_empty_files_from_folder(target_folder) # remove empty files that might have been left over from previous runs
    
    if " " in target_folder:
        logging.critical( "The third-party tools used by CoverageCheck don't accept spaces in folder names." )
        logging.critical( "Please replace the spaces in your bam folder with underscores.")
        logging.critical( "Sorry, aborting..." )
        
    ##############################################################################################
    # configure logging to both sys.stdout and a file 
    
    logging_filename = "CoverageCheck_log.txt"
    print "Log messages printed to %s" % (logging_filename)
    
    # set up logging to file
    logging.basicConfig(level=logging.INFO,
                  format='%(asctime)s %(name)-8s %(levelname)-8s %(message)s',
                  datefmt='%d-%m-%y %H:%M',
                  filename=logging_filename)
    
    # set up logging to console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    logging.info( "-" * 100 )
    logging.info( "Minimum accepted coverage per base: %sX" % (min_dp) )
    logging.info( "Maximum accepted coverage ratio between strands: %s:1" % (max_strand_ratio) )
        
    #####################################################
    # setup variables and report information 
     
    bed = check_bed_filetype(bed)
    if gene_alias_filename:
        bed = fix_gene_names_in_bedfile(bed, gene_alias_filename)
    
    # parse all exons and all proper amplicons into DataFrames
    exon_filename = run_bedtools_intersect(bed)
    exons, exons_per_gene = parse_exons_into_dataframe_and_dict(exon_filename)

    amplicons = parse_amplicons_into_df(bed, exons)
        
    # setup the class that collects stats on each sample    
    global sampleinfo
    sampleinfo = CoverageCheckClasses.SampleInfo()
    
    # parse list of expected variants
    if whitelist_filename:
        variant_no = len([x for x in open(whitelist_filename).readlines() if x[0] != "#"])
        logging.info( "Found %s expected variants to check in %s" % (variant_no, whitelist_filename) )
        
        if variant_no == 0:
            logging.critical( "ERROR: expected variant file doesn't contain variants." )
            sys.exit()
    else:
        logging.info( "No expected variants specified." )

    # parse list of allowed bams
    target_folder = target_folder + "/" # just to be on the safe side    
    bams = [x for x in os.listdir(target_folder) if x.endswith(".bam")]    
    
    logging.info( "Found %s bams in the target folder %s" % (len(bams), target_folder))    
    if target_bams_filename:
        target_bams = [x.strip() for x in open(target_folder + target_bams_filename).readlines()]
        target_bams = [x+".bam" for x in target_bams if not x.endswith(".bam")]
        logging.info( "Parsed list of %s allowed bams in %s" % (len(target_bams), target_bams_filename) )
        bams = [x for x in bams if x in target_bams]
        logging.info( "%s bams left to process in the target folder." % len(bams) )    

    logging.info( "-" * 100 )

    ##############################################################################################
    # process each samples separately
    
    for bam in bams:
        
        samplename = bam.split(".")[0]
        sample_output_folder = "%s/%s_results/" % (target_folder, samplename)
        
        sampleinfo.add_sample(samplename)
        
        if not os.path.exists(sample_output_folder):
            os.makedirs(sample_output_folder)
        
        ##############################################################################################
        # create the coverage file (if necessary) and parse it into a pandas DF
        
        bedtools_output = target_folder + samplename + "_coverage.csv"
        
        if not os.path.exists(bedtools_output):
            logging.info( "Running bedtools coverage for: " + bam )
            run_bedtools_coverage(target_folder + bam, output = bedtools_output, bed=bed)
        
        coverage_matrix = parse_coverage_file_into_dataframe(bedtools_output)
        
        ##############################################################################################
        # Initialize the class instance thats collects information on each expected variant
        
        whitelist = None
        if whitelist_filename:
            whitelist = CoverageCheckClasses.ExpectedVariants(whitelist_filename, samplename=samplename, folder=sample_output_folder,
                                        dp_cutoff = min_dp, strandbias = max_strand_ratio, exon_dict = exons_per_gene)
        
        ##############################################################################################
        # run
        
        if output_undercovered_regions:
            whitelist = find_bad_positions(coverage_matrix, target_folder = sample_output_folder, trait = "coverage",
                               samplename = samplename, trait_cutoff = min_dp, whitelist=whitelist)
        
        if output_strandbiased_regions:
            whitelist = find_bad_positions(coverage_matrix, target_folder = sample_output_folder, trait = "strandbias",
                               samplename = samplename, trait_cutoff = max_strand_ratio, whitelist=whitelist)
        
        if output_undercovered_regions or output_strandbiased_regions:
            if whitelist_filename:
                whitelist.print_output()
            
        ##############################################################################################
        # plot
        
        if output_coverage_plots:
            sample_df = plotting.create_all_coverage_plots(bedtools_output, exons=exons, exons_per_gene = exons_per_gene,
                                        target_folder = sample_output_folder, whitelist=whitelist, amplicons=amplicons)
        logging.info(  "-" * 100 )
        
    ##############################################################################################
    # create summaries across all samples
    # this has to happen after the individual bam treatment because it assumes that coverage files for each bam are present
    
    sampleinfo.print_output()
    
    all_sample_filename = "all_samples.csv"
    UniteCoverage.unite(all_sample_filename, target_folder = target_folder)
    
    byte_size = os.path.getsize(target_folder+all_sample_filename)
    mb_size = byte_size/1.049e+6
    
    if output_all_sample_summary:
        if mb_size > 2000:
            logging.error( "The united coverage size of all samples is %sMB. Skipping summary plots." % (round(mb_size,2)) )     
        else:
            all_sample_plotting.plot_exon_coverage(all_sample_filename, exons=exons, exons_per_gene = exons_per_gene,
                                                   target_folder = target_folder, whitelist=whitelist)
        
##################################################################################

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bed", help="region file in bed or Illumina Manifest format",
                        required=True)

    parser.add_argument("-x", "--expected_variants", help="expected variants in vcf or HGMD BioMart export format")
    parser.add_argument("-a", "--alias", help="gene name aliases, format: [old new]")
    parser.add_argument("-w", "--whitelist", help="list of bams to analyse")

    parser.add_argument("-c", "--min_coverage",  type=int, help="minimum coverage (Default: 50X)", default=default_min_coverage)
    parser.add_argument("-s", "--max_strandratio", type=float, help="maximum strand ratio (Default: 5)", default=default_max_strandbias)
    
    args = parser.parse_args()
    pprint.pprint(args)
       
    ##################################################################################
    
    bed = args.bed
    whitelist_filename = args.expected_variants
    gene_alias_filename = args.alias
    target_bams_filename = args.whitelist
    min_dp = args.min_coverage
    max_strand_ratio = args.max_strandratio
    
    if "/" in bed:
        target_folder = "/".join(bed.split("/")[:-1]) + "/"
    else:
        target_folder = "./"
    
    run(bed,
        target_folder,
        min_dp,
        max_strand_ratio,
        whitelist_filename=whitelist_filename,
        gene_alias_filename=gene_alias_filename,
        target_bams_filename=target_bams_filename)

    
    