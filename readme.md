#CoverageCheck

##Introduction

This script was mainly written for a clinical diagnostics setting. The use case is where a group of patients was sequenced using a panel. It is crucial to distinguish between 

A) high-confidence positions with a reference allele and

B) low-confidence positions where a variant call is impossible.

Without coverage information, A) and B) would look the same on a clinical report: "no variants found". For clinical decision making though, it's quite a difference if a patient is wildtype on TP53 or if we can't make a statement about TP53.

High-confidence positions are defined as having a coverage above the minimum accepted coverage (default 50X) and a ratio of coverage on minus vs plus strand (default 5x)
below the maximum accepted strand ratio.

#####################################################################################

##Usage

A GUI is started via a Desktop icon or directly from the command line

> *"python CoverageCheck_GUI_tkinter.py".*

CoverageCheck can also be started directly from the command line:

> *"python CoverageCheck.py regions.bed expected_variants.vcf gene_aliases.txt accepted_bams.txt"*

##Input

CoverageCheck assumes that all input files are located in the same folder.
The minimum accepted coverage and the maximum accepted strand ratio are set in the GUI.

### Required
- location of the bam folder

- region file (Bed format or Illumina Manifest)
CoverageCheck parses the gene names from amplicon names of type 'SMAD4_Exon_(2023461)_2777650' (as is the standard for Illumina custom panels).

### Optional
- expected variants (vcf format or export format of HGMD mutation mart)
Specifiying this will lead to the creation of additional output files listing the 'theoretical discoverability' for each variant for the given cutoff values. Expected variants per exon will also be listed in the plots.

- gene alias file (txt file, one gene per row, format: old_name new_name)
This option is a simple hack to deal with amplicons that have gene names that differ from the offical names listed in the HumanExons_V75(...) file.  

- allowed bams file (txt file, one bam per row)
This option allows to analyze not all bams in the target folder by explicitly specifiying the bams to be analyzed. 

##Output

For each sample, CoverageCheck creates a results folder with the following output files:

1. [samplename]_coverage_raw_exon_coverage.png

    A graph with a dot for each base, showing the coverage across exons per gene. Exons plotted on their actual chromosome location.

2. [samplename]_coverage_raw_exon_coverage.png

    A graph with a boxplot for each exon, showing the coverage distribution across exons per gene. Exons are simply plotted next to each other.

3. [samplename]_failed_regions_coverage.tsv

    A list of bases that fail the coverage cutoff. Continous streches of bases are summarized as a region. 

4. [samplename]_failed_regions_strandbias.tsv

    A list of bases that fail the strandbias cutoff. Continous streches of bases are summarized as a region. 

If a list of expected variants is provided, two additional output files are produced:

5. [samplename]_expected_variant_coverage.tsv

    A list of all variant position in the original vcf, with columns added for coverage, strandbias and whether the position passes the cutoffs.

5. [samplename]_expected_variant_coverage_per_gene.tsv

    A list of all genes that contain at least one expected variant. For each gene, the total expected variants, the variants that have sufficient quality to be         theoretically discoverable and the ratio fo discoverable to total is listed.

#####################################################################################

##Dependencies

All bams to be analyzed need to be in a folder, together with the  file specifying the sequenced regions (in bed or Illumina manifest format). 

CoverageCheck expects a Linux system with Python 2.7 and bedtools installed. The following Python packages are expected as well: numpy, pandas, matplotlib and seaborn.

The bam files need to be indexed (i.e. have a *.bam.bai file in the same folder). If missing, this can be created with "samtools index sample.bam"