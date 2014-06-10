#####################################################
# Constant input files

#expected format:
#chr17	7590695	7590856	TP53	1	-
human_exon_file = "HumanExons_Ensembl_v75_refseqs.bed"

#expected format:
#chr17	7565097	7590856	TP53
human_gene_location_file = "human_gene_locations.txt"

#####################################################
# Default values

default_min_coverage = 100
default_max_strandbias = 5.0

#####################################################
# Flags

# create the output tables for failed region?
output_undercovered_regions = True
output_strandbiased_regions = False

# create plots at all?
output_coverage_plots = True

# create the individual plots for each gene?
output_per_gene_plots = True

# create the (still ugly) summary plots for all samples?
# this might be slowing the computer down
output_all_sample_summary = True