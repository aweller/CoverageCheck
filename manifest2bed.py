# Copyright (C) 2014 Andreas M. Weller <andreas.m.weller@gmail.com>
#
# convert an Illumina manifest file to a bed file

import sys

def convert_manifest(manifest_file):
    
    stem = manifest_file.split(".")[0]
    all_out_file = stem + ".bed"
    plusminus_out_file = stem + "_separate_strands.bed"
    
    all_out = open(all_out_file, "w")
    plusminus = open(plusminus_out_file, "w")
    
    with open(manifest_file) as handle:
        for row in handle:
            f = row.strip().split("\t")
            
            if row[0] == "[": # new section
                section = row.strip().strip("[").strip("]")
            
            elif "Chromosome" in row: # header for this section
                header = f
            
            else:
                #print section, row,
                if section == "Targets":
                    name, chrom, start, stop, actual_strand = f[1], f[3], f[4], f[5], f[6]
                    name =  name.replace(" ", "_")
                    
                    bed_row = "\t".join([chrom, start, stop, name, ".", actual_strand])
                    all_out.write(bed_row + "\n")
                    
                    plus_row = "\t".join([chrom, start, stop, name, ".", "+"])
                    plusminus.write(plus_row + "\n")
                    
                    minus_row = "\t".join([chrom, start, stop, name, ".", "-"])
                    plusminus.write(minus_row + "\n")
                    
    all_out.close()
    plusminus.close()
    
    return all_out_file, plusminus_out_file

def main():
    manifest_file = sys.argv[1]
    convert_manifest(manifest_file)

##################################################################################

if __name__ == "__main__":
    main()