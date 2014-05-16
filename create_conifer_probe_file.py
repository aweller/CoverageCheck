# convert a bed file into a conifer input file

import sys

bed = sys.argv[1]

print "\t".join(["chr","start","stop","name"])

blacklist = ["2","6","8", "14", "18"]

with open(bed) as handle:
    for row in handle:
        f = row.strip().split("\t")
        
        gene = f[3].split("_")[0]
        
        out = f[:3]
        out.append(gene)
        out = "\t".join(out)
        
        if f[0][3:] not in blacklist:
            print out.replace("chr", "")