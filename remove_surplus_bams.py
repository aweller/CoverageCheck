# Copyright (C) 2014 Andreas M. Weller <andreas.m.weller@gmail.com>
#
# remove all files from current folder that don't start with samplename defined in sysargv file

import os
import sys

valid_samples = open(sys.argv[1]).readlines()
valid_samples = [x.strip() for x in valid_samples]

print valid_samples

folder = [x for x in os.listdir(".") if "bam" in x]

for filename in folder:
    sample = filename.split(".")[0].split("_")[0]
    
    if not sample in valid_samples:
        print "removing " + filename
        os.remove(filename)
    else:
        print "keeping " + filename
