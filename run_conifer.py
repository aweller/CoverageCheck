"""
# create RPKM for sample
python ~/bin/conifer_v0.2.2/conifer.py rpkm
--probes bams_14MiSeqRCP/OHTC_CAT_Manifest_TC0023302-CAT_conifer_input.txt
--input bams_14MiSeqRCP/G158265J_S5.bam
--output confier_14MiSeqRCP/rpkm/G158265J_S5.rpkm.txt


python ~/bin/conifer_v0.2.2/conifer.py analyze \
--probes bams_14MiSeqRCP/OHTC_CAT_Manifest_TC0023302-CAT_conifer_input.txt \
--rpkm_dir confier_14MiSeqRCP/rpkm/ \
      --output analysis.hdf5 \
      --svd 6 \
      --write_svals singular_values.txt \
      --plot_scree screeplot.png \
      --write_sd sd_values.txt

python conifer.py analyze \
      --probes probes.txt \
      --rpkm_dir ./RPKM/ \
      --output analysis.hdf5 \
      --svd 6 \
      --write_svals singular_values.txt \
      --plot_scree screeplot.png \
      --write_sd sd_values.txt

python conifer.py call \
      --input analysis.hdf5 \
      --output calls.txt

python conifer.py plot \
      --input analysis.hdf5 \
      --region chr#:start-stop \
      --output image.png \
      --sample sampleID
"""

#run conifer on a set of samples (bams in a folder)

import sys
import os
import subprocess as sp
import multiprocessing as mp

target_folder = sys.argv[1]
target_probe = sys.argv[2]
output_folder = sys.argv[3]

bams = [x for x in os.listdir(target_folder) if x.endswith(".bam")]

###################################################################################
# calculate RPKM values

def call_cmd(cmd):
    sp.call(cmd, shell=True)

parallel = True

if parallel:
    pool = mp.Pool(4)

for bam in bams:
    sample = bam[:-4]
    output = "%s/rpkm/%s.rpkm.txt" % (output_folder, sample)
    
    rpkm_cmd = """python ~/bin/conifer_v0.2.2/conifer.py rpkm \
                    --probes %s \
                    --input %s/%s \
                    --output %s""" % (target_probe, target_folder, bam, output)
                    
    
    if os.path.exists(output):
        print "%s exists already, skipping..." % sample
    else:
        print rpkm_cmd
        print "Starting %s calculation..." % sample
        
        if parallel:
            pool.apply_async(call_cmd, args = [rpkm_cmd, ])
        else:
            call_cmd(rpkm_cmd)
        
if parallel:
    pool.close()
    pool.join()
    
###################################################################################
# main analysis

analysis_cmd = """python ~/bin/conifer_v0.2.2/conifer.py analyze \
                    --probes %s \
                    --rpkm_dir %s \
                    --output %s/analysis.hdf5 \
                    --svd 4 \
                    --min_rpkm 1 \
                    --write_svals singular_values.txt \
                    --plot_scree screeplot.png \
                    --write_sd sd_values.txt""" % (target_probe, output_folder + "/rpkm/", output_folder)

print "-" * 100
print analysis_cmd
sp.call(analysis_cmd, shell=True)

###################################################################################
# CNV calling

call_cdm = """python ~/bin/conifer_v0.2.2/conifer.py call \
                --input %s/analysis.hdf5 \
                --output %s/calls.txt \
                --threshold 0.4""" % (output_folder, output_folder)

print "-" * 100
print call_cdm
sp.call(call_cdm, shell=True)

###################################################################################
# plotting

plot_cmd = """python ~/bin/conifer_v0.2.2/conifer.py plotcalls \
                --input %s/analysis.hdf5 --outputdir %s/call_imgs/ \
                --calls %s/calls.txt """ % (output_folder, output_folder, output_folder)
                
print "-" * 100
print plot_cmd
sp.call(plot_cmd, shell=True)               

