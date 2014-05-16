import easygui as eg
import sys        
import os
import time
import CoverageCheck as CC

def run_CC():
    min_dp = int(fieldValues[0])
    max_strand_ratio = int(fieldValues[1])
    target_folder = fieldValues[2]
    bed = fieldValues[3]

    CC.run(bed, target_folder, min_dp, max_strand_ratio)


###################################################################################
# collect values

msg         = "Please enter CoverageCheck parameters"
title       = ""
fieldNames  = ["Min. coverage", "Max. strandbias", "Bam folder", "Target region file"]
fieldValues = [50, 5, "/home/andreas/bioinfo/projects/nhs_coverage_tool/data/Alignment",
               "/home/andreas/bioinfo/projects/nhs_coverage_tool/data/Alignment/CLL_Dis_V1_TruSeq_CAT_Manifest_TC0031836-CAT.bed"]  # we start with blanks for the values
fieldValues = eg.multenterbox(msg,title, fieldNames, values=fieldValues)

###################################################################################
# look for errors in the returned values


while 1:  # do forever, until we find acceptable values and break out

    errmsg = ""
    
    for i in range(len(fieldNames)):
        if fieldValues[i].strip() == "":
            errmsg = errmsg + ('"%s" is a required field.\n\n' % fieldNames[i])
    
    if not os.path.exists( fieldValues[2] ):
        errmsg += ('"%s" is not a valid folder.\n\n' % fieldNames[2])
    
    if not os.path.exists( fieldValues[3] ):
        errmsg += ('"%s" not found.\n\n' % fieldNames[3])        


    if errmsg == "":
        run_CC()
        break
    else:
        # show the box again, with the errmsg as the message    
        fieldValues = eg.multenterbox(errmsg, title, fieldNames, fieldValues)
        
    #    choice_msg = "Minimum coverage: %s \nMaximum strand bias: %s \nBam folder: %s \nRegion file: %s \n" % (fieldValues[0], fieldValues[1], fieldValues[2], fieldValues[3])
    #    title = "Start CoverageCheck?"
    #    
    #    if eg.ccbox(choice_msg, title):     # show a Continue/Cancel dialog
    #        break  # user chose Continue
    #    else:  # user chose Cancel
    #        sys.exit(0)
