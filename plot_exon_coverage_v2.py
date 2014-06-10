import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import scipy.stats
from IPython.html.widgets import interact
import os
import logging

from CoverageCheckConfig import *

# seaborn default
blue = (0.2980392156862745, 0.4470588235294118, 0.6901960784313725)
black = [ 0.22196078,  0.1945098 ,  0.27843137,  1.        ] # from sns.dark_palette("MediumPurple")
red = (0.7686274509803922, 0.3058823529411765, 0.3215686274509804)
sns.set_style("darkgrid")
sns.set_context("talk")

##################################################################################################################
##################################################################################################################

def find_matching_exon_from_dict(row, *args):
    exons_per_gene = args[0]
    
    gene = row["gene"]
    pos = row["pos"]

    exon_no = 0
    
    target_exons = exons_per_gene.get(gene)
    
    if target_exons:
        for region in target_exons:
            if region[0] < pos < region[1]:
                exon_no = region[2]
                break
    
    return exon_no     
            
def plot_exon(row):
    start = int(row["exon_start"])
    stop = int(row["exon_stop"])
    size = stop - start    
    
    #print start, stop
    
    rectangle = plt.Rectangle((start, -20), size, 10, fc='red')
    plt.gca().add_patch(rectangle)

def plot_exon_relpos(row):
    start = pos_dict[int(row["exon_start"])]
    stop = pos_dict[int(row["exon_stop"])]
    size = stop - start    
        
    rectangle = plt.Rectangle((start, -1000), size, 2000, fc=blue)
    plt.gca().add_patch(rectangle)
    
def plot_amplicons_relpos(row):
    #ax = plt.subplot(111)
    
    start = int(row["start"])
    stop = int(row["stop"])
    
    color = red
    
    global switch
    switch = switch == False

    height = 500
    if switch:
        y_start = -520
    else:
        y_start = 20
    
    #try:
    if True:
        try:
            relstart = pos_dict[start]
        except:
            key, relstart = min(pos_dict.items(), key=lambda (v, _): abs(v - start))
            #print "Replaced start %s with the closest: %s" % (start, key)
            #color = "cyan"
        try:
            relstop = pos_dict[stop]
        except:
            key, relstop = min(pos_dict.items(), key=lambda (v, _): abs(v - stop))
            #print "Replaced stop %s with the closest: %s" % (stop, key)
            #color = "cyan"
        
        if relstart != relstop:
            size = stop - start
            relsize = relstop - relstart
            
            rectangle = plt.Rectangle((relstart, y_start), relsize, height, fc=red, ec=black)
            plt.gca().add_patch(rectangle)
            
def plot_raw_exon_coverage(gene, save_per_gene=False, df=None,
                           exons=None, plot_no=None, all_genes=None,
                           amplicons=None, per_gene_output_folder=None, whitelist=None):
    
    pdf = df[df.gene == gene]
        
    #########################################################
    # prepare gene exons and the mapping from chrom pos to adjusted plot pos
    # introns are shortened if they are longer than the offset
    
    global gexons
    gexons = exons[exons.gene == gene]
    
    if not len(gexons) > 0: # this is an empty placeholder
        logging.debug( "Exon list empty for " + gene )
        if not save_per_gene:
            plt.subplot(len(all_genes), 1, plot_no)
            plt.ylabel(gene)
            plt.xticks([],[])
            plt.yticks([],[])
            plt.text( 0.5, 0.5, 'No exons defined')

        return
    
    plot_range = gexons.exon_stop.max() - gexons.exon_start.min()
    if plot_range < 10000:
        offset = 1000
    elif plot_range < 30000:
        offset = 200
    else:
        offset = 100
    
    gexons["min_start"] = gexons.exon_start - offset
    gexons["max_stop"] = gexons.exon_stop + offset
    
    global pos_dict
    pos_dict = {}
    rel_count = 0
    for index, row in gexons.iterrows():
        start = row["min_start"]
        stop = row["max_stop"]    
        for i in range(start, stop+1):
            if not pos_dict.get(i):
                rel_count += 1
                pos_dict[i] =  rel_count
            else:
                pass
                #print i, rel_count, pos_dict.get(i)
    
    pdf["relpos"] = pdf.apply(lambda x: pos_dict.get(x["pos"]), axis=1)
    pdf = pdf[pdf.relpos > 0]
    
    #########################################################
    # setup plot area
            
    if save_per_gene:
        plt.figure(figsize=(20, 5))
    else:
        plt.subplot(len(all_genes), 1, plot_no)

    plt.ylabel(gene)
    
    target_chrom = gexons.chrom.unique()[0]
    x_start = min(pos_dict.values())
    x_stop = max(pos_dict.values())
    plt.xlim(x_start, x_stop)

    global gamplicons
    gamplicons = amplicons[(amplicons.gene == gene) | (amplicons.altgene == gene)]
    
    global switch
    switch = True
    gexons.apply(plot_exon_relpos, axis =1)
    gamplicons.apply(plot_amplicons_relpos, axis=1)

    ylim = max( [max(pdf.plus_dp), max(pdf.minus_dp)] )*1.2
    if ylim < 5000:
        ylim = 5000    
    plt.ylim(-ylim,ylim)

    ############################################
    # create actual plot content
    
    center = 1000

    y = pdf.plus_dp
    y = y + center
    x = pdf.relpos

    _ = plt.plot(x, y, color=black)
    plt.fill_between(x, y, 1000,color=black)

    y = pdf.minus_dp
    y = y + center
    y = y * -1
    x = pdf.relpos

    _ = plt.plot(x, y, color=black)
    plt.fill_between(x, y, -1000,color=black)

    plt.axhline(center, color=black)
    plt.axhline(-center, color=black)
    
    ############################################
    # set xticks
    
    gexons["center"] = (gexons["exon_start"] + gexons["exon_stop"])*0.5
    
    xlabels = []
    xlocs = []
    
    for index, row in gexons.iterrows():
        location = int(row["center"])
        rel_loc = pos_dict[location]
        xlocs.append(rel_loc)

        exon_no = row["exon_no"]        
        if whitelist:
            var_count = whitelist.get_variants_per_exon(gene, exon_no)
            label = "%s.\n(%s)" % (exon_no, var_count)
            #
            #covered_exon_count = 0
            #for i in range(1,int(gexons.exon_no.max())+1):
            #    covered_exon_count += 1                    
            #    var_count = whitelist.get_variants_per_exon(gene, i)
            #    label = "%s.\n(%s)" % (covered_exon_count, var_count)
            
        else:
            label = "%s." % exon_no
        
        xlabels.append(label)
    
    plt.xticks(xlocs, xlabels)
    
    ############################################
    # set yticks
            
    locs,labels = plt.yticks()
    
    upper = int(ylim)
    if upper < 2000:
        stepsize = 500
    elif upper < 5000:
        stepsize = 1000
    else:
        stepsize = 2000    
    
    custom_locs =  range(1000, int(ylim), stepsize)
    custom_labels = [x-1000 for x in custom_locs]
    
    custom_locs_lower =  [-x for x in custom_locs]
    custom_labels_lower = custom_labels
    
    plt.yticks(custom_locs+custom_locs_lower, custom_labels+custom_labels_lower)    
    
    if save_per_gene:
        
        ############################################
        # save
        
        title = per_gene_output_folder + gene + " raw exon coverage"
        plt.savefig(title.replace(" ", "_")+".png", dpi=300)

def plot_mean_exon_coverage(gene, save_per_gene=False, df=None,
                           exons=None, plot_no=None, all_genes=None,
                           amplicons=None, per_gene_output_folder=None, whitelist=None):
        
        if save_per_gene:
            plt.figure(figsize=(20, 5))
        else:
            plt.subplot(len(all_genes), 1, plot_no)
            
        #try:
        plt.ylabel(gene)
        
        lower = 0
        upper = np.percentile(df[df.dp > 0].dp, 95)
        plt.ylim(lower, upper)
        
        ######################################################################
        
        gene= gene.upper()
        gexons = exons[(exons.gene_upper == gene)]
        
        gene_df = df[(df.gene.str.upper() == gene)]
        gene_df["dp_capped"] = gene_df.apply(lambda x : upper if x["dp"] > upper else x["dp"], axis=1)
        
        pdf = df[(df.gene.str.upper() == gene)].sort(columns = ["pos"]) 
                
        xs = []
        ys = [] 

        if len(gexons.exon_no.values) > 0:
            for i in range(1,int(gexons.exon_no.max())+1):
                if len(gexons[gexons.exon_no == i]) > 0: # exons have matching amplicons
         
                    #if whitelist:
                    #    logging.debug( "Exon %s, %s variants" % (i, whitelist.get_variants_per_exon(gene, i)))
                    
                    data = list(pdf[pdf.exon_no == i].dp)
                    mean_x = pdf[pdf.exon_no == i].pos.mean()
                    
                    if data and mean_x:
                        #print i, np.mean(data[0]), mean_x
                        ys.append(data)
                        xs.append(i)
                    else:
                        ys.append([0,0,0])
                        xs.append(i)

                    
                else: # empty exon that can't have coverage
                    midlevel = (lower + upper)/2.0
                    plt.text(i+0.25, 0, "X", fontsize=16)
                    
                    ys.append([0,0,0,0,0])
                    xs.append(i) 
        
        if len(xs) > 0:
            
            ######################################################################
            # plot bars
            
            plt.bar(xs, [np.mean(y) for y in ys], color=black)
            
            ######################################################################
            # create x ticks

            labels = [] 
            
            if whitelist:                                    
                covered_exon_count = 0
                for i in range(1,int(gexons.exon_no.max())+1):
                    covered_exon_count += 1                    
                    var_count = whitelist.get_variants_per_exon(gene, i)
                    label = "%s.\n(%s)" % (covered_exon_count, var_count)
                    labels.append(label)
            
            else:
                for i in range(1,int(gexons.exon_no.max())+1):
                    labels.append("%s." % i)
                    
            locs = [x + 0.5 for x in xs]
            
            plt.xticks(locs, labels)
        
        else: # empty placeholder
            plt.ylim(0,1)
            plt.xlim(0,1)
            plt.ylabel(gene)
            plt.xticks([],[])
            plt.yticks([],[])
            plt.text( 0.5, 0.5, 'No exons defined')
            return
                
        ############################################
        # set yticks
                
        locs,labels = plt.yticks()
        
        if upper < 2000:
            stepsize = 500
        elif upper < 5000:
            stepsize = 1000
        else:
            stepsize = 2000
        
        custom_locs =  range(0, int(upper), stepsize)
        custom_labels = custom_locs
        plt.yticks(custom_locs, custom_labels)
        
        ############################################
        # save
        
        if save_per_gene:            
            title = per_gene_output_folder + gene + " mean exon coverage"
            plt.savefig(title.replace(" ", "_")+".png", dpi=300)

##################################################################################################################
##################################################################################################################
    
def create_all_coverage_plots(filename, exons=None, exons_per_gene=None,
                       target_folder = None, whitelist=None, amplicons=None):
    
    
    sample = filename.split("/")[-1][:-13]
    header = ["chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]
    
    rawdf = pd.read_csv(filename, sep="\t", names=header)
    rawdf["pos"] = rawdf.start +rawdf.amplicon_pos
    rawdf["chrompos"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["pos"]) ]), axis = 1 ) 
    rawdf["name"] = rawdf.amplicon 
    rawdf["gene"] = rawdf.apply(lambda x : x["amplicon"].split("_")[0].upper(), axis = 1 ) 
    
    #######################################################################################
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
    df = df.fillna(0)
    df["gene"] = df.gene.str.upper()
    df["exon_no"] = df.apply(find_matching_exon_from_dict, axis=1, args=(exons_per_gene, ))    
        
    #######################################################################################    
    # setup variables

    all_genes = df.gene.unique()
    all_genes.sort()
    figure_length = len(all_genes)
    
    per_gene_output_folder = target_folder + "/gene_plots/"
    
    if not os.path.exists(per_gene_output_folder):
        os.makedirs(per_gene_output_folder)
    
    #######################################################################################
    # fancy plot for coverage vs actual chrom location
    
    plt.figure(figsize=(15, figure_length*2))
    logging.debug( "Creating overview plot of raw coverages." )
    
    plot_no = 0
    for gene in all_genes:
        plot_no += 1 
        plot_raw_exon_coverage(gene, df=df, exons=exons, all_genes=all_genes, amplicons=amplicons, per_gene_output_folder=per_gene_output_folder, whitelist=whitelist, 
                               plot_no=plot_no)
        
    plt.tight_layout()
    sample = filename.split("/")[-1].split(".")[0]
    title = target_folder + sample + " raw exon coverage"
    plt.savefig(title.replace(" ", "_")+".png", dpi=300)
    plt.close()
    
    ###################################
    # the same again, one plot per gene
    
    if output_per_gene_plots:
        logging.debug( "Creating per-gene plot of raw coverages." )
        for gene in all_genes:
            plot_raw_exon_coverage(gene, df=df, exons=exons, all_genes=all_genes, amplicons=amplicons, per_gene_output_folder=per_gene_output_folder, whitelist=whitelist, 
                                   save_per_gene = True)
            plt.close()
    
    #######################################################################################    
    # plot boxplots vs exon counts
     
    plt.figure(figsize=(15, figure_length*1.2))
    logging.debug( "Creating overview plot of mean coverages." )
   
    plot_no = 0
    for gene in all_genes:
        plot_no += 1 
        plot_mean_exon_coverage(gene, df=df, exons=exons, all_genes=all_genes, amplicons=amplicons, per_gene_output_folder=per_gene_output_folder, whitelist=whitelist, 
                               plot_no=plot_no)
        
    plt.tight_layout()
    
    sample = filename.split("/")[-1].split(".")[0]
    title = target_folder + sample + " summarized exon coverage"
    plt.savefig(title.replace(" ", "_")+".png", dpi=300)
    plt.close()

    ###################################
    # the same again, one plot per gene
    
    if output_per_gene_plots:
        logging.debug( "Creating per-gene plot of mean coverages." )    
        for gene in all_genes:
            plot_mean_exon_coverage(gene, df=df, exons=exons, all_genes=all_genes, amplicons=amplicons, per_gene_output_folder=per_gene_output_folder, whitelist=whitelist, 
                                   save_per_gene = True)
            plt.close()

    return df

##################################################################################################################
##################################################################################################################

if __name__ == "__main__":
    
    files = [x for x in os.listdir(".") if x.endswith("_coverage.tsv")]
    
    for filename in files:
        print "Starting " + filename
        create_all_coverage_plots(filename, target_folder = "./")


