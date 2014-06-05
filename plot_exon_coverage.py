import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import scipy.stats
from IPython.html.widgets import interact
import os
import logging

def plot_exon_coverage(filename, exons=None, exons_per_gene = None, target_folder = None, whitelist=None):
    
    sample = filename.split("/")[-1][:-13]
    header = ["chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]
    
    rawdf = pd.read_csv(filename, sep="\t", names=header)
    rawdf["pos"] = rawdf.start +rawdf.amplicon_pos
    rawdf["chrompos"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["pos"]) ]), axis = 1 ) 
    rawdf["name"] = rawdf.amplicon 
    rawdf["gene"] = rawdf.apply(lambda x : x["amplicon"].split("_")[0].upper(), axis = 1 ) 
    
    ########################################################################################
    ## parse the complete list of human exons
    #
    ##exon_filename = "/home/andreas/bioinfo/core/general/data/HumanExons_Ensembl_v65_merged.tsv"
    #exon_filename = "/home/andreas/bioinfo/core/general/data/HumanExons_Ensembl_v75_merged.tsv"
    #
    #header = ["chrom", "exon_start", "exon_stop", "gene", "strand", "exon_no"]
    #exons = pd.read_csv(exon_filename, sep="\t", names=header)
    #exons["gene_upper"] = exons.gene.str.upper()
    #exons = exons.sort(columns = ["gene", "exon_start", "exon_stop"])
    #    
    #exons_per_gene = {}
    #for _, row in exons.iterrows():
    #    gene = row["gene"].upper()
    #    start, stop = int(row["exon_start"]), int(row["exon_stop"])
    #    exon_no = row["exon_no"]
    #    
    #    if not exons_per_gene.get(gene):
    #        exons_per_gene[gene] = []
    #        
    #    exons_per_gene[gene].append((start, stop, exon_no))
    
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
    df["gene"] = df.gene.str.upper()
    
    #######################################################################################
    
    def find_matching_exon_from_dict(row):    
        gene = row["gene"]
        pos = row["pos"]
    
        exon_no = 0
        try:
            for region in exons_per_gene[gene]:
                if region[0] < pos < region[1]:
                    exon_no = region[2]
                    break
        except:
            pass
        
        return exon_no     
        
    df["exon_no"] = df.apply(find_matching_exon_from_dict, axis = 1)
    
    #######################################################################################
    
    def plot_exon(row):
        start = int(row["exon_start"])
        stop = int(row["exon_stop"])
        size = stop - start    
        
        #print start, stop
        
        rectangle = plt.Rectangle((start, -20), size, 10, fc='red')
        plt.gca().add_patch(rectangle)
    
    
    #######################################################################################    
    # plot individual bases vs chromosome locations

    all_genes = df.gene.unique()
    all_genes.sort()
    
    plot_no = 0
    
    figure_length = len(all_genes)
    plt.figure(figsize=(15, figure_length))
    
    upper = np.percentile(df[df.dp > 0].dp, 90)
    lower = -np.percentile(df[df.dp > 0].dp, 40)
    plt.ylim(lower, upper)
        
    for gene in all_genes:
        plot_no += 1
        plt.subplot(len(all_genes), 1, plot_no)
        
        try:
            plt.ylabel(gene)
            gene= gene.upper()
            
            ######################################################################
        
            gene_exons = exons[(exons.gene_upper == gene)]
            
            gene_exons.apply(plot_exon, axis =1 )
        
            x_start = int(gene_exons.head(1).exon_start) - 1000
            x_stop = int(gene_exons.tail(1).exon_start) + 1000 
            
            plt.xlim(x_start, x_stop)
            plt.ylim(-20,210)
            
            ######################################################################
        
            gene_df = df[(df.gene.str.upper() == gene)]
            gene_df["dp_capped"] = gene_df.apply(lambda x : 200 if x["dp"] > 200 else x["dp"], axis=1)
            
            pdf = gene_df[gene_df.dp > 10].sort(columns = ["pos"])
            y = pdf.dp_capped
            x = pdf.pos
                
            plt.scatter(x,y, c="black", s=10)
            #plt.gray()
            
            plt.axhline(0)
            
            ######################################################################
            # create x labels
            
            locs = []
            labels = [] 
            
            for i in range(1,int(gene_exons.exon_no.max())+1):                
                if len(gene_exons[gene_exons.exon_no == i]) > 0: # exons have matching amplicons
                    start = gene_exons[gene_exons.exon_no == i].exon_start.min()
                    stop = gene_exons[gene_exons.exon_no == i].exon_stop.max()
                    
                    locs.append( np.mean([start, stop]))
                    
                    label = str(i)
                    #if i % 2 == 0:
                    #    label = "\n" + str(i)
                    #else:
                    #    label = str(i) + "\n"
                    
                    labels.append(label)
                    
            plt.xticks(locs, labels)
        
        except:
            logging.warning( "Gene %s in sample %s is not plotted due to an error. Is it in the exon list?" % (gene, sample))
            
            
    plt.tight_layout()
    
    sample = filename.split("/")[-1].split(".")[0]
    title = target_folder + sample + " raw exon coverage"
    plt.savefig(title.replace(" ", "_")+".png", dpi=300)
    plt.close()
    
    #######################################################################################    
    # plot boxplots vs exon counts
    
    all_genes = df.gene.unique()
    all_genes.sort()
    #all_genes = ["EFNA5", ]
    
    plot_no = 0
    
    figure_length = len(all_genes)
    plt.figure(figsize=(15, figure_length))
    
    upper = np.percentile(df[df.dp > 0].dp, 90)
    lower = -np.percentile(df[df.dp > 0].dp, 40)
    
    for gene in all_genes:
        plot_no += 1
        plt.subplot(len(all_genes), 1, plot_no)
        
        #try:
        plt.ylabel(gene)
        plt.axhline(0)
        plt.ylim(lower, upper)
        
        ######################################################################
        
        gene= gene.upper()
        gene_exons = exons[(exons.gene_upper == gene)]
        
        gene_df = df[(df.gene.str.upper() == gene)]
        gene_df["dp_capped"] = gene_df.apply(lambda x : upper if x["dp"] > upper else x["dp"], axis=1)
        
        pdf = df[(df.gene.str.upper() == gene)].sort(columns = ["pos"]) 
        
        xs = []
        ys = [] 
        
        logging.debug( "-" * 150 )
        logging.debug( gene )

        if len(gene_exons.exon_no.values) > 0:
            for i in range(1,int(gene_exons.exon_no.max())+1):
                if len(gene_exons[gene_exons.exon_no == i]) > 0: # exons have matching amplicons
         
                    if whitelist:
                        logging.debug( "Exon %s, %s variants" % (i, whitelist.get_variants_per_exon(gene, i)))
                    
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
                    midlevel = lower + (upper-lower)/2
                    plt.text(i, midlevel, "X", fontsize=16)
                    
                    ys.append([0,0,0])
                    xs.append(i)
        
        if len(xs) > 0:
            
            plt.boxplot(ys, positions=xs)
            
            if whitelist:
                    
                ######################################################################
                # create x labels
                        
                locs = []
                labels = [] 
            
                covered_exon_count = 0
                for i in range(1,int(gene_exons.exon_no.max())+1):
                    #if len(gene_exons[gene_exons.exon_no == i]) > 0: # exons have matching amplicons
                    covered_exon_count += 1
                    locs.append( i )
                    
                    var_count = whitelist.get_variants_per_exon(gene, i)
                    label = "%s.\n(%s)" % (covered_exon_count, var_count)
                    labels.append(label)
                    #else: # empty exon that can't have coverage
                    #    locs.append( i )
                    #    label = ""
                    #    labels.append(label)
                    
                plt.xticks(locs, labels)
        
    plt.tight_layout()
    
    sample = filename.split("/")[-1].split(".")[0]
    title = target_folder + sample + " summarized exon coverage"
    plt.savefig(title.replace(" ", "_")+".png", dpi=300)
    plt.close()


####################################################################################

if __name__ == "__main__":
    
    files = [x for x in os.listdir(".") if x.endswith("_coverage.tsv")]
    
    for filename in files:
        print "Starting " + filename
        plot_exon_coverage(filename, target_folder = "./")


