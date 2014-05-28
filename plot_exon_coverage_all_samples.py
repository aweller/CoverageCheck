import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import scipy.stats
from IPython.html.widgets import interact

def plot_exon_coverage(filename, target_folder = "./"):
    
    header = ["sample", "chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]
    
    df = pd.read_csv(target_folder + filename, sep="\t", names=header)
    df["pos"] = df.start +df.amplicon_pos
    df["chrompos"] = df.apply(lambda x : "_".join([str(x["chr"]), str(x["pos"]) ]), axis = 1 ) 
    df["name"] = df.amplicon 
    df["gene"] = df.apply(lambda x : x["amplicon"].split("_")[0], axis = 1 ) 
    
    ##################################################################################################
    
    exon_filename = "/home/andreas/bioinfo/core/general/data/HumanExons_Ensembl_v65_merged.tsv"
    header = ["exon_start", "exon_stop", "gene", "exon_no"]
    exons = pd.read_csv(exon_filename, sep="\t", names=header)
    
    exons["gene_upper"] = exons.gene.str.upper()
    exons = exons.sort(columns = ["gene", "exon_start", "exon_stop"])
    
    ##################################################################################################
    
    exon_summary = pd.DataFrame(dict(mean_dp = df.dp.groupby(df.gene).mean().round().astype(int), median_dp = df.dp.groupby(df.gene).median())).reset_index()
    title = "all samples exon summary"
    exon_summary.to_csv(target_folder + title.replace(" ", "_")+".tsv", sep = "\t", index=False)
    
    ##################################################################################################
    
    exon_dict = {}
    
    for _, row in exons.iterrows():
        gene = row["gene"]
        start, stop = int(row["exon_start"]), int(row["exon_stop"])
        exon_no = row["exon_no"]
        
        if not exon_dict.get(gene):
            exon_dict[gene] = []
            
        exon_dict[gene].append((start, stop, exon_no))
    
    
    def find_matching_exon_from_dict(row):    
        gene = row["gene"]
        pos = row["pos"]
    
        exon_no = 0
        try:
            for region in exon_dict[gene]:
                if region[0] < pos < region[1]:
                    exon_no = region[2]
                    break
        except:
            pass
        
        return exon_no     
        
    df["exon_no"] = df.apply(find_matching_exon_from_dict, axis = 1)
    
    
    ##################################################################################################
    
    def plot_exon(row):
        start = int(row["exon_start"])
        stop = int(row["exon_stop"])
        size = stop - start    
        
        #print start, stop
        
        rectangle = plt.Rectangle((start, -20), size, 10, fc='red')
        plt.gca().add_patch(rectangle)
    
    
    ##################################################################################################
    # plot a boxplot for samples per exon
    
    all_genes = df.gene.unique()
    all_genes.sort()
    #all_genes = ["MIR147", ]
    
    plot_no = 0
    
    plt.figure(figsize=(15, 40))
    
    upper = np.percentile(df[df.dp > 0].dp, 95)
    lower = -np.percentile(df[df.dp > 0].dp, 40)
    
    for gene in all_genes:
        plot_no += 1
        plt.subplot(len(all_genes), 1, plot_no)
    
        plt.ylabel(gene)
        plt.axhline(0)
        plt.ylim(lower, upper)
    
        ######################################################################
    
        gene = gene.upper()
        gene_exons = exons[(exons.gene_upper == gene)]
        gene_df = df[(df.gene.str.upper() == gene)].sort(columns = ["pos"]) 
    
        xs = []
        ys = [] 
        wids = []
    
        for i in range(1,int(gene_df.exon_no.max())+1):
    
            data = list(gene_df[gene_df.exon_no == i].dp.groupby(gene_df.sample).mean())
            mean_x = gene_df[gene_df.exon_no == i].pos.mean()
    
            if data and mean_x:
                #print gene, i, np.mean(data[0]), mean_x
                ys.append(data)
                xs.append(i)
    
            else:
                ys.append([0,0.1,0])
                xs.append(i)
    
        if len(xs) + len(ys) > 1:
            plt.boxplot(ys, positions=xs)
        ######################################################################
    
    
    plt.tight_layout()
    
    title = "all samples summarized exon coverage boxplots"
    plt.savefig(target_folder + title.replace(" ", "_")+".png", dpi=300)
    plt.close()
    
    ##################################################################################################
    # plot a line with the mean per exon per sample
    
    all_genes = df.gene.unique()
    all_genes.sort()
    #all_genes = ["MIR147", ]
    
    all_samples = df.sample.unique()
    
    plot_no = 0
    
    plt.figure(figsize=(15, 40))
    
    upper = np.percentile(df[df.dp > 0].dp, 95)
    lower = -np.percentile(df[df.dp > 0].dp, 40)
    
    for gene in all_genes:
        plot_no += 1
        plt.subplot(len(all_genes), 1, plot_no)
    
        plt.ylabel(gene)
        plt.axhline(0, linewidth=1, color='black')
        plt.ylim(lower, upper)
    
        ######################################################################
    
        gene = gene.upper()
        gene_df = df[(df.gene.str.upper() == gene) & (df.exon_no != 0)].sort(columns = ["pos"]) 
        
        if len(gene_df) > 0:
            for sample in all_samples:
                sample_df = gene_df[gene_df.sample == sample]
                x = []
                y = []
    
                for i in range(1,int(gene_df.exon_no.max())+1):
                    try:
                        mean_dp = int(sample_df[sample_df.exon_no == i].dp.median())
    
                    except:
                        mean_dp = 0
                    y.append(mean_dp)
    
                x = range(1, len(y)+1)
    
                plt.plot(x, y)
            
        plt.xlim(0, gene_df.exon_no.max()+1)
    
        ######################################################################
    
    plt.tight_layout()
    
    title = "all samples summarized exon coverage lineplots"
    plt.savefig(target_folder + title.replace(" ", "_")+".png", dpi=300)
    plt.close()
