
# coding: utf-8

# In[2]:

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import scipy.stats
from IPython.html.widgets import interact

get_ipython().magic(u'pylab inline')
get_ipython().magic(u'script')

project = "coverage_tool"
dp_cutoff = 50
pass_cutoff = 0.5
strand_diff_cutoff = 0.5


# In[55]:

#filename = "/home/andreas/vm/projects/nhs_coverage_tool/data/Alignment/ADM-00169_S3_coverage.tsv"
#filename = "/home/andreas/datadump/201-2547_S10_coverage.tsv"
filename = "/home/andreas/bioinfo/projects/nhs_coverage_tool/data/test/ARC-00031_S11_coverage.tsv"
header = ["chr", "start", "stop", "amplicon", "na", "strand",  "amplicon_pos", "dp"]

rawdf = pd.read_csv(filename, sep="\t", names=header)
rawdf["above_cutoff"] = rawdf.dp > dp_cutoff
rawdf["pos"] = rawdf.start +rawdf.amplicon_pos
rawdf["chrompos"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["pos"]) ]), axis = 1 ) 
rawdf["name"] = rawdf.amplicon 
#rawdf["amplicon"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["start"]), str(x["stop"]) ]), axis = 1 ) 
rawdf["gene"] = rawdf.apply(lambda x : x["amplicon"].split("_")[0], axis = 1 ) 


# In[56]:

_ = plt.hist(list(rawdf.dp), bins = 100)


# In[57]:

# df: per base

df = dict(chrom = rawdf.chr.groupby(rawdf.chrompos).min(), 
          pos = rawdf.pos.groupby(rawdf.chrompos).min(),
          gene = rawdf.gene.groupby(rawdf.chrompos).min(),
          minus_dp = rawdf[rawdf.strand == "-"].dp.groupby(rawdf.chrompos).max(),    
          plus_dp = rawdf[rawdf.strand == "+"].dp.groupby(rawdf.chrompos).max(),  
          dp = rawdf.dp.groupby(rawdf.chrompos).sum(),    
)

df = pd.DataFrame(df).reset_index()
df.head()


# In[58]:

def get_strand_ratios(row):
    minus = row["minus_dp"]
    plus = row["plus_dp"]
    
    if plus + minus == 0:
        return 0
    if plus == 0:
        return minus
    elif minus == 0:
        return plus
    else:
        return minus/float(plus)
    
df["strand_ratio"] = df.apply(get_strand_ratios, axis =1 )


# In[53]:

max_strand_ratio = 5

failed = df[(df.strand_ratio > max_strand_ratio) & (df.dp > 10)]


# In[38]:

_ = plt.hist(list(df.strand_ratio), bins = 50)


# In[40]:

x = df.minus_dp
y = df.plus_dp

plot = sns.jointplot(x,y, xlim=(-1000,10000), ylim=(-1000,10000))


# In[ ]:



