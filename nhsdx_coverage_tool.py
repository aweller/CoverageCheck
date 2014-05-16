
# coding: utf-8

# In[49]:

import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

project = "coverage_tool"
dp_cutoff = 50
pass_cutoff = 0.5
strand_diff_cutoff = 0.5


# In[50]:

#filename = "/home/andreas/vm/projects/nhs_coverage_tool/data/Alignment/ADM-00169_S3_coverage.tsv"
filename = "/home/andreas/datadump/201-2547_S10_coverage.tsv"
header = ["chr", "start", "stop", "amplicon", "na", "strand",  "pos", "dp"]

rawdf = pd.read_csv(filename, sep="\t", names=header)
rawdf["above_cutoff"] = rawdf.dp > dp_cutoff


# In[51]:

rawdf["name"] = rawdf.amplicon 
rawdf["amplicon"] = rawdf.apply(lambda x : "_".join([str(x["chr"]), str(x["start"]), str(x["stop"]) ]), axis = 1 ) 


# In[52]:

amplicon_sizes = dict(rawdf.pos.groupby(rawdf.amplicon).max())
def get_relative_pos(row):
    size = amplicon_sizes[row["amplicon"]]
    rel_size = row["pos"]/float(size)
    return round(rel_size,2)
    
rawdf["rel_pos"] = rawdf.apply( get_relative_pos, axis =1)


# In[76]:

covered_amplicons = len( rawdf[rawdf.dp > 0].amplicon.groupby(rawdf.amplicon) )
total_amplicons = len( rawdf.amplicon.groupby(rawdf.amplicon) )
covered_ratio = round( (covered_amplicons/float(total_amplicons)) * 100, 2)
covered_ratio


# In[80]:

plt.plot(rawdf.dp.groupby(rawdf.rel_pos).median())
plt.ylabel("Median coverage across amplicons")
plt.xlabel("Relative position across amplicon [%]")


# In[54]:

plt.plot(rawdf.dp.groupby(rawdf.pos).median())
plt.ylabel("Median coverage across amplicons")
plt.xlabel("Absolute position across amplicon [bp]")


# In[55]:

#rawdf[(rawdf.amplicon == "DNAJB6_+_DNAJB6_UserDefined_(23297017)_17919981") & (rawdf.pos == 2)]


# In[56]:

adf_dict= dict(name = rawdf.name.groupby([rawdf.amplicon, rawdf.strand]).max(),
               size = rawdf.pos.groupby([rawdf.amplicon, rawdf.strand]).max(),
               total_pass = rawdf[rawdf.above_cutoff].dp.groupby([rawdf.amplicon, rawdf.strand]).count(),
               median_dp = rawdf.dp.groupby([rawdf.amplicon, rawdf.strand]).median(),
               )

adf = pd.DataFrame(adf_dict) 
adf = adf.reset_index()
adf = adf.replace(NaN, 0)

adf["pass_ratio"] = adf.total_pass / adf.size.astype(float)


# In[57]:

adf.describe()


# In[59]:

df_dict= dict(name = adf.name.groupby(adf.amplicon).max(),
              size = adf.size.groupby(adf.amplicon).max(),
                total_pass = adf.total_pass.groupby(adf.amplicon).sum(),
               minus_pass_ratio = adf[adf.strand == "-"].pass_ratio.groupby(adf.amplicon).max(),
                plus_pass_ratio = adf[adf.strand == "+"].pass_ratio.groupby(adf.amplicon).max(),
                minus_median_dp = adf[adf.strand == "-"].median_dp.groupby(adf.amplicon).max(),
                plus_median_dp = adf[adf.strand == "+"].median_dp.groupby(adf.amplicon).max(),
               )

df = pd.DataFrame(df_dict) 
df = df.reset_index()
df.replace(NaN, 0)


df["total_pass_ratio"] =  df.total_pass / (df.size.astype(float) * 2.0)
#df["strand_pass_ratio"] =  abs(df.minus_pass_ratio - df.plus_pass_ratio)
df["strand_dp_ratio"] =  abs(1- ( df.minus_median_dp / df.plus_median_dp.astype(float) ))

df["pass_ratio_ok"] = df.total_pass_ratio > pass_cutoff
df["strand_ratio_ok"] = df.strand_dp_ratio < strand_diff_cutoff

df = df.replace(NaN, 0)
df = df.replace(inf, 1)

df.head()


# In[60]:

plt.hist(df[df.strand_dp_ratio < 1].strand_dp_ratio, bins = 40)
plt.xlabel("Difference between median coverage on left vs right strand ")
plt.ylabel("Amplicon count")


# In[61]:

plt.hist(df.total_pass_ratio, bins =20)
plt.xlabel("Ratio of bases above coverage cutoff")
plt.ylabel("Amplicon count")


# In[62]:

x = df.plus_median_dp
y = df.minus_median_dp

plt.scatter(x,y) 
plt.xlabel("Median coverage of bases on PLUS strand")
plt.ylabel("Median coverage of bases on MINUS strand")


# In[63]:

x = df.plus_pass_ratio
y = df.minus_pass_ratio

plt.scatter(x,y) 
plt.xlabel("Ratio of bases above cutoff on PLUS strand")
plt.ylabel("Ratio of bases above cutoff on MINUS strand")


# In[63]:



