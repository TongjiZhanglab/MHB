#!/usr/bin/python

# Time-stamp: <2020-01-08 Zhaowei Yu>
"""Description: Locate the location of the motif in the genome
@author: Zhaowei Yu
@contact: edwardbioyu@gmail.com
"""

# -----------
#  modules
# -----------

import os,sys
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing

# ----------------
#  misc functions
# ----------------

def ave_methyl(group):
        """ methylation level = methylC / totalC """
        chrom = group["chrom"].values[0]
        start = group["start"].values[0]
        end = group["end"].values[0]
        dmr = group["dmr"].values[0]
        totalC = group["totalC"].sum()
        methylC = group["methylC"].sum()
        ratio = float(methylC) / float(totalC)
        ser = pd.Series([chrom,start,end,dmr,ratio])
        return(ser)

def applyParallel(dfGrouped, func):
        retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group) for name, group in dfGrouped)
        return pd.concat(retLst)

# --------
#  main
# --------

inf = sys.argv[1]
outf = sys.argv[2]

df_in = pd.read_table(inf, header = None, sep = "\t")
df_raw = df_in[df_in.columns[0:10]]     # cut the top 10 columns
df_raw.columns = ["chrom","start","end","dmr","chrom2","start2","end2","ratio","totalC","methylC"]
df_raw[len(df_raw.columns)] = df_raw["chrom"] + "_" + df_raw["start"].astype("str") + "_" + df_raw["end"].astype("str") + "_" + df_raw["dmr"]
df_raw.columns = ["chrom","start","end","dmr","chrom2","start2","end2","ratio","totalC","methylC","label"]
df_filter = df_raw[df_raw["totalC"]>=3]
grouped = df_filter.groupby("label", sort = False)      # grouped rows by dmr name
df_out = grouped.apply(ave_methyl)
# df_out = applyParallel(grouped, ave_methyl)
df_out.to_csv(outf, sep = "\t",  index = False, header = False)
