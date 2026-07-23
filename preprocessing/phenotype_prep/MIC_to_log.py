import pandas as pd 
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="convert MIC into log")
parser.add_argument("-i","--input",required=True,help="file with Genome Id and MIC")
parser.add_argument("-o","--output", required=True,help="output file")
args=parser.parse_args()



df=pd.read_csv(args.input,sep=None,engine="python")

df["MIC"]=df['trait'].astype(float)
df["Phenotype"]=np.log2(df["MIC"])
df["Genome ID"]=df["taxon"]

df[["Genome ID","Phenotype"]].to_csv(args.output,sep="\t",index=False)
