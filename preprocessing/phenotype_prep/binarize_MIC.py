import pandas as pd
df =pd.read_csv("phenotype_srr.tsv",sep='\t') 
df=pd.DataFrame({ 
    "GenomeID":df["taxon"], 
    "Phenotype":(df["trait"]>1).astype(int)
    })
df.to_csv("phenotype_binary.csv",index=False)
