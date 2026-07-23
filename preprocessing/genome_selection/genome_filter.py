import pandas as pd
import argparse
import os

def filter(df,mic_col,acc_col,sus_thresh,res_thresh,n=None):
    df=df[[acc_col,mic_col]].dropna()
    df[mic_col]=pd.to_numeric(df[mic_col],errors='coerce')
    
    sus=df[df[mic_col] <= sus_thresh]\
        .sort_values(by=mic_col, ascending=True)\
        .drop_duplicates(subset=acc_col)
    
    res=df[df[mic_col] >= res_thresh]\
        .sort_values(by=mic_col,ascending=True)\
        .drop_duplicates(subset=acc_col)
    
    if n:
        sus=sus.head(n)
        res=res.head(n)
    
    return sus[acc_col],res[acc_col]

def main():
    parser=argparse.ArgumentParser(description="select sus and res genomes based on MIC values")
    parser.add_argument('--file','-f', required=True)
    parser.add_argument('--sus_threshold','-s', type=float,default=4,help="susceptibility threshold")
    parser.add_argument('--res_threshold','-r', type=float,default=16,help="resistance threshold")
    parser.add_argument('--number_of_genomes','-n',type=int, help='optional:if not specified,this will take all res and sus genomes')
    parser.add_argument('--output_dir','-o',default='.')
    
    args=parser.parse_args()
    df=pd.read_excel(args.file ,engine='openpyxl')
    
    sus,res=filter(df, 'Measurement Value', 'Genome ID', args.sus_threshold,args.res_threshold,args.number_of_genomes)
    
    sus.to_csv(os.path.join(args.output_dir, 'genome_Sus_list.txt'),index=False,header=False)
    res.to_csv(os.path.join(args.output_dir,'genome_Res_list.txt'),index=False,header=False)
    
    print(f'written {len(sus)} susceptible and {len(res)} resistant genomes IDs')
    

if __name__=='__main__':
    main()