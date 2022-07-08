"""Main module."""
import sys
import os
import pysam
import numpy as np
import pandas as pd
__version__ = 'v1.0'
   
def main():
    args = parse_args()
    genome_ref, annotation_ref = config_getter('config.ini')
    try:
        bamfile = pysam.AlignmentFile(args.input, 'rb', require_index = True)
    except FileNotFoundError:
        try:
            print(' bam index file is needed')
        except FileNotFoundError:
            print(f'There is a problem opening bam file at: {args.input}')
    
    chrms = ['chr1','chr2','chr3','chr4', 'chr5',
         'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11','chr12', 'chr13', 'chr14', 'chr15',
         'chr16','chr17', 'chr18', 'chr19', 'chr20',
         'chr21', 'chr22', 'chrX','chrY', 'chrM']

    ce_total=[]
    for chrm in chrms:
        print(f'Finding cryptic exon in {chrm}')
        sys.stdout.flush()
        ce  = ce_caller(bamfile,
                    annotation_ref,
                    chrm,
                    args.stranded,
                    args.mapq)
        ce_total.extend(ce)

    ce_filtered=filter_ce(ce_total, ao_min=1, PSI_min=0)
        
    out_file_name = args.out
    if not out_file_name:
        prefix = os.path.splitext(os.path.basename(bamfile.filename.decode('UTF-8')))[0]
        out_file_name = f'{prefix}.ce'
    print(f'Finished ce calling and filtering. Printing to {out_file_name}')
 
    df = pd.Series(ce_filtered).reset_index()   
    df.columns = ['chrom','D','A','start','end','ao1','ao2','ao3','a_count','PSI','strand','gene_id','gene_name','p_uniform']   
    df.sort_values(by=['chrom','D'],inplace=True)
    df.to_csv(out_file_name,index=False,sep='\t')
 
            

if __name__ == '__main__':

    main()
    
