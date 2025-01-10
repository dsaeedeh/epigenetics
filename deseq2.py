
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
import pandas as pd
from glob import glob
import os
import argparse

pd.options.display.float_format = '{:.2f}'.format

def parse_args():
    parser = argparse.ArgumentParser(description='Differential expression analysis using DESeq2')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory of bam files')
    parser.add_argument('-g', '--gtf_file', required=True, help='GTF file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    args = parser.parse_args()
    return args

# Count features of all bam files in input directory using featureCounts
def count_features(input_dir, gtf_file, output_dir):
    # Mapped regions can be classified as exons, introns, or intergenic regions. 
    # The program featureCounts can be used to count the number of reads that map to each of these regions.
    # -T: number of threads
    # -a: input annotation file (gtf)
    # -o: output file name
    # input_dir: input directory of bam files
    # *_dedup.bam: deduplicated bam files

    os.system(f"featureCounts -T 8 -p -t exon -g gene_id -a {gtf_file} -o {output_dir}/exon_count.out {input_dir}/*_dedup.bam")
    os.system(f"featureCounts -T 16 -p -a {gtf_file} -o {output_dir}/gene_count.out {input_dir}/*_dedup.bam")
 
if __name__ == '__main__':
    args = parse_args()
    input_dir = args.input_dir
    gtf_file = args.gtf_file
    output_dir = args.output_dir

    count_features(input_dir, gtf_file, output_dir)
    print("Counting features is done!")
