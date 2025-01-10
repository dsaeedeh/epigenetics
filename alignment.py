import os
import argparse
import pandas as pd
from Bio import SeqIO
import pysam
from tqdm import tqdm

#http://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
#http://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

def parse_args():
    parser = argparse.ArgumentParser(description='Alignment pipeline')
    parser.add_argument('-i', '--input_dir', required=True, help='Path to QC-result directory')
    parser.add_argument('-r', '--reference_genome', required=False, help='Path to reference genome.fa file')
    parser.add_argument('-a', '--assembly_file', required=False, help='Path to assembly file')
    parser.add_argument('-g', '--gtf_file', required=False, help='Path to gtf file')
    parser.add_argument('-o', '--output_dir', required=True, help='Path to output directory')
    parser.add_argument('-m', '--model_name', required=True, help='model name', choices=['HISAT-2', 'STAR'])
    return parser.parse_args()

# Align all paired-end FASTQ files in a directory to a reference genome using HISAT2
def align_fastq(input_dir, output_dir, reference_genome, model_name):
    for sample_file_1 in tqdm(os.listdir(input_dir)):
        if sample_file_1.endswith('_1_val_1.fq.gz'):
            sample_name_prefix = sample_file_1.replace('_1_val_1.fq.gz', '')
            sample_file_2 = os.path.join(f"{sample_name_prefix}_2_val_2.fq.gz")
            sample_file_1 = os.path.join(input_dir, sample_file_1)
            sample_file_2 = os.path.join(input_dir, sample_file_2)
            if model_name == 'HISAT-2':
                genome = '/'.join(reference_genome.split('/')[:-1]) + '/genome'
                if not os.path.exists(f"{output_dir}/{sample_name_prefix}_sorted.bam"):
                    os.system(f"hisat2 -p 112 -x {genome} -1 {sample_file_1} -2 {sample_file_2} | samtools view -bSh >{output_dir}/{sample_name_prefix}.bam")
                    os.system(f"samtools sort -o {output_dir}/{sample_name_prefix}_sorted.bam {output_dir}/{sample_name_prefix}.bam")
                    os.system(f"samtools index {output_dir}/{sample_name_prefix}_sorted.bam")
                    os.system(f"rm {output_dir}/{sample_name_prefix}.bam")
                    print(f"Alignment of {sample_name_prefix} is done!")
                else:
                    print(f"Alignment of {sample_name_prefix} is already existed!")
            elif model_name == 'STAR':
                genome = "genome-STAR"
                os.system(f"STAR --runMode alignReads --genomeDir {genome} --outSAMtype BAM SortedByCoordinate --readFilesIn {sample_file_1} {sample_file_2} --runThreadN 20 --readFilesCommand zcat --outFileNamePrefix {output_dir}/{sample_name_prefix} ulimit -n 5000")

def remove_duplicates(input_dir):
    #read all bam files in input_dir and remove duplicates using picard and save it to input_dir with suffix _dedup.bam
    for file_name in tqdm(os.listdir(input_dir)):
        if file_name.endswith('.bam'):
            input_bam = os.path.join(input_dir, file_name)
            output_bam = os.path.join(input_dir, file_name.replace('.bam', '_dedup.bam'))
            os.system(f"java -jar picard/build/libs/picard.jar MarkDuplicates I={input_bam} O={output_bam} M={output_bam.replace('.bam', '.txt')} REMOVE_DUPLICATES=true")
            os.system(f"samtools index {output_bam}")
            #os.system(f"rm {input_bam}")
            print(f"remove duplicates of {file_name} is done!")

# Process all pairs of FASTQ files in the input directory and return a dictionary of summary statistics
def process_fastq(input_dir):

    file_name_list = []
    total_reads_list = []
    total_map_list = []
    unique_map_list = []
    multi_map_list = []
    read1_map_list = []
    read2_map_list = []
    positive_map_list = []
    negative_map_list = []
    splice_map_list = []
    unsplice_map_list = []
    proper_map_list = []

    # Replace 'your_file.bam' with the path to your BAM file
    for file_name in tqdm(os.listdir(input_dir)):
        if file_name.endswith('_dedup.bam'):
            input_bam = os.path.join(input_dir, file_name)
            # Initialize counters for each metric
            total_reads = 0
            total_map = 0
            unique_map = 0
            multi_map = 0
            read1_map = 0
            read2_map = 0
            positive_map = 0
            negative_map = 0
            splice_map = 0
            unsplice_map = 0
            proper_map = 0

            # Create a generator to filter reads in memory
            def filtered_reads_generator():
                with pysam.AlignmentFile(input_bam, 'rb') as bam_file:
                    for read in bam_file:
                        if read.flag & 256 == 0:  # Check the bitwise flag to exclude reads with flag 256
                            yield read

            # Iterate over the filtered reads using the generator
            for read in tqdm(filtered_reads_generator()):
                total_reads += 1
                total_map += 1 if not read.is_unmapped else 0
                try:
                    nh_tag = read.get_tag("NH")  # Try to get the NH tag
                    if nh_tag == 1:
                        unique_map += 1
                    elif nh_tag > 1:
                        multi_map += 1
                except KeyError:
                    # Handle the case where the NH tag is not present
                    pass

                if read.is_read1:
                    read1_map += 1
                if read.is_read2:
                    read2_map += 1
                if not read.is_unmapped:
                    if read.is_reverse:
                        negative_map += 1
                    else:
                        positive_map += 1
                    if "N" in read.cigarstring:
                        splice_map += 1
                    else:
                        unsplice_map += 1
                    if read.is_proper_pair:
                        proper_map += 1
            file_name = file_name.replace('Aligned.sortedByCoord.out_dedup.bam', '')
        
            # Append the results to the lists
            file_name_list.append(file_name)
            total_reads_list.append(total_reads)
            total_map_list.append(f"{total_map}({total_map / total_reads * 100:.2f}%)")
            unique_map_list.append(f"{unique_map}({unique_map / total_map * 100:.2f}%)")
            multi_map_list.append(f"{multi_map}({multi_map / total_map * 100:.2f}%)")
            read1_map_list.append(f"{read1_map}({read1_map / total_map * 100:.2f}%)")
            read2_map_list.append(f"{read2_map}({read2_map / total_map * 100:.2f}%)")
            positive_map_list.append(f"{positive_map}({positive_map / total_map * 100:.2f}%)")
            negative_map_list.append(f"{negative_map}({negative_map / total_map * 100:.2f}%)")
            splice_map_list.append(f"{splice_map}({splice_map / total_map * 100:.2f}%)")
            unsplice_map_list.append(f"{unsplice_map}({unsplice_map / total_map * 100:.2f}%)")
            proper_map_list.append(f"{proper_map}({proper_map / total_map * 100:.2f}%)")
                     
    # Create a DataFrame from the results
    result_df = pd.DataFrame({
        "Sample Name": file_name_list,
        "Total Reads": total_reads_list,
        "Total Mapped Reads": total_map_list,
        "Unique Mapped Reads": unique_map_list,
        "Multi-Mapped Reads": multi_map_list,
        "Read 1 Mapped": read1_map_list,
        "Read 2 Mapped": read2_map_list,
        "Positive Strand Mapped": positive_map_list,
        "Negative Strand Mapped": negative_map_list,
        "Splice Mapped Reads": splice_map_list,
        "Unsplice Mapped Reads": unsplice_map_list,
        "Properly Mapped Reads": proper_map_list
    })

    # Save the DataFrame to a CSV file
    result_df.to_csv(os.path.join(input_dir, 'alignment_results.csv'), index=False)
    print("alignment_results.csv is created!")


if __name__ == '__main__':
    args = parse_args()
    input_dir = args.input_dir
    reference_genome = args.reference_genome
    assembly_file = args.assembly_file
    gtf_file = args.gtf_file
    output_dir = args.output_dir
    model_name = args.model_name

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # if genome index does not exist, create it
    if model_name == 'HISAT-2' and not os.path.exists('/'.join(reference_genome.split('/')[:-1]) + '/genome.1.ht2'):
        print("Creating genome index...")
        os.system(f"hisat2-build {reference_genome} genome")
        print("genome index is created!")
    elif model_name == 'STAR' and not os.path.exists('genome-STAR'):
        print("Creating genome index...")
        os.system(f"STAR --runMode genomeGenerate --genomeDir genome-STAR --genomeFastaFiles {assembly_file} --sjdbGTFfile {gtf_file} --runThreadN 20")
        print("genome index is created!")

    print("Aligning FASTQ files...")
    #align_fastq(input_dir, output_dir, reference_genome, model_name)

    print("remove duplicates...")
    #remove_duplicates(output_dir)

    print("Processing FASTQ files...")
    process_fastq(output_dir)

    