import os
import gzip
import argparse
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


# Define a function to calculate the size in GB
def bytes_to_gb(byte_count):
    return byte_count / (1000 ** 3)

# Calculate the GC content percentage for a sequence
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count

# Create a function to process a pair of samples
def process_sample_pair(sample_path_1, sample_path_2, output_dir):
    sample_name = os.path.basename(sample_path_1).replace('_1.fastq.gz', '')
    raw_read_count = 0
    raw_base_count = 0
    clean_read_count = 0
    clean_base_count = 0
    phred_20_count = 0
    phred_30_count = 0
    gc_percentage = 0  # Initialize GC percentage
    clean_reads = []

    for sample_path in [sample_path_1, sample_path_2]:
        with gzip.open(sample_path, 'rt') as fastq_file:
            for record in SeqIO.parse(fastq_file, 'fastq'):
                raw_read_count += 1
                raw_base_count += len(record.seq)
    
    sample_path_1 = os.path.join(output_dir, sample_name + '_1_val_1.fq.gz')
    sample_path_2 = os.path.join(output_dir, sample_name + '_2_val_2.fq.gz')
    
    for sample_path in [sample_path_1, sample_path_2]:
        with gzip.open(sample_path, 'rt') as fastq_file:
            for record in SeqIO.parse(fastq_file, 'fastq'):
                clean_read_count += 1
                clean_base_count += len(record.seq)
                clean_reads.append(record)

                # Calculate Phred scores
                phred_20_count += sum(1 for q in record.letter_annotations["phred_quality"] if q > 20)
                phred_30_count += sum(1 for q in record.letter_annotations["phred_quality"] if q > 30)

                # Calculate GC content percentage
                gc_percentage += calculate_gc_content(record.seq)

    # Calculate error rate
    error_rate = 1 - (clean_read_count / raw_read_count)

    # Calculate percentages
    phred_20_percent = (phred_20_count / clean_base_count) * 100
    phred_30_percent = (phred_30_count / clean_base_count) * 100
    gc_percentage = (gc_percentage / (clean_base_count + 1)) * 100  # +1 to avoid division by zero

    # Create a dictionary to store the results
    result_dict = [sample_name, raw_read_count,
        bytes_to_gb(raw_base_count),
        clean_read_count,
        bytes_to_gb(clean_base_count),
        error_rate,
        phred_20_percent,
        phred_30_percent,
        gc_percentage
    ]

    return result_dict

def process_fastq(input_dir, output_dir):
    sample_name_list = []
    raw_read_count_list = []
    raw_base_size_list = []
    clean_read_count_list = []
    clean_base_size_list = []
    error_rate_list = []
    phred_20_percent_list = []
    phred_30_percent_list = []
    gc_content_percent_list = []
    # run trimgalore for each pair of samples in the input directory and save the results in the output directory
    for sample_file_1 in tqdm(os.listdir(input_dir)):
        if sample_file_1.endswith('_1.fastq.gz'):
            sample_name_prefix = sample_file_1.replace('_1.fastq.gz', '')
            sample_file_2 = os.path.join(f"{sample_name_prefix}_2.fastq.gz")
            os.system(f"trim_galore --cores 64 --trim-n --max_n 15 --paired {sample_file_1} {sample_file_2} -o {output_dir}")
            print(f"Trimming of {sample_name_prefix} is done!")

    # Get a list of _1.fastq.gz files in the input directory
    sample_files_1 = [f for f in os.listdir(input_dir) if f.endswith('_1.fastq.gz')]

    # Process each pair of samples and append the results to the DataFrame and CSV file
    for sample_file_1 in tqdm(sample_files_1):
        sample_file_1 = os.path.join(input_dir, sample_file_1)
        sample_name_prefix = sample_file_1.replace('_1.fastq.gz', '')
        sample_file_2 = os.path.join(f"{sample_name_prefix}_2.fastq.gz")
        sample_name = os.path.basename(sample_file_1).replace('_1.fastq.gz', '')
        result_dict = process_sample_pair(sample_file_1, sample_file_2, output_dir)
        sample_name_list.append(result_dict[0])
        raw_read_count_list.append(result_dict[1])
        raw_base_size_list.append(result_dict[2])
        clean_read_count_list.append(result_dict[3])
        clean_base_size_list.append(result_dict[4])
        error_rate_list.append(result_dict[5])
        phred_20_percent_list.append(result_dict[6])
        phred_30_percent_list.append(result_dict[7])
        gc_content_percent_list.append(result_dict[8])

    # Create an empty DataFrame to store the results
    result_df = pd.DataFrame({
        "Sample Name": sample_name_list,
        "Raw Read Count": raw_read_count_list,
        "Raw Base Size (GB)": raw_base_size_list,
        "Clean Read Count": clean_read_count_list,
        "Clean Base Size (GB)": clean_base_size_list,
        "Error Rate": error_rate_list,
        "Phred > 20 (%)": phred_20_percent_list,
        "Phred > 30 (%)": phred_30_percent_list,
        "GC Content (%)": gc_content_percent_list
    })

    result_df.to_csv(os.path.join(output_dir, 'sample_results.csv'), index=False)
    print("sample_results.csv is created!")


# main function
if __name__ == 'main':
    parser = argparse.ArgumentParser(description="Process FASTQ files and append summary statistics to a CSV file.")
    parser.add_argument("input_dir", help="Input directory containing Raw FASTQ files.")
    parser.add_argument("output_dir", help="Output directory for saving clean FASTQ files.")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    
    # if output directory does not exist, create it
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    
    print("Processing FASTQ files...")
    process_fastq(input_dir, output_dir)

    
