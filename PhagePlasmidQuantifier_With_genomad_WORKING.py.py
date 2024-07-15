import pandas as pd
import sys
import subprocess
import re
import os
import logging

# Set up logging
logging.basicConfig(filename='script_errors.log', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')

# Column names for HMMsearch output
column_names = [
    "hmm_target",
    "hmm_target_accession",
    "query_protein",
    "query_protein_accession",
    "full_evalue",
    "full_score",
    "full_bias",
    "domain_evalue",
    "domain_score",
    "domain_bias",
    "exp",
    "reg",
    "clu",
    "ov",
    "env",
    "dom",
    "rep",
    "inc",
    "description",
]

column_names_container = {
    "hmm_target": 0,
    "hmm_target_accession": 1,
    "query_protein": 2,
    "query_protein_accession": 3,
    "full_evalue": 4,
    "full_score": 5,
    "full_bias": 6,
    "domain_evalue": 7,
    "domain_score": 8,
    "domain_bias": 9,
    "exp": 10,
    "reg": 11,
    "clu": 12,
    "ov": 13,
    "env": 14,
    "dom": 15,
    "rep": 16,
    "inc": 17,
    "description": 18,
}

def extract_columns(line):
    line = line.rstrip()
    return re.split(r"\s+", line, maxsplit=18)

def parse_hmmtbl_line(line):
    columns = extract_columns(line)
    parsed_data = [
        columns[column_names_container[column_name] - 1]
        for column_name in column_names_container.keys()
    ]
    return parsed_data

def parse_hmmtbl(file_path):
    parsed_data_list = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parsed_line = parse_hmmtbl_line(line)
                parsed_data_list.append(parsed_line)
    return parsed_data_list

def write_to_tbl(parsed_data, output_file):
    with open(output_file, 'w') as file:
        header = "\t".join(column_names_container.keys()) + "\n"
        file.write(header)
        for row in parsed_data:
            row_data = "\t".join(row) + "\n"
            file.write(row_data)

def extract_seq_name(hmm_target_accession):
    return '_'.join(hmm_target_accession.split('_')[:2])

def process_hmmsearch(input_file):
    print(f"Processing HMMsearch output from {input_file}...")
    parsed_data = parse_hmmtbl(input_file)
    df = pd.DataFrame(parsed_data, columns=column_names)
    df['seq_name'] = df['hmm_target_accession'].apply(extract_seq_name)
    filtered_df = df[['seq_name', 'query_protein_accession']]
    unique_scaffolds = filtered_df.drop_duplicates()
    unique_counts = unique_scaffolds.groupby('seq_name')['query_protein_accession'].nunique()
    result_df = unique_counts.reset_index()
    result_df.columns = ['seq_name', 'unique_protein_hits']
    print("HMMsearch processing completed.")
    return result_df

def merge_csv(input_csv1, input_csv2, output_csv):
    print(f"Merging {input_csv1} and {input_csv2} into {output_csv}...")
    df1 = pd.read_csv(input_csv1, sep='\t')
    df2 = pd.read_csv(input_csv2)
    merged_df = pd.merge(df1, df2, on='seq_name')
    merged_df['proteins_over_genes'] = merged_df['unique_protein_hits'] / merged_df['n_genes']
    output_df = merged_df[['seq_name', 'n_genes', 'unique_protein_hits', 'proteins_over_genes']]
    output_df.to_csv(output_csv, index=False)
    print(f"Output saved to {output_csv}.")

def run_hmmsearch_and_process(genomad_output_fasta, hmm_folder, prefix):
    print(f"Running HMMsearch for {prefix}...")

    # Check if the fasta file is empty
    if os.path.getsize(genomad_output_fasta) == 0:
        print(f"Skipping HMMsearch for {prefix} because {genomad_output_fasta} is empty.")
        return None

    # Iterate over HMM files in the specified folder
    for hmm_file in os.listdir(hmm_folder):
        if hmm_file.endswith('.hmm'):
            hmm_path = os.path.join(hmm_folder, hmm_file)
            if os.path.getsize(hmm_path) == 0:
                print(f"Skipping empty HMM file: {hmm_path}")
                continue

            hmmsearch_output = f'{prefix}_{os.path.splitext(hmm_file)[0]}_hmmsearch_output.tblout'
            try:
                subprocess.run(['hmmsearch', '--tblout', hmmsearch_output, hmm_path, genomad_output_fasta], check=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running HMMsearch for {prefix} with {hmm_file}: {str(e)}")
                continue  # Skip this HMM and continue with the next one
            print(f"HMMsearch completed successfully for {hmm_file}.")

            hmmsearch_results = process_hmmsearch(hmmsearch_output)
            temp_file = f'{prefix}_{os.path.splitext(hmm_file)[0]}_PhageToPlasmidHits.csv'
            hmmsearch_results.to_csv(temp_file, index=False)
            print(f"Temporary file saved to {temp_file}.")
    return temp_file

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 PhagePlasmidQuantifier.py genomad_input.fasta genomad_db virus_hmm_folder plasmid_hmm_folder splits threads output_folder")
        sys.exit(1)

    genomad_input = sys.argv[1]
    genomad_db = sys.argv[2]
    virus_hmm_folder = sys.argv[3]
    plasmid_hmm_folder = sys.argv[4]
    splits = sys.argv[5]
    threads = sys.argv[6]
    output_folder = sys.argv[7]

    os.makedirs(output_folder, exist_ok=True)

    fasta_basename = os.path.splitext(os.path.basename(genomad_input))[0]
    genomad_output_dir = os.path.join(output_folder, 'genomad_output')

    virus_output_fasta = os.path.join(genomad_output_dir, f'{fasta_basename}_summary/{fasta_basename}_virus_proteins.faa')
    plasmid_output_fasta = os.path.join(genomad_output_dir, f'{fasta_basename}_summary/{fasta_basename}_plasmid_proteins.faa')
    virus_output_csv = os.path.join(genomad_output_dir, f'{fasta_basename}_summary/{fasta_basename}_virus_summary.tsv')
    plasmid_output_csv = os.path.join(genomad_output_dir, f'{fasta_basename}_summary/{fasta_basename}_plasmid_summary.tsv')

    print("Running Genomad...")
    try:
        subprocess.run(['genomad', 'end-to-end', '--splits', splits, '--threads', threads, genomad_input, genomad_output_dir, genomad_db], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running Genomad: {str(e)}")
        sys.exit(1)
    print("Genomad completed successfully.")

    for file_path in [virus_output_fasta, plasmid_output_fasta, virus_output_csv, plasmid_output_csv]:
        if not os.path.exists(file_path):
            logging.error(f"Genomad output file {file_path} not found. Please check if Genomad ran successfully.")
            sys.exit(1)

    virus_temp_file = run_hmmsearch_and_process(virus_output_fasta, virus_hmm_folder, 'virus')
    plasmid_temp_file = run_hmmsearch_and_process(plasmid_output_fasta, plasmid_hmm_folder, 'plasmid')

    if virus_temp_file:
        virus_final_output = os.path.join(output_folder, f'virus_{fasta_basename}_final_output.csv')
        merge_csv(virus_output_csv, virus_temp_file, virus_final_output)

    if plasmid_temp_file:
        plasmid_final_output = os.path.join(output_folder, f'plasmid_{fasta_basename}_final_output.csv')
        merge_csv(plasmid_output_csv, plasmid_temp_file, plasmid_final_output)

    print("Processing complete.")

