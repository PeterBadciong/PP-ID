import pandas as pd
import sys
import subprocess
import re

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
        # Write header
        header = "\t".join(column_names_container.keys()) + "\n"
        file.write(header)
        # Write data
        for row in parsed_data:
            row_data = "\t".join(row) + "\n"
            file.write(row_data)

# Function to extract scaffold names from hmm_target_accession
def extract_scaffold_name(hmm_target_accession):
    return '_'.join(hmm_target_accession.split('_')[:2])

# Function to process the HMMsearch results
def process_hmmsearch(input_file):
    parsed_data = parse_hmmtbl(input_file)
    df = pd.DataFrame(parsed_data, columns=column_names)

    # Add scaffold_name based on target_accession
    df['scaffold_name'] = df['hmm_target_accession'].apply(extract_scaffold_name)

    # Filter necessary columns
    filtered_df = df[['scaffold_name', 'query_protein_accession']]

    # Remove duplicate rows based on scaffold name and query protein
    unique_scaffolds = filtered_df.drop_duplicates()

    # Count unique protein hits for each scaffold
    unique_counts = unique_scaffolds.groupby('scaffold_name')['query_protein_accession'].nunique()

    # Convert the result to a DataFrame
    result_df = unique_counts.reset_index()
    result_df.columns = ['scaffold_name', 'unique_protein_hits']

    return result_df

# Function to merge data and calculate ratios
def merge_csv(input_csv1, input_csv2, output_csv):
    # Read CSV files into pandas DataFrames
    df1 = pd.read_csv(input_csv1)
    df2 = pd.read_csv(input_csv2)

    # Merge DataFrames on 'scaffold_name' column
    merged_df = pd.merge(df1, df2, on='scaffold_name')

    # Calculate proteins over genes ratio
    merged_df['proteins_over_genes'] = merged_df['unique_protein_hits'] / merged_df['n_genes']

    # Select desired columns for output
    output_df = merged_df[['scaffold_name', 'n_genes', 'unique_protein_hits', 'proteins_over_genes']]

    # Write to output CSV file
    output_df.to_csv(output_csv, index=False)
    print(f"Output saved to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 PhagePlasmidQuantifier.py hmm_db.hmm input_sequences.fasta genomad_output.csv final_output.csv")
        sys.exit(1)

    hmm_db = sys.argv[1]
    input_sequences = sys.argv[2]
    genomad_input = sys.argv[3]
    output_file = sys.argv[4]

    # Run HMMsearch
    hmmsearch_output = 'hmmsearch_output.tblout'
    subprocess.run(['hmmsearch', '--tblout', hmmsearch_output, hmm_db, input_sequences], check=True)

    # Process HMMsearch results
    hmmsearch_results = process_hmmsearch(hmmsearch_output)

    # Save the intermediate results to a temporary file
    temp_file = 'PhageToPlasmidHits.csv'
    hmmsearch_results.to_csv(temp_file, index=False)

    # Merge the results with genomad output
    merge_csv(genomad_input, temp_file, output_file)

