import pandas as pd
import sys

# Function to extract scaffold names from hmm_target_accession
def extract_scaffold_name(hmm_target_accession):
    return '_'.join(hmm_target_accession.split('_')[:2])

# Function to process the HMMsearch results
def process_hmmsearch(input_file):
    # Read the CSV file
    df = pd.read_csv(input_file, delimiter=',')

    # Adjust column names based on the actual data
    if 'hmm_target_accession' in df.columns:
        df['scaffold_name'] = df['hmm_target_accession'].apply(extract_scaffold_name)
    else:
        print("Error: 'hmm_target_accession' column not found in the DataFrame.")
        sys.exit(1)

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
    if len(sys.argv) != 4:
        print("Usage: python3 PhagePlasmidQuantifier.py input1.csv input2.csv output.csv")
        sys.exit(1)

    hmmsearch_input = sys.argv[1]
    genomad_input = sys.argv[2]
    output_file = sys.argv[3]

    # Process HMMsearch results
    hmmsearch_results = process_hmmsearch(hmmsearch_input)

    # Save the intermediate results to a temporary file
    temp_file = 'PhageToPlasmidHits.csv'
    hmmsearch_results.to_csv(temp_file, index=False)

    # Merge the results with genomad output
    merge_csv(genomad_input, temp_file, output_file)
