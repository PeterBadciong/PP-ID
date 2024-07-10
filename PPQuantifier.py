import pandas as pd
import sys

# Function to extract scaffold names from hmm_target_accession
def extract_scaffold_name(hmm_target_accession):
    return '_'.join(hmm_target_accession.split('_')[:2])

def main(input_file, output_file):
    # Read the CSV file
    df = pd.read_csv(input_file, delimiter=',')

    # Adjust column names based on the actual data
    # Assuming hmm_target_accession should be used for scaffold_name extraction
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

    # Save the result to a new CSV file
    result_df.to_csv(output_file, index=False)

    print("Unique protein hits have been saved to", output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 PhagePlasmidQuantifier.py inputFile outputFile")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    main(input_file, output_file)
