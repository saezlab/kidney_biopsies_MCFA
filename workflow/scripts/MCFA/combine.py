# 31.10.2023 Charlotte Boys

# Load libraries
import pandas as pd

# Get the list of input CSV files
input_files = snakemake.input

# Initialize an empty DataFrame to store the combined data
combined_data = pd.DataFrame()

# Iterate over the input files
for file in input_files:
    # Extract the 'nfactors' value from the filename
    nfactors = int(file.split('~')[-1].split('.')[0])
    # Read the CSV file and add the 'n_factors' column
    data = pd.read_csv(file)
    data['n_factors'] = nfactors

    # Append the data to the combined DataFrame
    combined_data = combined_data.append(data, ignore_index=True)

# Write the combined data to the output file
combined_data.to_csv(snakemake.output[0], index=False)
