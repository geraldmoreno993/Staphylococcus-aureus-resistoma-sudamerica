#Para correr el siguiente codigo, debes tener un combined.tsv proveniente de amrfinderplus
#Simplemente luego de eso escribir en la terminal python consolidate_amr_results.py donde se encuentra tu archivo combined.tsv y este pipeline

import os
import pandas as pd

# Define a function to parse AMRFinderPlus output
def parse_amrfinder_output(file):
    df = pd.read_csv(file, sep='\t', comment='#')
    return df['Gene symbol'].unique()

# Get a list of all AMRFinderPlus output files
amr_files = [f for f in os.listdir() if f.endswith('.amrfinder')]

# Initialize a dictionary to hold presence/absence data
presence_absence_dict = {}

# Populate the dictionary with data from each file
for file in amr_files:
    strain = os.path.splitext(file)[0]
    genes = parse_amrfinder_output(file)
    for gene in genes:
        if gene not in presence_absence_dict:
            presence_absence_dict[gene] = {}
        presence_absence_dict[gene][strain] = 1

# Create a DataFrame from the dictionary
presence_absence_df = pd.DataFrame.from_dict(presence_absence_dict, orient='index').fillna(0)

# Save the DataFrame to a TSV file
presence_absence_df.to_csv('presence_absence.tsv', sep='\t')

print("Consolidated TSV file created: presence_absence.tsv")

