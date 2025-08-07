import pandas as pd
import numpy as np
from scipy.optimize import lsq_linear
from Bio.Seq import Seq

# Read the CSV file
data = pd.read_csv('data/TNN_data.csv')
data = data.rename(columns={'sequence': 'sequence_A', 'sequence.1': 'sequence_B'})

# Replace U with T in both sequence columns
data['sequence_A'] = data['sequence_A'].str.replace('U', 'T')
data['sequence_B'] = data['sequence_B'].str.replace('U', 'T')

# Generate all possibletrimers
bases = ['A', 'C', 'G', 'T']
trimers = [Seq(b1 + b2 + b3) for b1 in bases for b2 in bases for b3 in bases]

# Create trimer index mapping
used_trimers = set()
trimer_to_index = {}

index = 0
for trimer in trimers:
    if trimer not in used_trimers:
        rc = trimer.reverse_complement()
        trimer_to_index[trimer] = index
        trimer_to_index[rc] = index
        index += 1
        used_trimers.add(trimer)
        used_trimers.add(rc)

def get_trimer_profile(sequence):
    # Initialize profile with zeros
    profile = np.zeros(len(set(trimer_to_index.values())))
    
    # Count trimers in sequence
    for i in range(len(sequence) - 2):
        trimer = Seq(sequence[i:i+3])
        if trimer in trimer_to_index:
            profile[trimer_to_index[trimer]] += 1
    
    return profile

# Calculate profiles for all sequences in sequence_A
sequence_profiles = data['sequence_A'].apply(get_trimer_profile)
sequence_profiles = np.array(sequence_profiles.tolist())

# Print sequence profiles with entries separated by commas
TNN_profile = lsq_linear(sequence_profiles, data['dG'],bounds=(-np.inf,0))['x']

# for i in range(0, len(TNN_profile), 4):
#     row_items = TNN_profile[i:min(i+4, len(TNN_profile))]
#     print(" ".join([f"{item}, " for item in row_items]))

# Print trimer to index mapping, 4 per row
items = sorted([(str(k), v) for k, v in trimer_to_index.items()], key=lambda x: x[1])
base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
for i in range(0, len(items), 4):
    row_items = items[i:min(i+4, len(items))]
    for i in range(len(row_items)):
        val = 0
        for c in row_items[i][0]:
            val = val << 2 | base_map[c]
        row_items[i] = (val, row_items[i][1])
    print(" ".join([f"({item[0]}, {item[1]})," for item in row_items]))



# np.savetxt('data/TNN_profile.csv', TNN_profile, delimiter=',')