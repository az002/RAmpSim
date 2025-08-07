import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

probes_df = pd.read_csv('data/probes.csv')
sequences = probes_df['Sequence']
genes = probes_df['Target Gene Accession']
records = []

for i,(id,s,g) in enumerate(probes_df[['ProbeID','Sequence','Target Gene Accession']].values):
    rec = SeqRecord(Seq(s),id=f'probe_{id}',description=str(g))
    records.append(rec)

with open('data/all_probes.fa','w') as f:
    SeqIO.write(records,f,'fasta')

# e_coli_records = []
# for i,(id,s,g) in enumerate(probes_df[['ProbeID','Sequence','Target Gene Accession']].values):
#     g = str(g)
#     if 'coli' in g or 'Coli' in g or 'Escherichia' in g or 'escherichia' in g:
#         rec = SeqRecord(Seq(s),id=f'probe_{id}',description=g)
#         e_coli_records.append(rec)


# with open('data/e_coli_probes.fa','w') as f:
#     SeqIO.write(e_coli_records,f,'fasta')