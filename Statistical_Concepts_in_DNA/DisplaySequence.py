# -*- coding: utf-8 -*-

import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import Bio

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

file_path ='C:/Biopython/genome.fasta'

try:
 sequence_record = SeqIO.read(file_path, 'fasta')
 print(f'Sequence ID:{sequence_record.id}')
 print(f'DNA sEQUENCE:{sequence_record.seq}')
except Exception as e:
 print("Error:", e)
 
from Bio.SeqUtils import gc_fraction
gc_content = gc_fraction(sequence_record.seq)
print("GC_content:", gc_content)

from collections import Counter
sequence_length = len(sequence_record)
nucleotide_counts = Counter(sequence_record.seq)
print('Sequence length:', sequence_record)
print('Nucleotide counts', nucleotide_counts)

import matplotlib.pyplot as plt
plt.bar(nucleotide_counts.keys(), nucleotide_counts.values())
plt.xlabel('Nucleotide')
plt.ylabel('Count')
plt.title('Nucleotide Distrubition')
plt.show()