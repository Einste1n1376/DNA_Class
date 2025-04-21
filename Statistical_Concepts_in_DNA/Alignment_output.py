# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt

# ---- Load sequences ----
query_seqs = list(SeqIO.parse("Statistical_Concepts_in_DNA/fasta_data/Database.fasta", "fasta"))
db_seqs = list(SeqIO.parse("Statistical_Concepts_in_DNA/fasta_data/Database.fasta", "fasta"))

# ---- Store all alignments for combined plotting ----
alignment_data = []

for query in query_seqs:
    for db in db_seqs:
        alignments = pairwise2.align.localms(query.seq, db.seq, 2, -1, -0.5, -0.1)
        if alignments:
            top_alignment = alignments[0]
            q_aln, db_aln, score, start, end = top_alignment

            print(f"\n--- Alignment: {query.id} vs {db.id} ---")
            print(f"Score: {score}")
            print(format_alignment(*top_alignment))

            # Match line: '|' for match, ' ' otherwise
            match_line = "".join('|' if q == d and q != '-' else ' ' for q, d in zip(q_aln, db_aln))

            # Store for later plotting
            alignment_data.append({
                'query_id': query.id,
                'db_id': db.id,
                'q_aln': q_aln,
                'match_line': match_line,
                'db_aln': db_aln,
                'score': score
            })

# ---- Create combined plot with subplots ----
n_alignments = len(alignment_data)
fig_height = max(2, n_alignments * 2)  # Auto adjust height
fig, axes = plt.subplots(n_alignments, 1, figsize=(15, fig_height))

# If only one alignment, wrap axes in a list
if n_alignments == 1:
    axes = [axes]

for ax, aln in zip(axes, alignment_data):
    q = aln['q_aln']
    m = aln['match_line']
    d = aln['db_aln']

    for i, (q_char, m_char, d_char) in enumerate(zip(q, m, d)):
        ax.text(i, 2, q_char, ha='center', va='center', fontsize=10, fontfamily='monospace')
        ax.text(i, 1, m_char, ha='center', va='center', fontsize=9, color='gray', fontfamily='monospace')
        ax.text(i, 0, d_char, ha='center', va='center', fontsize=10, fontfamily='monospace')

    ax.set_xlim(-1, len(q) + 1)
    ax.set_ylim(-0.5, 2.5)
    ax.set_title(f"{aln['query_id']} vs {aln['db_id']} | Score: {int(aln['score'])}", fontsize=12)
    ax.axis('off')

plt.tight_layout()
plt.show()