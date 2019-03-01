#!/usr/bin/env python
import sys
from Bio import SeqIO

if len(sys.argv)!=5:
    print "Usage: python build_contaminant.py rrna_fa trna_fa trna_species ofile"
    exit(1)

rrna_fa = sys.argv[1]
trna_fa = sys.argv[2]
trna_specie = sys.argv[3]
out_fa = sys.argv[4]

ofile = open(out_fa, 'w')
i = 0
trna_file = open(trna_fa)
for rec in SeqIO.parse(trna_file, "fasta"):
    if rec.id.startswith(trna_specie):
        SeqIO.write(rec, ofile, "fasta")
        i += 1
trna_file.close()
print "{0} tRNA seqs included".format(i)

j = 0
rrna_file = open(rrna_fa)
for rec in SeqIO.parse(rrna_file, "fasta"):
    theader = rec.description
    words = theader.split()
    for w in words:
        if w.startswith("gene_biotype:"):
            gbt = w.lstrip("gene_biotype:")
            if gbt == "rRNA":
                SeqIO.write(rec, ofile, "fasta")
                j += 1
rrna_file.close()
print "{0} rRNA seqs included".format(j)
ofile.close()

print "include {0} contaminated sequences in {1}".format(i+j, out_fa)


# tf = open("Homo_sapiens.GRCh38.ncrna.fa")
# gb = set([])
# tb = set([])
# for line in tf:
#     if line.startswith(">"):
#         words = line.split()
#         for w in words:
#             if w.startswith("gene_biotype:"):
#                 gb.add(w.lstrip("gene_biotype:"))
#             elif w.startswith("transcript_biotype:"):
#                 tb.add(w.lstrip("transcript_biotype:"))
# print gb
# print tb
