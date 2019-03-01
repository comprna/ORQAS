#!/usr/bin/env python
def build_aa_to_codon(stop):
    aa2codon = {
        'A': ('GCT', 'GCC', 'GCA', 'GCG'),
        'C': ('TGT', 'TGC'),
        'D': ('GAT', 'GAC'),
        'E': ('GAA', 'GAG'),
        'F': ('TTT', 'TTC'),
        'G': ('GGT', 'GGC', 'GGA', 'GGG'),
        'I': ('ATT', 'ATC', 'ATA'),
        'H': ('CAT', 'CAC'),
        'K': ('AAA', 'AAG'),
        'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
        'M': ('ATG',),
        'N': ('AAT', 'AAC'),
        'P': ('CCT', 'CCC', 'CCA', 'CCG'),
        'Q': ('CAA', 'CAG'),
        'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'T': ('ACT', 'ACC', 'ACA', 'ACG'),
        'V': ('GTT', 'GTC', 'GTA', 'GTG'),
        'W': ('TGG',),
        'Y': ('TAT', 'TAC'),
        stop: ('TAA', 'TAG', 'TGA'), #stop codon
    }
    return aa2codon

def convert_codon_to_aa(aa2codon):
    return {codon:aa for aa in aa2codon for codon in aa2codon[aa]}

def build_codon_to_aa(stop):
    aa2codon = build_aa_to_codon(stop)
    return convert_codon_to_aa(aa2codon)

def encode_peptide(tseq, codon2aa):
    return "".join([ codon2aa[ tseq[i:i+3] ] for i in xrange(0,len(tseq),3) ])

def hamming_dist(a,b):
    d = 0
    slen = min(len(a),len(b))
    for i in xrange(slen):
        if a[i]!=b[i]:
            d += 1
    return d

def is_start_codon(c):
    return (hamming_dist(c, 'ATG') < 2)

def check_mid_stop(pseq,stop_codon):
    """ return true if there are stop codon in the middle of the peptide sequence"""
    idx_stop = pseq.find(stop_codon)
    return (idx_stop != -1 and idx_stop != len(pseq)-1)

def check_tseq_dup(tseq,tid,tseq2tids):
    """ return true if transcript sequence is seen before"""
    tseq2tids.setdefault(tseq,[]).append(tid)
    return len(tseq2tids[tseq])>1

def peptide_len(tseq, codon2aa, stop_codon):
    l = len(tseq)/3
    if is_start_codon(tseq[:3]): l -= 1
    if codon2aa[tseq[-3:]] == stop_codon: l -= 1
    return l

def get_gene_list(fa_fname, sep = ' '):
    print "getting gene name list..."
    glist = {}
    tf = open(fa_fname)
    for line in tf:
        if line.startswith('>'):
            gid = line.lstrip('>').rstrip().split(sep)[0]
            glist[gid] = True
    tf.close()
    return glist

def filter_fasta(ifa,ofa, glist, sep=' '):
    from Bio import SeqIO
    print "filtering fasta..."
    ifile = open(ifa, "rU")
    ofile = open(ofa,'w')
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = rec.id.split(sep)[0]
        if not glist[tid]: continue
        SeqIO.write(rec,ofile, "fasta")
    ifile.close()
    ofile.close()

