#!/usr/bin/env python
import sys
import re
from Bio import SeqIO
from translation import *

def get_transcript_frame(gtfname):
    """
    gtf fields:
    0: seqname
    1: source
    2: feature
    3: start
    4: end
    5: score
    6: strand
    7: frame
    8: attribute
      gene_id ENSGXXXXXXXXXXX *
      transcript_id ENSTXXXXXXXXXXX *
      gene_type list of biotypes
      gene_status {KNOWN, NOVEL, PUTATIVE}
      gene_name string
      transcript_type list of biotypes
      transcript_status {KNOWN, NOVEL, PUTATIVE}
      transcript_name string
      level 1 (verified loci), 2 (manually annotated loci), 3 (automatically annotated loci)
    """
    print "getting frame"
    gtfile = open(gtfname)
    tid2frame= {}
    for line in gtfile:
        if line.startswith("#"): continue
        words = line.rstrip('\n').split('\t')
        if words[2] == "CDS":
            frame = int(words[7])
            assert words[8].find("exon_number") != -1
            attributes = words[8].split("; ")
            for a in attributes:
                aname, aval = a.lstrip().split(" ")
                aval = aval.lstrip('"').rstrip('"')
                if aname == "transcript_id":
                    tid = aval
                if aname == "exon_number":
                    enum = int(aval)
            if tid in tid2frame:
                e_pre, f_pre = tid2frame[tid]
                if enum < e_pre:
                    tid2frame[tid] = (enum, frame)
            else:
                tid2frame[tid] = (enum, frame)
    gtfile.close()
    return {tid:frame for tid,(enum,frame) in tid2frame.iteritems()}

def get_gencode_cds_range(transcript_fa, tid2frame, sep=' '):
    """
    transcript.fa header:
    0 transcript-id|
    1 gene-id|
    2 Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
    3 Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
    4 transcript-name|
    5 gene-name|
    6 sequence-length|
    7 5'-UTR (3'-UTR if reverse strand) location in the transcript|
    8 CDS location in the transcript|
    9 3'-UTR (5'-UTR if reverse strand) location in the transcript
    """
    print "getting cds range from gencode header..."
    tid2cds = {}
    tid2theader = {}
    tf = open(transcript_fa)
    for line in tf:
        if line.startswith('>'):
            theader = line.lstrip('>').rstrip()
            tid = theader.split(sep)[0]
            twords = theader.split('|')
            if len(twords)<8: continue
            for w in twords[7:]:
                if w.startswith("CDS:"):
                    start, stop = map(int, w.lstrip("CDS:").split('-'))
                    frame = tid2frame[tid]
                    start += frame
                    r = (stop-start+1)%3
                    stop -= r
                    tid2cds[tid] = (start-1, stop)
                    tid2theader[tid] = theader
    tf.close()
    return tid2cds, tid2theader

def compare_gencode_pseq(pconvert, pseq):
    if pseq.startswith('X'): pseq = pseq[1:]
    pconvert.rstrip(stop_codon)
    l = min(len(pconvert), len(pseq))
    cset = set([])
    for i in xrange(1,l):
        if pconvert[i] != pseq[i]:
            cset.add((pconvert[i], pseq[i]))
    return cset

def find_tid_in_pheader(s, prefix=("ENST", "ENSMUST")):
    words = s.split('|')
    for w in words:
        if w.startswith(prefix):
            return w
    return ''

def gencode_codon_check(tfname, pfname, glist, tid2cds, sep=' '):
    """
    exclude :
    1) transcripts that I cannot correctly encode to peptides
    2) transcripts with stop codon in the middle
    3) transcripts with duplicated sequences
    4) peptide sequence length less than 3 after getting rid of start and end of theseq
    """
    Mcnt = 0
    Ucnt = 0
    tcnt = 0
    ms_tid = []
    short_tid = []
    dup_tid = []
    failed_tid = []
    ccset = set([])
    tseq2tids = {}
    print "checking transcripts..."
    tfile = open(tfname, "rU")
    pfile = open(pfname, "rU")
    for trec, prec in zip(SeqIO.parse(tfile, "fasta"), SeqIO.parse(pfile, "fasta")):
        twords = trec.id.split('|')
        pwords = prec.id.split('|')
        tid = twords[0]
        tidp = find_tid_in_pheader(prec.id)
        # make sure transcript id between two files match
        assert tid == tidp
        # make sure sequence length in header match the real length
        assert len(trec) == int(twords[6]) 
        assert len(prec) == int(pwords[-1])
        tcnt += 1
        if len(twords)<8:
            glist[tid] = False
            continue
        start, stop = tid2cds[tid]
        tseq = str(trec.seq)
        pconvert = encode_peptide(tseq[start: stop], codon2aa)
        if is_start_codon(tseq[start:start+3]): Mcnt += 1
        if pconvert.endswith(stop_codon): Ucnt += 1
        if check_mid_stop(pconvert, stop_codon):
            glist[tid] = False
            ms_tid.append(tid)
        elif check_tseq_dup(tseq, tid, tseq2tids):
            glist[tid] = False
            dup_tid.append(tid)
        elif peptide_len(tseq[start: stop], codon2aa, stop_codon)<3:
            glist[tid] = False
            short_tid.append(tid)
        else: 
            pseq = str(prec.seq)
            cset = compare_gencode_pseq(pconvert, pseq)
            if cset:
                glist[tid] = False
                failed_tid.append(tid)
                ccset |= cset
                print twords
                print pwords
                print "frame: {0}, pconvert_len: {1}, pseq_len: {2}, CDS:{3}-{4}".format(tid2frame[tid], len(pconvert), len(pseq), start, stop)
                print pconvert
                print pseq
                print cset
    tfile.close()
    pfile.close()
    print "# peptides with start codon: {0} ({1})".format(Mcnt, Mcnt/float(tcnt))
    print "# peptides ends with stop codon: {0} ({1})".format(Ucnt, Ucnt/float(tcnt))
    print "# middle stop codon: {0} ({1})".format(len(ms_tid), len(ms_tid)/float(tcnt))
    print "# duplicate transcripts: {0} ({1})".format(len(dup_tid), len(dup_tid)/float(tcnt))
    print "# tiny peptides (len<3) : {0} ({1})".format(len(short_tid), len(short_tid)/float(tcnt))
    print "{0} transcripts failed to encode correctly".format(len(failed_tid))
    print "codons that are encoded wrong: ", ccset
    print "total transcripts: {0}".format(tcnt)
    return glist

def filter_gencode_pfasta(ifa, ofa, glist):
    """
    peptide.fa header
    0 peptide-id -- newer version is like this
    1 transcript-id|
    2 gene-id|
    3 Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
    4 Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
    5 transcript-name|
    6 gene-name|
    7 sequence-length
    """
    print "filtering peptide fasta..."
    ifile = open(ifa, "rU")
    ofile = open(ofa,'w')
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = find_tid_in_pheader(rec.id)
        if not glist[tid]: continue
        SeqIO.write(rec,ofile, "fasta")
    ifile.close()
    ofile.close()

def record_gencode_cds_range(ofn, tid2cds, tid2theader, glist):
    print "writing gencode cds range file..."
    ofile = open(ofn, 'w')
    for tid in tid2cds:
        if glist[tid]:
            start, stop = tid2cds[tid]
            ofile.write("{0}\t{1}\t{2}\n".format(tid2theader[tid], start, stop))
    ofile.close()

def build_new_name(old_name, suffix):
    return old_name[:old_name.rfind('.')]+suffix

def main():
    if len(sys.argv)!=4:
        print "Usage: python filter_gencode_transcript.py gtf_fname transcript_fa peptide_fa"
        exit(1)
    global stop_codon, codon2aa, tid2frame
    stop_codon = 'U'
    codon2aa = build_codon_to_aa(stop_codon)
    gtf_fn = sys.argv[1]
    tfa = sys.argv[2]
    pfa = sys.argv[3]
    tid2frame = get_transcript_frame(gtf_fn)
    glist = get_gene_list(tfa, '|')
    tid2cds, tid2theader = get_gencode_cds_range(tfa, tid2frame, '|')
    glist = gencode_codon_check(tfa, pfa, glist, tid2cds, '|')
    print "total included transcripts: {0}".format(sum(glist.values()))
    # output filterred results
    cds_fn = build_new_name(tfa, "_cds.txt")
    otfa = build_new_name(tfa, "_filter.fa")
    opfa = build_new_name(pfa, "_filter.fa")
    record_gencode_cds_range(cds_fn, tid2cds, tid2theader, glist)
    filter_fasta(tfa, otfa, glist, '|')
    filter_gencode_pfasta(pfa, opfa, glist)
    gencode_codon_check(otfa, opfa, glist, tid2cds, '|')

if __name__ == "__main__": main()
