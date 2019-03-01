#! /usr/bin/python
import sys
import re
import glob
import linecache
import math
from argparse import ArgumentParser, RawTextHelpFormatter


description = "Description:\n\n" + \
              "Obtain uniformity (PME) and periodicity values for each ORF \n" + \
              "given the counts per base files (.base) from the Ribomap output"
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=False)
parser.add_argument("-i", "--input-files", nargs="+", required = True,
                    help = "spaced separated list of files containing the counts per base (.base) from Ribomap output")
parser.add_argument("-o", "--output-file",
                    help="Name for the output file")
parser.add_argument("-c", "--cds",
                    help="File path to the TXTtoCDS file")





def main():
    #Recovering the parsed options
    args = parser.parse_args()
    samples = args.input_files
    out = args.output_file
    merged = args.cds

    gene_t = {}
    genes_c = {}
    for line in open(merged):
        for trans in line.split("\t")[1].rstrip("\n").split(":"):
            gene_t[trans] = [line.split("\t")[0],line.split("\t")[1].rstrip("\n").split(":")]
            try:
                genes_c[line.split("\t")[0]] = genes_c[line.split("\t")[0]] + 1
            except:
                genes_c[line.split("\t")[0]] = 1

    s_counts_base = {}
    s_counts_rna = {}
    for sample in samples:
        s_counts_base[sample] = {}
        s_counts_rna[sample] = {}
        for line in open(sample): #####Base output RIBOMAP
            if "tid:" in line:
                trans = line.split(" ")[1].rstrip("\n").split(".")[0]
            elif ("ribo profile:" in line) and (line.startswith("ribo")):
                s_counts_base[sample][trans] = line.split(": ")[1].rstrip("\n").split(" ")
                del s_counts_base[sample][trans][-1]
                s_counts_base[sample][trans] = list(map(float, s_counts_base[sample][trans]))
            elif "mRNA profile:" in line:
                s_counts_rna[sample][trans] = line.split(": ")[1].rstrip("\n").split(" ")
                del s_counts_rna[sample][trans][-1]
                s_counts_rna[sample][trans] = list(map(float, s_counts_rna[sample][trans]))

    counts_base = {}
    counts_rna = {}
    for sample in samples:
        for trans in s_counts_base[sample]:
            if trans not in gene_t:
                continue

            if len(s_counts_base[sample][trans]) <= 1:
                continue

            if not (trans in counts_base):
                counts_base[trans] = []
                for i,pos in enumerate(s_counts_base[sample][trans]):
                    counts_base[trans].append(0)
                    counts_base[trans][i] = counts_base[trans][i] + pos

            else:
                for i,pos in enumerate(s_counts_base[sample][trans]):
                    counts_base[trans][i] = counts_base[trans][i] + pos

        for trans in s_counts_rna[sample]:
            if trans not in gene_t:
                continue

            if len(s_counts_rna[sample][trans]) <= 1:
                continue

            if not (trans in counts_rna):
                counts_rna[trans] = []
                for i,pos in enumerate(s_counts_rna[sample][trans]):
                    counts_rna[trans].append(0)
                    counts_rna[trans][i] = counts_rna[trans][i] + pos

            else:
                for i,pos in enumerate(s_counts_rna[sample][trans]):
                    counts_rna[trans][i] = counts_rna[trans][i] + pos


    #Merge CDS
    new_counts_rna = {}
    for trans in counts_rna:
        if not trans in gene_t:
            continue

        if len(gene_t[trans][1]) > 1:
            lists = []
            for t in gene_t[trans][1]:
                if t in counts_rna:
                    lists.append(counts_rna[t])
            if not ":".join(gene_t[trans][1]) in new_counts_rna:
                new_counts_rna[":".join(gene_t[trans][1])] = [sum(sublist) for sublist in zip(*lists)]
        else:
            new_counts_rna[trans] = counts_rna[trans]

    new_counts_base = {}
    for trans in counts_base:
        if not trans in gene_t:
            continue

        if len(gene_t[trans][1]) > 1:
            lists = []
            for t in gene_t[trans][1]:
                if t in counts_base:
                    lists.append(counts_base[t])
            if not ":".join(gene_t[trans][1]) in new_counts_base:
                new_counts_base[":".join(gene_t[trans][1])] = [sum(sublist) for sublist in zip(*lists)]
        else:
            new_counts_base[trans] = counts_base[trans]

    #Uniformity
    pme = {}
    for trans in new_counts_base:
        new_counts_ribo = []
        for n2,pos in enumerate(new_counts_base[trans]):
            if n2 % 3 == 0:
                new_counts_ribo.append(pos)
            else:
                new_counts_ribo[-1] = new_counts_ribo[-1] + pos

        l = len(new_counts_ribo)
        n = sum(new_counts_ribo)
        if n < 10:
            pme_t = 0
            hm = 0
        else:
            if n > l:
                rl = float(1)
            else:
                rl = float(round(l/n))
            c = 0
            nr = 0
            h = 0
            hm = 0
            for i,count in enumerate(new_counts_ribo):
                nr = nr + count
                c += 1
                if c == rl:
                    p = float(nr)/float(n)
                    if p > 0:
                        h = h + (p*math.log(p,2))
                    pm = (n/(l/rl))/n
                    hm = hm + (pm*math.log(pm,2))
                    c = 0
                    nr = 0

            pme_t = h/hm
        pme[trans] = pme_t

    lex = [0,0,0]
    counts = {}
    out1 = open(out, "w+")

    out1.write("cds\tgene\tn_cds\tcov_ribo\tcov_rna\tf1\tf2\tribo_reads\tpme\n")
    for trans in new_counts_base:
        try:
            cov_ribo = float(len(new_counts_base[trans]) - new_counts_base[trans].count(0)) / len(new_counts_base[trans]) * 100
        except:
            cov_ribo = 0
        try:
            cov_rna = float(len(new_counts_rna[trans]) - new_counts_rna[trans].count(0))  / len(new_counts_rna[trans]) * 100
        except:
            cov_rna = 0

        if not trans.split(":")[0] in counts:
            counts[trans.split(":")[0]] = 0

        f = [0, 0, 0, 0]
        for i,count in enumerate(new_counts_base[trans]):
            if float(count) > 0:
                if i % 3 == 0:
                    f[1] = f[1] + count
                elif i % 3 == 1:
                    f[2] = f[2] + count
                elif i % 3 == 2:
                    f[3] = f[3] + count

        if (f[1] + f[2] + f[3]) >= 10:
            f1 = f[1]/(f[1] + f[2] + f[3])
            f2 = f[2]/(f[1] + f[2] + f[3])
        else:
            f1 = 0
            f2 = 0

        out1.write(trans + "\t" + gene_t[trans.split(":")[0]][0] + "\t" + str(genes_c[gene_t[trans.split(":")[0]][0]]) + "\t" + str(cov_ribo) + "\t" + str(cov_rna) + "\t" + str(f1) + "\t" + str(f2) + "\t" + str(f[1] + f[2] + f[3]) + "\t" + str(pme[trans]) + "\n")

    out1.close()


if __name__ == "__main__":
    main()
