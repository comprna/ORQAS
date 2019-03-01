import sys
import re
import logging
from argparse import ArgumentParser, RawTextHelpFormatter

description = \
    "Description:\n\n" + \
    "This script takes a CDS annotation (fasta format) and identifies transcripts that encode for the same cds\n" + \
    "It also provides a list of the genes that have a single possible CDS\n" + \
    "It can be imported to use the functions or run in the command line to write the results in an outputfile.\n"

parser = ArgumentParser(description = description, formatter_class=RawTextHelpFormatter, add_help=False)
parser.add_argument("-i", "--input", required = True,
                    help = "input CDS annotation file (fasta with Ensembl formatted headers)")
parser.add_argument("-o", "--output", required=True,
                    help = "prefix of the output files")


def transcripttoCDS(inFile):
    dictionary = {} #key = Gene_name; value=list of ENST encoding the same transcript.
    sequence = str()
    tid = re.compile('ENS.{0,6}T[0-9]+')
    gid = re.compile('ENS.{0,6}G[0-9]+')
    with open(inFile, 'r') as f:
        for l in f:
            if l.startswith(">"):
                if sequence:
                    if ENSG not in dictionary:
                        dictionary.setdefault(ENSG, {})
                        dictionary[ENSG][sequence] = [ENST]
                    else:
                        if sequence in dictionary[ENSG].keys():
                            dictionary[ENSG][sequence].append(ENST)
                        else:
                            dictionary[ENSG][sequence] = [ENST]

                sequence = str()
                line = l.rstrip('\n')

                #Parses different headers with ENST and ENSG identifiers
                ENST = re.search(tid, line).group()
                ENSG = re.search(gid, line).group()

            else:
                sequence += l.rstrip('\n')

        if sequence:
                    if ENSG not in dictionary:
                        dictionary.setdefault(ENSG, {})
                        dictionary[ENSG][sequence] = [ENST]
                    else:
                        if sequence in dictionary[ENSG].keys():
                            dictionary[ENSG][sequence].append(ENST)
                        else:
                            dictionary[ENSG][sequence] = [ENST]

    logging.info("File %s closed." % inFile)
    #We are not interested in the sequence of the CDS so we can retrieve only the matching isoforms
    #dictionary.update((x[])]) for x in dictionary.items())
    return(dictionary)


def writeOuput(outPref, dictionary):
    f = open(outPref + ".ENSTtoCDS.txt", "w")
    for key, value in [(x,y) for x, y in dictionary.items()]:
        for cds in value:
            f.write(key + "\t" + (":").join(dictionary[key][cds]) + "\n")
    f.close()


def main():
    try:
        args = parser.parse_args()
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - INFO: %(message)s')
        #Recovering the parsed options
        inFile = args.input
        CDSdict = transcripttoCDS(inFile)
        outPref = args.output
        writeOuput(outPref, CDSdict)
        logging.info("Process ended succesfully.")

    except BaseException as err:
         logging.error('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
         logging.error("%s" % err)
         sys.exit(1)


if __name__ == "__main__":
    main()
