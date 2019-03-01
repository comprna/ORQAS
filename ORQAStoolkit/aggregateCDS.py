import sys
import re
import logging
from argparse import ArgumentParser, RawTextHelpFormatter

description = "Description:\n\n" + \
              "Obtain aggregated CDS TPM or count (or both) for transcripts encoding the same CDS taking as a reference \n" + \
              "a CDS annotation or a txt of aggregate transcripts in the following format: ENST1;ENST2"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=False)
parser.add_argument("-i", "--input-file", required=True,
                    help="input file of TPM, OPM or count values where each column is a sample")
parser.add_argument("-s", "--separator", default="\t",
                    help="field separator")
parser.add_argument("-f", "--fields", nargs="+", type=int, default=None,
                    help="fields from the input file that will be considered for analysis (takes all the fields by default)")
parser.add_argument("-o", "--output-file",
                    help="Name for the output file")
parser.add_argument("-c", "--cds",
                    help="File path to the TXTtoCDS file")
parser.add_argument("-a","--annot",
                    help="File path to the Ensembl cds annotation")


def readQuantFile(quantFile, seen_header, separator, fields):
    dictionary = {}  # key = keyField; value=list of fields to join from different files.
    with open(quantFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(separator)
            if not seen_header:
                dictionary["header"] = line
                seen_header = True
            else:
                #dictionary.setdefault(line[(keyField-1)], {})
                #Load all the fields from the file
                dictionary[line[0]] = []
                if fields == None:
                    dictionary[line[0]] += [float(i) for i in line[1:]]
                else:
                    for f in fields:
                        dictionary[line[0]].append(float(line[f-1]))
    logging.info("File %s closed." % quantFile)
    return(dictionary)


def sumvalues(quantDict, ENSTlist, samplenum):
    LoL = []
    for id in ENSTlist:
        try:
            LoL.append(quantDict[id])
        except:
            LoL.append([0]*samplenum)
    CDSvalue = [sum(sample) for sample in zip(*LoL)] #unpack and zip
    return(CDSvalue)


def calcCDSexpr(quantDict, CDSinfo):
    with open(CDSinfo, 'r') as f:
        dictionary = {}
        samplenum = len(list(quantDict["header"])) - 1
        for l in f:
            line = l.rstrip('\n').split("\t")
            ENST = line[1]
            ENSTlist = ENST.split(":")
            if len(ENSTlist) > 1:
                dictionary[ENST] = sumvalues(quantDict, ENSTlist, samplenum)
            else:
                if ENSTlist[0] in quantDict.keys():
                    dictionary[ENST] = quantDict[ENSTlist[0]]
                else:
                    dictionary[ENST] = [0]*(samplenum+1)
    logging.info("File %s closed." %CDSinfo)
    return(dictionary)

def printOutFiles(values, head, outFile):
    f = open(outFile, 'w')
    #Writing the header in the output file
    headStr = ("\t").join(head)
    f.write("%s\n" %(headStr))

    for ENST,tpm in values.items():
        f.write("%s" %(ENST))
        for value in tpm:
            f.write("\t%.6f" %(float(value)))
        f.write("\n")
    f.close()
    logging.info("Output file %s closed." % outFile)


def main():
    try:
        args = parser.parse_args()
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - INFO: %(message)s')
        #Recovering the parsed options
        quantFile = args.input_file
        fields = args.fields
        separator = args.separator
        outputFile = args.output_file
        if args.cds:
            ENSTFile = args.cds
        elif args.annot:
            AnnotFile = args.annot
        else:
            logging.error("You should provide an annotation file or the aggregated transcript ids file.")

        seen_header = False

        logging.info('Reading file %s' %quantFile)
        quantDict = readQuantFile(quantFile, seen_header, separator, fields)
        logging.info("Retrieving aggregated TPM values")
        newTPMs= calcCDSexpr(quantDict, ENSTFile)

        ##WRITING OUTPUT
        header = quantDict["header"]
        logging.info("Writing output to %s" %outputFile)
        printOutFiles(newTPMs, header, outputFile)
        logging.info("Process ended succesfully.")


    except BaseException as err:
         logging.error('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
         logging.error("%s" % err)
         sys.exit(1)

if __name__ == "__main__":
    main()


