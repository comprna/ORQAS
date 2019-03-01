import sys
import glob
import os
import logging
from math import log, exp
from argparse import ArgumentParser, RawTextHelpFormatter

description = \
    "Description:\n\n" + \
    "This script takes ribosome counts form multiple Ribomap output .stats file\n" + \
    "and returns them transformed into TPMs in a tab formated file with as columns as input files\n" + \
    "mereged according to a common identifier.\n" + \
    "If a identifier do not appear in a file, the corresponding field will\n" + \
    "be set up to 0."

parser = ArgumentParser(description = description, formatter_class=RawTextHelpFormatter, add_help=False)
parser.add_argument("-i", "--input-dirs", nargs="+", required = True,
                    help = "spaced separated list of directories containing both the output and sm_quant directories with the respective Ribomap and Salmon outputs (as it is generated when running ORQAS Pipeline)")
parser.add_argument("-s", "--separator", default = "\t",
                    help = "field separator")
parser.add_argument("-o", "--output-file", required=True,
                    help = "prefix the output file.")
parser.add_argument("-k", "--key-field", default="1",
                    help = "id field in Salmon and Effective length files")
parser.add_argument("-p", "--pattern", default="",
                    help = "pattern to identify .stats file in case there are multiple in the same directory (the pattern should be other than the .stats extension)")
parser.add_argument("--no-header", action="store_true",
                    help = "use it if the file has no header.")
parser.add_argument("-l", "--ef-len", nargs="+", default=False,
                    help="space separated list of files with the CDS effective length "
                        "in case it is different from the Salmon file - MUST BE IN THE SAME ORDER AS THE INPUT FILES")
parser.add_argument("-lf", "--ef-len-field", default=2, help="field of the effective length in the --ef-len files")


def readQuantFile(quantFile, seen_header, separator):
    dictionary = {} #key = keyField; value=list of fields to join from different files.
    with open(quantFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(separator)
            if not seen_header:
                seen_header = True
                continue
            else:
                #dictionary.setdefault(line[(keyField-1)], {})
                #Load all the fields from the file
                key = line[0].split(".")
                dictionary[key[0]] = [line[2]]
    with open(quantFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(separator)
            if not seen_header:
                seen_header = True
                continue
            else:
                #dictionary.setdefault(line[(keyField-1)], {})
                #Load all the fields from the file
                key = line[0].split(".")
                if key[0] in dictionary:
                #Append Effective Length, TPM value and Read Counts from salmon
                    dictionary[key[0]].extend(line[slice(3,5)])

        logging.info("File %s closed." % quantFile)
        return(dictionary)


def readEffLenFile(quantFile, seen_header, separator):
    dictionary = {} #key = keyField; value=list of fields to join from different files.
    with open(quantFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(separator)
            if not seen_header:
                seen_header = True
                continue
            else:
                #dictionary.setdefault(line[(keyField-1)], {})
                #Load all the fields from the file
                key = line[0].split(".")
                dictionary[key[0]] = [line[2]]
                #Append Effective Length, TPM value and Read Counts from salmon
    with open(quantFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(separator)
            if not seen_header:
                seen_header = True
                continue
            else:
                #dictionary.setdefault(line[(keyField-1)], {})
                #Load all the fields from the file
                key = line[0].split(".")
                if key[0] in dictionary:
                #Append Effective Length, TPM value and Read Counts from salmon
                    dictionary[key[0]].extend(line[slice(3,5)])

        logging.info("File %s closed." % quantFile)
        return(dictionary)


def readRiboFile(riboFile):
    riboDict = {} #key = keyField; value=list of fields to join from different files.
    with open(riboFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(" ")
            if line[0] == "tid:":
                #riboDict.setdefault(line[1], {})
                keyname = line[1].split(".")
                riboDict[keyname[0]] = []
            elif line[0] == "rabd:":
                #Append Riboseq counts
                riboDict[keyname[0]].append(line[1])

        logging.info("File %s closed." % riboFile)
        return(riboDict)


def calcTPM(riboDict, quantDict):
    rateList = []
    for key in riboDict.keys():
        counts = float(riboDict[key][0])
        if counts == 0:
            riboDict[key].append(0)
        else:
            efflen = float(quantDict[key][0])
            rate = log(counts) - log(efflen)
            rateList.append(float(rate))
            riboDict[key].append(rate)
    erateList = map(lambda x: exp(x), rateList)
    denom = log(sum(erateList))
    riboDict.update((x, [y[0], exp(y[1] - denom + log(1e6))]) for x, y in riboDict.items() if y[1] != 0)
    return(riboDict)


def getOutDict(outDict, riboDict, quantDict, dictID):
    logging.info("Calculating TPMs for Ribomap counts...")
    TPMriboDict = calcTPM(riboDict, quantDict)

    for key in quantDict.keys():
        if key not in outDict:
            outDict.setdefault(key, {})
            outDict[key][dictID] = []
        else:
            outDict[key][dictID] = []

        ##Introduce RC values in the output dict
        outDict[key][dictID].append(str(quantDict[key][2]))
        if key in TPMriboDict:
            outDict[key][dictID].append(str(riboDict[key][0]))
        else:
            outDict[key][dictID].append(str(0))

        ##Introduce TPM values in the output dict
        outDict[key][dictID].append(str(quantDict[key][1]))
        if key in TPMriboDict:
            outDict[key][dictID].append(str("{:.6f}".format(TPMriboDict[key][1])))
        else:
            outDict[key][dictID].append(str(0))

    return(outDict)


def printOutFiles(outDict,outputFile, ids):
        #Generating new header
        outputHeader = []
        for ID in ids:
            outputHeader.extend(((ID + "_rna"), (ID + "_ribo")))

        f = open(outputFile + ".readCounts.txt", 'w')
        g = open(outputFile + ".Abundance.txt", 'w')
        #Writing the header in the output file
        f.write("\t".join(outputHeader) + "\n")
        g.write("\t".join(outputHeader) + "\n")
        #Looping through all the unique IDs avoiding the "header" .
        for key, value in [(x,y) for x, y in outDict.items()]:
            abline = []
            abline.append(key)    #Adding the common id
            rcline = []
            rcline.append(key)    #Adding the common id

            for ID in ids:
                rcline += value.get(ID)[slice(0,2)]
                abline += value.get(ID)[slice(2,4)]

            f.write("\t".join(rcline) + "\n")
            g.write("\t".join(abline) + "\n")
        f.close()
        g.close()
        logging.info("Output files in %s closed." % outputFile)


def main():
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - INFO: %(message)s')
    try:
        #Recovering the parsed options
        dirList = args.input_dirs
        keyField = int(args.key_field)
        separator = args.separator
        outputFile = args.output_file
        pattern = args.pattern
        l = args.ef_len
        lf = args.ef_len_field

        if args.no_header:
            seen_header = True
        else:
            seen_header = False

        IDList = []
        outdict = {} #key = keyField; value=list of fields to join from different files.

        #BUFFERING INPUT
        #Loop through all the files
        n=1
        for inputDir in dirList:

            logging.info("Accessing dir: %s" %inputDir)

            dictID = os.path.basename(inputDir.rstrip("/"))
            IDList.append(dictID)

            if l == False:
                quantFile = inputDir+"sm_quant/quant.sf"
                quantdict = readQuantFile(quantFile, seen_header, separator)
            else:
                lenFile = l[n]
                quantdict = (lenFile, lf, seen_header, separator)

            [statsFile] = glob.glob(inputDir+"outputs/*"+pattern+"*.stats")
            ''.join(statsFile)
            ribodict = readRiboFile(statsFile)

            outdict = getOutDict(outdict, ribodict, quantdict, dictID)

        n += 1

        #WRITING OUTPUT
        logging.info("Writing output to %s" %outputFile)
        printOutFiles(outdict, outputFile, IDList)
        logging.info("Process ended succesfully.")

    except BaseException as err:
         logging.error('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
         logging.error("%s" % err)
         sys.exit(1)

if __name__ == "__main__":
    main()


