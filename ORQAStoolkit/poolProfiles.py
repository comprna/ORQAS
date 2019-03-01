import sys
import glob
import os
import logging
from math import log, exp, nan
from argparse import ArgumentParser, RawTextHelpFormatter, SUPPRESS

description = \
    "Description:\n\n" + \
    "This script takes ribosome profiles per codon or base for multiple samples and pools them together\n"

parser = ArgumentParser(description = description, formatter_class=RawTextHelpFormatter, add_help=False)
parser.add_argument("-i", "--input-dirs", nargs="+", required = True,
                    help = "spaced separated list of files to pool")
parser.add_argument("-c", "--codon", action="store_true",
                    help = "aggregated CDS file")
parser.add_argument("-b", "--base", action="store_true",
                    help = "aggregated CDS file")
parser.add_argument("--cds", required = True,
                    help = "aggregated CDS file")
parser.add_argument("-o", "--output-pref", required=True,
                    help = "prefix for the output file.")
parser.add_argument("--pooled", action="store_true", default=False,
                    help = "use it to obtain a single file with all the pooled samples")




def readProfFile(profFile, riboDict):
    with open(profFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(":")
            fields = line[1].strip()
            if line[0] == "tid":
                keyname = fields.split(".")
                if keyname[0] not in riboDict:
                    riboDict.setdefault(keyname[0], {"ribo":[], "norm_ribo":[], "mrna":[]})
            elif line[0] == "ribo profile":
                values = [float(x) for x in fields.split(" ")]
                riboDict[keyname[0]]["ribo"] = values

            elif line[0] == "normalized ribo profile":
                values = [float(x) for x in fields.split(" ")]
                riboDict[keyname[0]]["norm_ribo"] = values

            elif line[0] == "mRNA profile":
                values = [float(x) for x in fields.split(" ")]
                riboDict[keyname[0]]["mrna"] = values

        logging.info("File %s closed." % profFile)
        return(riboDict)

# def fillDictionary(profDict, num):
#     for k,v in profDict.items():
#         for k1, v1 in v.items():
#             if len(v1) < num:
#                 rep = num - len(v1)
#                 zeros = tuple([0]*len(v1[0]))
#                 for i in range(1, rep+1):
#                     profDict[k][k1].append(zeros)
#     return(profDict)


def sumvalues(quantDict, field, ENSTlist):
    LoT = []
    for id in ENSTlist:
        try:
            LoT.append(quantDict[id][field])
        except:
            continue
    CDSvalue = [sum(sample) for sample in zip(*LoT)] #unpack and zip
    return(CDSvalue)


def mergeCDS(profDict, cdsFile):
    with open(cdsFile, 'r') as f:
        dictionary = {}
        for l in f:
            line = l.rstrip('\n').split("\t")
            ENST = line[1]
            ENSTlist = ENST.split(":")
            for k, v in profDict.items():
                dictionary.setdefault(k, {})
                for k1, v1 in v.items():
                    if len(ENSTlist) > 1:
                        dictionary[k][k1] = sumvalues(profDict, k1, ENSTlist)
                    else:
                        try:
                            dictionary[k][k1] = profDict[ENSTlist[0]]
                        except:
                            continue
    logging.info("File %s closed." %CDSinfo)
    return(dictionary)

#def poolSamples:
#    dafhds


def main():
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - INFO: %(message)s')
    try:
        #Recovering the parsed options
        dirList = args.input_dirs
        outputPref = args.output_pref
        pooled = args.pooled
        cdsFile = args.cds
        c = args.codon
        b = args.base
        if c == False and b == False:
            logging.error("You must specify base (-b) or codon (-c) profiles")

        #Loop through all the files
        profDict = {}
        inputDir = dirList[0]
        logging.info("Accessing dir: %s" %inputDir)
        sampleID = os.path.basename(inputDir.rstrip("/"))
        if c:
            logging.info("Getting output codon profiles")
            profFile = glob.glob(inputDir+"outputs/*.codon")[0]
            profDict = readProfFile(profFile, profDict)


        #filledprofDict = fillDictionary(profDict, len(IDList))
        aggprofDict = mergeCDS(profDict, cdsFile)
        writeSingle(aggprofDict, outputPref)


    except BaseException as err:
         logging.error('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
         logging.error("%s" % err)
         sys.exit(1)
