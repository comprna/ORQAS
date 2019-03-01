import ribomaptoOPM as OPMcalculator
import validate_iso as validateORF
import TXTtoCDS
import aggregateCDS
import poolProfiles
import logging
import argparse
import sys


description = "Description:\n" + \
              "riboTools is developed in order to easily process Ribomap output. \n" \

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers()

# TXTtoCDS parser
TXTtoCDSSubparser = subparsers.add_parser(
    "TXTtoCDS", parents=[TXTtoCDS.parser],
    help="Get an index of the CDS that share the same sequence")
TXTtoCDSSubparser.set_defaults(which="TXTtoCDS")

# TPMcalculator parser
OPMcalculatorSubparser = subparsers.add_parser(
    "OPMcalculator", parents=[OPMcalculator.parser],
    help="Obtain OPM values using the ribosome abundance and the effective length from salmon quant")
OPMcalculatorSubparser.set_defaults(which="OPMcalculator")

# aggregateCDS parser
aggregateCDSSubparser = subparsers.add_parser(
     "aggregateCDS", parents=[aggregateCDS.parser],
     help="Obtain aggregated CDS TPM, OPM or count for transcripts encoding the same CDS if an annotation is provided")
aggregateCDSSubparser.set_defaults(which="aggregateCDS")

# poolProfiles parser
poolProfilesSubparser = subparsers.add_parser(
     "poolProfiles", parents=[poolProfiles.parser],
     help="Obtain pooled profiles (codon or base) by replicates to calculate coverage/periodicity")
poolProfilesSubparser.set_defaults(which="poolProfiles")

# validate_iso parser
validateORFSubparser = subparsers.add_parser(
     "validateORF", parents=[validateORF.parser],
     help="Get uniformity (PME) and periodicity values for each ORF")
validateORFSubparser.set_defaults(which="validateORF")


def main():
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    try:
        args = parser.parse_args()
        if args.which == "TXTtoCDS":
            TXTtoCDS.parser = parser  # Setting the module aparser
            TXTtoCDS.main()
        if args.which == "OPMcalculator":
            OPMcalculator.parser = parser  # Setting the module aparser
            OPMcalculator.main()
        elif args.which == "aggregateCDS":
            aggregateCDS.parser = parser  # Setting the module parser
            aggregateCDS.main()
        elif args.which == "poolProfiles":
            poolProfiles.parser = parser  # Setting the module parser
            poolProfiles.main()
        elif args.which == "validateORF":
            validateORF.parser = parser  # Setting the module parser
            validateORF.main()

    except BaseException as err:
        logging.error('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
        logging.error("%s" % err)
        sys.exit(1)


if __name__ == "__main__":
    main()
