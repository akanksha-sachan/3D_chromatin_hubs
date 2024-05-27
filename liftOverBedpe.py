import argparse
import os

#####################################################################################################
#                                                                                                   #
#                                            Functions                                              #
#                                                                                                   #
#####################################################################################################


def split_bedpe(bedpe, tmp1="tmp1.bed", tmp2="tmp2.bed", header="F", verbose="F"):
    """Define a function that splits up bedpe files into 2 temp files"""
    with open(tmp1, "w") as t1, open(tmp2, "w") as t2, open(bedpe, "r") as bedpein:
        if header == "T":
            headerline = bedpein.readline()
        name = 0
        for l in bedpein:
            name += 1
            e = l.strip().split("\t")
            t1.write("\t".join(e[0:3] + ["name" + str(name)] + [e[6]]) + "\n")
            t2.write("\t".join(e[3:6] + ["name" + str(name)] + [e[6]]) + "\n")


def do_lift_over(lift_over, chain, infile, verbose="F"):
    """Define a function that implements liftOver"""
    cmd = " ".join([lift_over, infile, chain, infile + ".success", infile + ".failure"])
    if verbose == "T":
        print(cmd)
    os.system(cmd)


def merge_lift_over(f1, f2, outfile, verbose="F"):
    """Define a function that merges liftOver"""
    readdict = {}
    with open(f1, "r") as f:
        for l in f:
            e = l.strip().split("\t")
            readdict[e[3]] = e[:6]

    with open(f2, "r") as f, open(outfile, "w") as o:
        for l in f:
            e = l.strip().split("\t")
            if e[3] in readdict:
                r1 = readdict[e[3]]
                r2 = e
                o.write("\t".join(r1[:3] + r2[:3] + [r2[4]]) + "\n")


#####################################################################################################
#                                                                                                   #
#                                        Parse arguments                                            #
#                                                                                                   #
#####################################################################################################

parser = argparse.ArgumentParser(
    description="Wrapper for liftOver to accommodate bedpe files"
)

# required arguments
parser.add_argument("--lift", dest="liftOver", help="Path to liftOver")
parser.add_argument("--chain", dest="chain", help="(e.g. hg19ToHg18.over.chain)")
parser.add_argument("--i", dest="infile", help="Input file in bedpe format")
parser.add_argument("--o", dest="outfile", help="Output file")
parser.add_argument("--v", dest="verbose", help="Verbose", default="F")
parser.add_argument(
    "--h", dest="header", help="T / F if there is a header line", default="F"
)

# parse arguments
args = parser.parse_args()

# read in args
LO = args.liftOver
chain = args.chain
bedpeIN = args.infile
bedpeOUT = args.outfile
TMP1 = "tmp1.bed"
TMP2 = "tmp2.bed"
header = args.header
verbose = args.verbose

#####################################################################################################
#                                                                                                   #
#                                        Run the Code                                               #
#                                                                                                   #
#####################################################################################################

# break up the files
split_bedpe(bedpeIN, TMP1, TMP2, header, verbose)

# perform liftOver
do_lift_over(LO, chain, TMP1, verbose)
do_lift_over(LO, chain, TMP2, verbose)

# merge liftOvered files
merge_lift_over(TMP1 + ".success", TMP2 + ".success", bedpeOUT, verbose)

# remove tmp files
os.remove(TMP1)
os.remove(TMP2)
os.remove(TMP1 + ".success")
os.remove(TMP1 + ".failure")
os.remove(TMP2 + ".success")
os.remove(TMP2 + ".failure")
