#!/Softs/tagc/Python-2.7.9/binaries/bin/python2.7

"""
Description:  Compute fragment length distribution from a PE bam file.
Developer: Denis Puthier
Last modifications: Wed Jun 29 12:05:57 CEST 2016
Version: {v}
"""

import re  # Import regular expression
import sys
import argparse
import pysam
import os.path
from collections import defaultdict
__version__ = 0.1



def make_parser():
    """A function that returns the program parser.
    """

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp_main = parser.add_argument_group('Arguments')

    parser_grp_main.add_argument

    parser_grp_main.add_argument("inputfiles",
                                  type=argparse.FileType("r"),
                                  help="The bam files.",
                                  nargs='+')

    parser_grp_main.add_argument("-o",
                                  "--out-dir",
                                  type=str,
                                  help="The output directory",
                                  default=".",
                                  metavar="\b",
                                  required=False)

    parser_grp_main.add_argument("-l",
                                  "--label",
                                  type=str,
                                  help="Labels for the bam files. To be used as alternative header. Comma-separated list.",
                                  default=None,
                                  metavar="\b",
                                  required=False)

    parser_grp_main.add_argument("-m",
                                  "--min",
                                  type=int,
                                  help="The min value for distribution. Lower values won't be stored.",
                                  default=1,
                                  metavar="\b",
                                  required=False)

    parser_grp_main.add_argument("-M",
                                  "--max",
                                  type=int,
                                  help="The max value for distribution. Higher values won't be stored",
                                  default=400,
                                  metavar="\b",
                                  required=False)

    parser_grp_main.add_argument("-q",
                                  "--mapq",
                                  type=int,
                                  help=" Only keep read with mapping quality equal or greater than q.",
                                  default=0,
                                  metavar="\b",
                                  required=False)


    return parser


def bam_frag_len_pe(inputfiles=None,
                    out_dir=".",
                    min=1,
                    max=400,
                    mapq=10,
                    label=None,
                    tmp_dir=None,
                    verbosity=None,
                    keep=False,
                    logger_file=None):

    # check labels
    if label is not None:
        label_list = label.split(",")
        if len(label_list) != len(inputfiles):
            sys.stderr.write("Error: the number of bam files and labels should be the same.\n")
            sys.exit()
    else:
        label_list = [x.name for x in inputfiles]

    label_dict = dict()
    # Dict of labels
    if label is not None:
        for i in range(len(label_list)):
            label_dict[inputfiles[i].name] = label_list[i]
    else:
        for i in range(len(label_list)):
            label_dict[inputfiles[i].name] = inputfiles[i].name

    # Dict of fragment length
    frag_len = defaultdict(dict)




    # This dict stores the results

    for lb in label_dict:
        lb = label_dict[lb]
        for i in range(min, max + 1):
            frag_len[lb][i] = 0


    nb_read_1 = 0


    for bam in inputfiles:

        if not os.path.isfile(bam.name + ".bai"):
            bam_fn_test = re.sub(".bam", "", bam.name, re.IGNORECASE)
            if not os.path.isfile(bam_fn_test + ".bai"):
                sys.stderr.write("Can't find any bam index "
                                 "for file : " + bam.name + "\n")
                sys.exit()
        bam.close()

        samfile = pysam.Samfile(bam.name, "rb")

        sam_records = samfile.fetch()

        for read in sam_records:
            if read.is_read1:
                if read.is_paired:
                    if not read.mate_is_unmapped:
                        if read.is_proper_pair:
                            if not read.is_secondary:
                                if read.mapq >= mapq:
                                    cur_len = abs(int(read.tlen))
                                    if cur_len >= min and cur_len <= max:
                                        frag_len[label_dict[bam.name]][cur_len] += 1

    print("\t".join(["SIZE"] + label_list))

    for i in range(min, max + 1):
        tmp = [str(i)]
        for lb in label_list:
            tmp += [str(frag_len[lb][i])]

        print("\t".join(tmp))


if __name__ == '__main__':

    args = make_parser()
    args = args.parse_args()
    args = dict(args.__dict__)
    bam_frag_len_pe(**args)

else:

    from missgenomics.cmd_object import CmdObject

    cmd = CmdObject(name="bam_frag_len_pe",
                    message="Compute fragment length distribution from a PE bam file.",
                    parser=make_parser(),
                    fun=bam_frag_len_pe,
                    desc=__doc__.format(v=__version__))





