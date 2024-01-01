#!/usr/bin/env python3
##########################################################
# hmmscan parser for dbCAN meta server
#
# Based off the hmmscan parser used in the dbCAN server,
# written by Dr. Yin
#
# Written by Tanner Yohe under the supervision
# of Dr. Yin in the YinLab at NIU.
#
# Updated by Le Huang from tips the contributor WATSON Mick <mick.watson@roslin.ed.ac.uk>,
# Thank you!
#
# Modified by Alex Fraser to have a run() method that can be called and returns data for better integration with other
# scripts. This script also retains the ability to be called from shell and output to pipe redirection.
# This file had to be renamed from "hmmscan-parser.py" to "hmmscan_parser.py" because of python module import conventions.
# Modified on 07/06/22
#
# INPUT
# python hmmscan-parser-dbCANmeta.py [inputFile] [eval] [coverage]
# eval and coverage are optional, inputFile is required
# -updating info:
# -adds pid for every subprocess to make codes robust.
# Last updated: 1/10/19
###########################################################

import os
import sys
from subprocess import call


def run(input_file, eval_num=1e-15, coverage=0.35, verbose=False):
    """
    Executes a genomic data processing pipeline on a given input file.

    This function processes an input genomic file, applying various transformations and filtering criteria based on evaluation number and coverage. It first executes a series of shell commands to manipulate the data, then reads the processed data, applies additional filtering, and returns the result as a string.

    Parameters
    ----------
        input_file (str): The path to the input file containing genomic data.
        eval_num (float, optional): Evaluation number threshold for data filtering, defaults to 1e-15.
        coverage (float, optional): Coverage threshold for data filtering, defaults to 0.35.
        verbose (bool, optional): If True, outputs additional processing information, defaults to False.

    Returns
    -------
        str: A string containing processed and filtered genomic data. Each line in the string represents a row of data that meets the specified evaluation number and coverage thresholds.

    Side Effects:
        - Creates a temporary file named 'temp.<process_id>' for intermediate processing.
        - Deletes the temporary file after processing.
        - Outputs processing details to the standard output if verbose is True.

    Notes
    -----
        - The function relies on external shell commands and Perl for initial data processing.
        - Assumes a specific format for the input data, including expected columns and data types.
    """
    tmpfile = "temp." + str(os.getpid())

    call(
        "cat "
        + input_file
        + "  | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_.\"\n\";}}' > "
        + tmpfile,
        shell=True,
    )

    output = ""
    with open(tmpfile) as f:
        for line in f:
            row = line.rstrip().split("\t")
            row.append(float(int(row[6]) - int(row[5])) / int(row[1]))
            if float(row[4]) <= float(eval_num) and float(row[-1]) >= float(coverage):
                if verbose:
                    print("\t".join([str(x) for x in row]))
                output += "\t".join([str(x) for x in row]) + "\n"
    call(["rm", tmpfile])

    return output


if __name__ == "__main__":
    if len(sys.argv) > 3:
        file = sys.argv[1]
        eval_arg = float(sys.argv[2])
        coverage_arg = float(sys.argv[3])
        run(file, eval_arg, coverage_arg, verbose=True)
    if len(sys.argv) > 1:
        file = sys.argv[1]
        run(file, verbose=True)
    else:
        print("Please give a hmmscan output file as the first command")
        exit()
