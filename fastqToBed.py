#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# fastqToBed.py
# 
# (c) 2016 Molecular Systems Lab
# Wyss Institute for Biologically-Inspired Engineering
# Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# --------------------------------------------------------------------------

# Specific script name.
scriptName = 'fastqToBed'

# Specify script version.
Version = '1.7'

# Import module for handling input arguments.
import argparse

# Import Biopython mt module.
from Bio.SeqUtils import MeltingTemp as mt

def probeTm(seq1, saltConc, formConc):
    """Calculates the melting temperature of a given sequence under the
    specified salt and formamide conditions."""

    tmval = float(('%0.2f' % mt.Tm_NN(seq1, Na=saltConc)))
    fcorrected = ('%0.2f' % mt.chem_correction(tmval, fmd=formConc))
    return fcorrected


def convertFastqToBed(inputFile, saltConc, formConc, outNameVal):
    """Converts a given .fastq file to a .bed file."""

    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    # Open input file for reading.
    with open(inputFile, 'r') as f:
        file_read = [line.strip() for line in f]

    # Create list to hold output.
    outList = []

    # Parse .fastq and extract probe information.
    for i in range(0, len(file_read), 4):
        chrom = file_read[i].split(':')[0].split('@')[1]
        start = file_read[i].split(':')[1].split('-')[0]
        stop = file_read[i].split(':')[1].split('-')[1]
        seq = file_read[i+1]
        Tm = probeTm(seq, saltConc, formConc)
        outList.append('%s\t%s\t%s\t%s\t%s' % (chrom, start, stop, seq, Tm))

    # Determine the name of the output file.
    if outNameVal is None:
        outName = fileName
    else:
        outName = outNameVal

    # Create the output file.
    output = open('%s.bed' % outName, 'w')

    # Write the output file.
    output.write('\n'.join(outList))
    output.close()


def main():
    """Converts a .bed file to a .fastq file, taking the filenames as
    command line arguments."""
    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .fastq file containing chr, start, stop '
        'information in the sequence ID field for each entry in the format '
        '@chr:start-stop such as the .fastq files produced by blockParse. '
        'Returns a .bed file in the format chromosome <tab> start <tab> stop '
        '<tab> sequence <tab> Tm. Calculates the Tm of each probe based on -F '
        'and -s parameters.' % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .bed file to convert to .fastq')
    userInput.add_argument('-s', '--salt', action='store', default=390,
                           type=int,
                           help='The mM Na+ concentration, default is 390')
    userInput.add_argument('-F', '--formamide', action='store', default=50,
                           type=float,
                           help='The percent formamide being used, default is '
                                '50')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the name prefix of the output file')

    # Import user-specified command line values
    args = userInput.parse_args()
    inputFile = args.file
    saltConc = args.salt
    formConc = args.formamide
    outNameVal = args.output

    convertFastqToBed(inputFile, saltConc, formConc, outNameVal)


if __name__ == '__main__':
    main()
