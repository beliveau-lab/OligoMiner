#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# bedToFastq.py
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
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# --------------------------------------------------------------------------

# Specific script name.
scriptName = 'bedToFastq'

# Specify script version.
Version = '1.7'

# Import module for handling input arguments.
import argparse

def convertBedToFastq(inputFile, outNameVal):
    """Converts a .bed file to a .fastq file."""

    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    # Open input file for reading.
    with open(inputFile, 'r') as f:
        file_read = [line.strip() for line in f]

    # Create list to hold output.
    outList = []

    # A list to hold arbitrary quality scores for each base in the candidate
    # probe.
    quals = ['~' * len(file_read[i].split('\t')[3]) \
             for i in range(len(file_read))]

    # Parse out probe information and write into .fastq format.
    for i in range(0, len(file_read), 1):
        chrom = file_read[i].split('\t')[0]
        start = file_read[i].split('\t')[1]
        stop = file_read[i].split('\t')[2]
        probeSeq = file_read[i].split('\t')[3]
        outList.append('@%s:%s-%s\n%s\n+\n%s' \
                       % (chrom, start, stop, probeSeq, quals[i]))

    # Determine the name of the output file.
    if outNameVal is None:
        outName = fileName
    else:
        outName = outNameVal

    # Create the output file.
    output = open('%s.fastq' % outName, 'w')

    # Write the output file
    output.write('\n'.join(outList))
    output.close()


def main():
    """Converts a .bed file to a .fastq file, taking the filenames as
    command line arguments."""

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .bed file with first four columns in the '
        'format chromosome <tab> start <tab> stop <tab> sequence such as the '
        '.bed files produced by outputClean. Returns a .fastq file.' \
        % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .bed file to convert to .fastq')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the name prefix of the output file')

    # Import user-specified command line values
    args = userInput.parse_args()
    inputFile = args.file
    outNameVal = args.output

    convertBedToFastq(inputFile, outNameVal)


if __name__ == '__main__':
    main()
