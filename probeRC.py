#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# probeRC.py
#
# (c) 2017 Molecular Systems Lab
#
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
scriptName = 'probeRC'

# Specify script version.
Version = '1.7.1'

# Import module for handling input arguments.
import argparse

# Import Biopython modules.
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def createRCs(inputFile, outNameVal):
    """Creates a .bed file with the reverse complements of the given set of
    sequences."""

    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    # Open input file for reading.
    with open(inputFile, 'r') as f:
        file_read = [line.strip() for line in f]

    # Create list to hold output.
    outList = []

    # Parse out probe info, flip sequence to RC, and write to output list.
    for i in range(0, len(file_read), 1):
        chrom = file_read[i].split('\t')[0]
        start = file_read[i].split('\t')[1]
        stop = file_read[i].split('\t')[2]
        probeSeq = file_read[i].split('\t')[3]
        RevSeq = Seq(probeSeq, IUPAC.unambiguous_dna).reverse_complement()
        Tm = file_read[i].split('\t')[4]
        outList.append('%s\t%s\t%s\t%s\t%s' % (chrom, start, stop, RevSeq, Tm))

    # Determine the name of the output file.
    if outNameVal is None:
        outName = '%s_RC' % fileName
    else:
        outName = outNameVal

    # Create the output file.
    output = open('%s.bed' % outName, 'w')

    # Write the output file
    output.write('\n'.join(outList))
    output.close()


def main():
    """Produces a .bed file with the reverse complements of the given
    chromosome sequences."""

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .bed file with first four columns in the '
        'format chromosome <tab> start <tab> stop <tab> sequence <tab> Tm such '
	    ' as the .bed files produced by outputClean. Returns a .bed file that is '
        'identical to the input file except that the probe sequences have been '
        'replaced with their reverse complements.' % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .bed file containing the probe '
                                    'sequences to take the reverse complements '
                                    'of')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str, help='Specify the name prefix of the '
                                          'output file')

    # Import user-specified command line values
    args = userInput.parse_args()
    inputFile = args.file
    outNameVal = args.output

    createRCs(inputFile, outNameVal)

if __name__ == '__main__':
    main()
