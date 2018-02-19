#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# probeTm.py
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
scriptName = 'probeTm'

# Specify script version.
Version = '1.7'

# Import module for handling input arguments.
import argparse

# Import Biopython mt module.
from Bio.SeqUtils import MeltingTemp as mt

# Import regex library.
import re

def probeTm(seq1, conc1, conc2, saltConc, formConc):
    """Calculates the Tm of a given sequence."""
    tmval = float(('%0.2f' \
                   % mt.Tm_NN(seq1, Na=saltConc, dnac1=conc1, dnac2=conc2)))
    fcorrected = ('%0.2f' % mt.chem_correction(tmval, fmd=formConc))
    return fcorrected


def getTm(inputFile, saltConc, formConc, conc1, conc2, inputSeqVal, outNameVal):
    """Determines the melting temperatures of a given probe set."""

    # Iterate through input file, if present, and calculate the Tm of all input
    # sequences.
    if inputFile is not None:
        with open(inputFile, 'r') as f:
            file_read = [line.strip() for line in f]

        # Create list to hold output.
        outList = []

        # Iterate through probe file checking for predicted secondary structure.
        for i in range(0, len(file_read), 1):
            probeSeq = file_read[i].split('\t')[1]

            # Skip any sequences containing 'N' bases as these cannot be
            # processed.
            if len(re.findall('N', probeSeq, re.I)) > 0:
                print '\'N\' base(s) found in the sequence in row %d of the ' \
                      'input file...skipping this sequence' % i

            # Calculate Tm of all sequences not containing 'N' bases, add to
            # output list as new column.
            else:
                probeTmVal = probeTm(probeSeq, conc1, conc2, saltConc, formConc)
                outList.append(file_read[i] + '\t' + probeTmVal)

            # Determine the name of the output file.
            if outNameVal is None:
                # Determine the stem of the input filename.
                fileName = inputFile.split('.')[0]
                # Create standard output filename.
                outName = '%s_tm' % fileName
            else:
                # Or use user-specified filename.
                outName = outNameVal

            # Create the output file.
            output = open('%s.txt' % outName, 'w')

            # Write the output file.
            output.write('\n'.join(outList))
            output.close()

    # If no file is provided, get sequence from stdin or user input.
    else:
        # Take input sequence from stdin if -i is flagged.
        if inputSeqVal is not None:
            probeSeq = inputSeqVal

        # Prompt user input if no input file is present and '-i' is not flagged.
        else:
            probeSeq = raw_input('Please input your sequence: ')

        # Check input sequence for the presence of 'N' bases and alert
        # user if any are found.
        if len(re.findall('N', probeSeq, re.I)) > 0:
            print '\'N\' base(s) found in the sequence ... Tm calculation ' \
                  'cannot be performed'

        # Print Tm value of input sequence to terminal / stdout.
        else:
            print probeTm(probeSeq, conc1, conc2, saltConc, formConc)


def main():
    """Determines the melting temperatures of given sequences, provided either
    as a commandline argument are through stdin."""

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a two column input file in the format: '
        'sequence ID <tab> sequence. Returns a file in the format sequence ID '
        '<tab> sequence <tab> sequence Tm. Will prompt user for input if no '
        'input sequences are provided.' % (scriptName, Version))
    userInput.add_argument('-f', '--file', action='store',
                           help='The file to containing the sequences that Tm '
                                'calculation will be performed on. Providing a '
                                'file will override the \'-i\' flag.')
    userInput.add_argument('-s', '--salt', action='store', default=390,
                           type=int,
                           help='The mM Na+ concentration, default is 390')
    userInput.add_argument('-F', '--formamide', action='store', default=50,
                           type=float,
                           help='The percent formamide being used, default is '
                                '50')
    userInput.add_argument('-c', '--dnac1', action='store', default=25,
                           type=float,
                           help='Concentration of higher concentration strand '
                                '[nM] -typically the probe- to use for '
                                'thermodynamic calculations. Default is 25')
    userInput.add_argument('-C', '--dnac2', action='store', default=25,
                           type=float,
                           help='Concentration of lower concentration strand '
                                '[nM] -typically the target- to use for '
                                'thermodynamic calculations. Default is 25')
    userInput.add_argument('-i', '--inputSeq', action='store', default=None,
                           help='Use this to input a sequence directly on the '
                                'command line /stdin instead of providing an '
                                'in input file. User will be prompted for '
                                'input if no sequence is provided. Will print '
                                'result to terminal / stdout.')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the name prefix of the output file')

    # Import user-specified command line values.
    args = userInput.parse_args()
    inputFile = args.file
    saltConc = args.salt
    formConc = args.formamide
    conc1 = args.dnac1
    conc2 = args.dnac2
    inputSeqVal = args.inputSeq
    outNameVal = args.output

    # Assign concentration variables based on magnitude.
    if args.dnac1 >= args.dnac2:
        conc1 = args.dnac1
        conc2 = args.dnac2

    else:
        conc1 = args.dnac2
        conc2 = args.dnac1

    getTm(inputFile, saltConc, formConc, conc1, conc2, inputSeqVal, outNameVal)

if __name__ == '__main__':
    main()
