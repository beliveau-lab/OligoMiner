#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# structureCheck.py
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
scriptName = 'structureCheck'

# Specify script version.
Version = '1.7'

# Import module for handling input arguments.
import argparse

# Import os module.
import os

# Import timeit module and record start time. This provides a rough estimate of
# the wall clock time it takes to run the script.
import timeit

# Import numpy module.
import numpy as np

# Import subprocess module.
import subprocess
from subprocess import PIPE, Popen

# Import Biopython modules.
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class StructureChecker:
    def __init__(self, inputFile, formConc, saltConc, NUPACKmat, threshVal,
                 Temp, IDval, reportVal, debugVal, metaVal, tempDir,
                 outNameVal, startTime):
        self.inputFile = inputFile
        self.formConc = formConc
        self.saltConc = saltConc
        self.NUPACKmat = NUPACKmat
        self.threshVal = threshVal
        self.Temp = Temp
        self.IDval = IDval
        self.reportVal = reportVal
        self.debugVal = debugVal
        self.metaVal = metaVal
        self.tempDir = tempDir
        self.outNameVal = outNameVal
        self.startTime = startTime

        # Calculate T to use with NUPACK.
        self.CorrTemp = 0.65 * self.formConc + self.Temp
        # Calculate salt to use with NUPACK.
        self.CorrSalt = self.saltConc * 0.001
        # Determine the stem of the input filename.
        self.fileName = str(self.inputFile).split('.')[0]


    def prob_check(self, seq1, struct1):
        """Runs the NUPACK 'prob' command."""
        randomInt = np.random.random_integers(0, 1000000)
        NUPACK_input = open('%s_%s_%0.0f_%d_prob_temp.in'
                            % (self.fileName, self.IDval, self.Temp, randomInt),
                            'w')
        NUPACK_input.write('%s\n%s\n' % (seq1, struct1))
        NUPACK_input.close()
        p = Popen(['prob', '-T', str(self.CorrTemp), '-sodium',
                   str(self.CorrSalt), '-material', str(self.NUPACKmat),
                   '%s_%s_%0.0f_%d_prob_temp' \
                   % (self.fileName, self.IDval,self.Temp, randomInt),
                   '>%s_%s_%0.0f_%d_temp_prob.txt' \
                   % (self.fileName, self.IDval, self.Temp, randomInt)],
                   stderr=PIPE, stdout=PIPE, bufsize=1)
        stdout_lines = p.stdout.readlines()
        prob_val = None
        for i in range(len(stdout_lines)):
            if stdout_lines[i] == '% Probability:\n':
                prob_val = float(stdout_lines[i + 1])
        if prob_val == None:
            print '***********************************************************'
            print 'NUPACK ERROR: could not run prob command with these inputs.'
            print '***********************************************************'
        os.remove('%s_%s_%0.0f_%d_prob_temp.in' \
                  % (self.fileName, self.IDval, self.Temp, randomInt))
        return prob_val


    def run(self):
        """Runs the structureChecker with the given parameters."""
        # Create a randomized directory to hold temp files.
        dirName = '%s_%d' \
                  % (self.fileName, np.random.random_integers(0, 1000000))
        if not os.path.exists('%s/%s' % (self.tempDir, dirName)):
            os.makedirs('%s/%s' % (self.tempDir, dirName))

        # Open input file for reading.
        with open(self.inputFile, 'r') as f:
            file_read = [line.strip() for line in f]

        # Create list to hold output.
        outList = []
        probList = []

        # Make list to hold Report info if desired.
        if self.reportVal is True:
            reportList = []

        # Iterate through probe file checking for predicted secondary structure.
        for i in range(0, len(file_read), 1):
            probeSeq = file_read[i].split('\t')[3]
            struct_in = '.' * len(probeSeq)
            p = self.prob_check(probeSeq, struct_in)
            probList.append(p)
            if p >= self.threshVal:
                outList.append(file_read[i])
                if self.reportVal is True:
                    reportList.append('Candidate probe at %s:%s-%s is '
                                      'predicted to have a linear structure '
                                      'with p=%0.4f > %0.4f at %dC in %d mM '
                                      'Na+ and %d%% formamide, added to '
                                      'output' \
                                      % (file_read[i].split('\t')[0],
                                         file_read[i].split('\t')[1],
                                         file_read[i].split('\t')[2],
                                         p, self.threshVal, self.Temp,
                                         self.saltConc, self.formConc))
                if self.debugVal is True:
                    print('Candidate probe at %s:%s-%s is predicted to have a '
                          'linear structure with p=%0.4f > %0.4f at %dC in %d '
                          'mM Na+ and %d%% formamide, added to output' \
                          % (file_read[i].split('\t')[0],
                             file_read[i].split('\t')[1],
                             file_read[i].split('\t')[2],
                             p, self.threshVal, self.Temp,
                             self.saltConc, self.formConc))
            else:
                if self.reportVal is True:
                    reportList.append('Candidate probe at %s:%s-%s is '
                                      'predicted to have a linear structure '
                                      'with p=%0.4f < %0.4f %dC in %d mM Na+ '
                                      'and %d%% formamide, filtered from '
                                      'output' \
                                      % (file_read[i].split('\t')[0],
                                         file_read[i].split('\t')[1],
                                         file_read[i].split('\t')[2],
                                         p, self.threshVal, self.Temp,
                                         self.saltConc, self.formConc))
                if self.debugVal is True:
                    print('Candidate probe at %s:%s-%s is predicted to have a '
                          'linear structure with p=%0.4f < %0.4f %dC in %d mM '
                          'Na+ and %d%% formamide, filtered from output' \
                          % (file_read[i].split('\t')[0],
                             file_read[i].split('\t')[1],
                             file_read[i].split('\t')[2],
                             p, self.threshVal, self.Temp,
                             self.saltConc, self.formConc))

        # Remove the temporary directory.
        os.removedirs('%s/%s' % (self.tempDir, dirName))


        # Determine the name of the output file.
        if self.outNameVal is None:
            if self.IDval is None:
                outName = '%s_sC' % self.fileName
            else:
                outName = '%s_%s_sC' % (self.fileName, self.IDval)
        else:
            outName = self.outNameVal

        # Create the output file.
        output = open('%s.bed' % outName, 'w')

        # Write the output file.
        output.write('\n'.join(outList))
        output.close()

        # Print info about the results to terminal.
        candsNum = len(file_read)
        cleanNum = len(outList)
        if 'rna' in self.NUPACKmat:
            print '********************************************************************************'
            print 'NUPACK WARNING: No salt corrections available for RNA.  Using 1 M Na and 0 M Mg.'
            print '********************************************************************************'
        print('structureCheck predicted that %d of %d / %0.4f%% candidate '
              'probes are predicted to have a linear structure with p>%0.4f at '
              '%dC in %d mM Na+ and %d%% formamide' \
              % (cleanNum, candsNum, float(cleanNum) / float(candsNum) * 100,
                 self.threshVal, self.Temp, self.saltConc, self.formConc))

        # Write meta information to a .txt file if desired.
        if self.metaVal is True:
            metaText = open('%s_outputClean_meta.txt' % outName, 'w')
            metaText.write('%s\t%f\t%s\t%0.2f\t%d\t%d' \
                           % (self.inputFile,
                              timeit.default_timer() - self.startTime,
                              Version, self.Temp, cleanNum, candsNum))
            metaText.close()

        # If desired, create report file and tabulate stats.
        if self.reportVal is True:
            a = np.asarray(probList)
            bin1 = ((0 < a) & (a <= .0001)).sum()
            bin2 = ((.0001 < a) & (a <= .001)).sum()
            bin3 = ((.001 < a) & (a <= .01)).sum()
            bin4 = ((.01 < a) & (a <= .1)).sum()
            bin5 = ((.1 < a) & (a <= 1)).sum()
            reportOut = open('%s_structureCheck_log.txt' % outName, 'w')
            reportList.insert(0, 'Results produced by %s %s' \
                                 % (scriptName, Version))
            reportList.insert(1, '-' * 100)
            reportList.insert(2, 'structureCheck predicted that %d of %d / '
                                 '%0.4f%% candidate probes are predicted to '
                                 'have a linear structure with p>%0.4f at %dC '
                                 'in %d mM Na+ and %d%% formamide' \
                                 % (cleanNum, candsNum,
                                    float(cleanNum) / float(candsNum) * 100,
                                    self.threshVal, self.Temp, self.saltConc,
                                    self.formConc))
            if len(a) == 0:
                reportList.insert(3, '-' * 100)
            else:
                reportList.insert(3, '%d of %d / %0.4f%% of probes have a '
                                     'predicted to have linear structures with '
                                     '0 < prob. <= .0001' \
                                     % (bin1, len(a),
                                        100 * float(bin1) / float(len(a))))
                reportList.insert(4, '%d of %d / %0.4f%% of probes have a '
                                     'predicted to have linear structures with '
                                     '.0001 < prob. <= .001'
                                     % (bin2, len(a),
                                        100 * float(bin2) / float(len(a))))
                reportList.insert(5, '%d of %d / %0.4f%% of probes have a '
                                     'predicted to have linear structures with '
                                     '.001 < prob. <= .01'
                                     % (bin3, len(a),
                                        100 * float(bin3) / float(len(a))))
                reportList.insert(6, '%d of %d / %0.4f%% of probes have a '
                                     'predicted to have linear structures with '
                                     '.01 < prob. <= .1'
                                     % (bin4, len(a),
                                        100 * float(bin4) / float(len(a))))
                reportList.insert(7, '%d of %d / %0.4f%% of probes have a '
                                     'predicted to have linear structures with '
                                     '.1 < prob. <= 1' \
                                     % (bin5, len(a),
                                        100 * float(bin5) / float(len(a))))
                reportList.insert(8, '-' * 100)
            reportOut.write('\n'.join(reportList))
            reportOut.close()

def runStructureChecker(inputFile, formConc, saltConc, NUPACKmat, threshVal,
                        Temp, IDval, reportVal, debugVal, metaVal, tempDir,
                        outNameVal, startTime):
    """Creates and runs an instance of a StructureChecker, which scans probes
    and evaluates their structures using NUPACK."""
    sc = StructureChecker(inputFile, formConc, saltConc, NUPACKmat, threshVal,
                          Temp, IDval, reportVal, debugVal, metaVal, tempDir,
                          outNameVal, startTime)
    sc.run()


def main():
    """Uses NUPACK to check structures of a given probe set according to
    specified conditions."""

    startTime = timeit.default_timer()

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .bed file containing probe sequences in the '
        'fourth tab-separated column. Also requires NUPACK to be installed and '
        'in your path. Returns a .bed file in the same format as the input '
        'file filtered of probes predicted to have more thermodynamically '
        'stable secondary structures or homodimer structure than the desired '
        'probe-target duplex by NUPACK' % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .bed file to check probes in')
    userInput.add_argument('-F', '--formamide', action='store', default=50,
                           type=float,
                           help='The percent formamide being used, default is '
                                '50')
    userInput.add_argument('-s', '--salt', action='store', default=390,
                           type=float,
                           help='The mM Na+ concentration, default is 390. '
                                'NOTE: NUPACK\'s allowable range is 50-1100 mM')
    userInput.add_argument('-m', '--material', action='store',
                           default='dna1998', type=str,
                           help='The NUPACK material setting, default is '
                                'dna1998')
    userInput.add_argument('-t', '--threshold', action='store', type=float,
                           default=0.1,
                           help='The probability threshold for a probe having '
                                'a linear structure to use for filtering, '
                                'default=0.1')
    userInput.add_argument('-T', '--hybTemp', action='store', default=47,
                           type=float,
                           help='The temperature at which you want to '
                                'hybridize your probes')
    userInput.add_argument('-I', '--ID', action='store', type=str, default=None,
                           help='Specify an ID to be associated with temporary '
                                'file names & the default output name. Can be '
                                'useful for parallelization. Null by default. '
                                'Will not be automatically included in output '
                                'file name if -o is flagged')
    userInput.add_argument('-R', '--Report', action='store_true', default=False,
                           help='Write a Report file detailing the results of '
                                'the secondary structure filtering. Off by '
                                'default. Note, selecting this option will '
                                'slow the script.')
    userInput.add_argument('-D', '--Debug', action='store_true', default=False,
                           help='The same as -Report, but prints info to '
                                'terminal instead of writing a log file. Off '
                                'by default')
    userInput.add_argument('-M', '--Meta', action='store_true', default=False,
                           help='Write a text file containing meta '
                                'information. Off by default. Reports input '
                                'file <tab> estimated runtime <tab> '
                                'structureCheck version ]<tab> hybridization '
                                'temperature used  <tab> number of input '
                                'probes passing secondary structure filter '
                                '<tab> number of input probes')
    userInput.add_argument('-d', '--temp', action='store', type=str,
                           default='/tmp',
                           help='the directory where temporary files will be a '
                                'written, default="/tmp"')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the name prefix of the output file')

    # Import user-specified command line values.
    args = userInput.parse_args()
    inputFile = args.file
    formConc = args.formamide
    saltConc = args.salt
    NUPACKmat = args.material
    threshVal = args.threshold
    Temp = args.hybTemp
    IDval = args.ID
    reportVal = args.Report
    debugVal = args.Debug
    metaVal = args.Meta
    tempDir = args.temp
    outNameVal = args.output

    # Run the structure checker.
    runStructureChecker(inputFile, formConc, saltConc, NUPACKmat, threshVal,
                        Temp, IDval, reportVal, debugVal, metaVal, tempDir,
                        outNameVal, startTime)

    # Print wall-clock runtime to terminal.
    print 'Program took %f seconds' % (timeit.default_timer() - startTime)


if __name__ == '__main__':
    main()
