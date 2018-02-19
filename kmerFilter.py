#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# kmerFilter.py
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
scriptName = 'kmerFilter'

# Specify script version.
Version = '1.7'

# Import module for handling input arguments.
import argparse

# Import timeit module and record start time. This provides a rough estimate of
# the wall clock time it takes to run the script.
import timeit

# Import numpy module.
import numpy as np

# Import subprocess module.
import subprocess

def runFilter(inputFile, outNameVal, merLengthVal, jellyfishFile, kVal, IDval,
              reportVal, debugVal, metaVal, startTime):
    """Runs Jellyfish to screen probes sequence from a .bed file for high
    abundance k-mers."""

    # Determine the stem of the input filename.
    fileName = inputFile.split('.')[0]

    # Make list to hold Report info if desired.
    if reportVal is True:
      reportList = []
      failVals = []

    # Open input file for reading.
    with open(inputFile, 'r') as f:
      file_read = [line.strip() for line in f]

    # Make list to hold the sequences of the imported probes.
    seqList = []

    # Parse out sequence info from each probe in the input .bed
    # and add to the newly created list.
    for i in range(0, len(file_read), 1):
      seqList.append('>' + '\n' + file_read[i].split('\t')[3])

    # Generate a randomized number to add to the temporary input/output files
    # required to run Jellyfish.
    randomInt = np.random.random_integers(0, 1000000)

    # Create a temporary input .fa file for jellyfish.
    jfFasta = open('%s_%d_%d_%d_temp.fa' \
                   % (fileName, merLengthVal, kVal, randomInt), 'w')

    # Write probe sequences to the .fa file.
    jfFasta.write('\n'.join(seqList))
    jfFasta.close()

    # Create index of kmer positions within each probe to reference later.
    indexList = []
    for i in range(0, len(seqList), 1):
      probeLength = len(seqList[i]) - 2
      merWindows = probeLength - merLengthVal + 1
      for j in range(0, merWindows, 1):
          indexList.append('%d' % i)

    # Call jellyfish and count kmers in the probe sequences.
    subprocess.call(['jellyfish', 'query', '%s' % jellyfishFile, '-s',
                     '%s_%d_%d_%d_temp.fa' % (fileName, merLengthVal,
                                              kVal, randomInt),
                     '-o', '%s_%d_%d_%d_temp.txt' % (fileName, merLengthVal,
                                                     kVal, randomInt)],
                    stderr=None, shell=False)

    # Open jellyfish output.
    with open('%s_%d_%d_%d_temp.txt' \
              % (fileName, merLengthVal, kVal, randomInt), 'r') as jfData:
      jf_read = [line.strip() for line in jfData]

    # Iterate through jellyfish output, checking each probe for the presence of
    # kmers occurring too frequently.
    excludeList = []
    for i in range(0, len(jf_read), 1):
      countVal = int(jf_read[i].split(' ')[1])
      if countVal >= kVal:
          if len(excludeList) > 0 and excludeList[-1] == int(indexList[i]):
              pass
          else:
              excludeList.append(int(indexList[i]))
              if reportVal or debugVal is True:
                  rChrom = file_read[int(indexList[i])].split('\t')[0]
                  rStart = file_read[int(indexList[i])].split('\t')[1]
                  rStop = file_read[int(indexList[i])].split('\t')[2]
                  if reportVal is True:
                      reportList.append('Candidate probe at %s:%s-%s filtered '
                                        'due to the presence of a %smer '
                                        'occurring %d times' \
                                        % (rChrom, rStart, rStop,merLengthVal,
                                           countVal))
                      failVals.append(countVal)
                  if debugVal is True:
                      print('Candidate probe at %s:%s-%s filtered due to the '
                            'presence of a %smer occurring %d times' \
                            % (rChrom, rStart, rStop, merLengthVal, countVal))

    # Convert exclude list to set.
    excludeSet = set()
    for x in excludeList:
      excludeSet.add(x)

    # Create a list to hold output.
    outList = []

    # Add probes passing kmer filter threshold to output list.
    for i in range(0, len(file_read), 1):
      if i not in excludeSet:
          outList.append(file_read[i])
          if reportVal or debugVal is True:
              rChrom = file_read[i].split('\t')[0]
              rStart = file_read[i].split('\t')[1]
              rStop = file_read[i].split('\t')[2]
              if reportVal is True:
                  reportList.append('Candidate probe at %s:%s-%s passed  %smer '
                                    'filtering using an occurrence threshold '
                                    'of %d, added to output' \
                                    % (rChrom, rStart, rStop, merLengthVal,
                                       kVal))
          if debugVal is True:
              print('Candidate probe at %s:%s-%s passed  %smer filtering '
                    'using an occurrence threshold of %d, added to output'
                    % (rChrom, rStart, rStop, merLengthVal, kVal))

    # Remove temporary files.
    subprocess.call(['rm', '%s_%d_%d_%d_temp.fa' \
                           % (fileName, merLengthVal, kVal, randomInt)],
                    stderr=None, shell=False)
    subprocess.call(['rm', '%s_%d_%d_%d_temp.txt' \
                           % (fileName, merLengthVal, kVal, randomInt)],
                    stderr=None, shell=False)

    # Determine the name of the output file.
    if outNameVal is None:
      outName = '%s_%d_%d' % (fileName, merLengthVal, kVal)
    else:
      outName = outNameVal

    # Create the output file.
    output = open('%s.bed' % outName, 'w')

    # Write the output file.
    output.write('\n'.join(outList))
    output.close()

    # Print info about the results to terminal.
    candsNum = len(file_read)
    cleanNum = len(outList)
    print 'kmerFilter identified %d of %d candidate probes / %0.4f%% as ' \
        'containing only %dmers occurring < %d times' \
        % (cleanNum, candsNum, float(cleanNum)/float(candsNum) * 100,
           merLengthVal, kVal)

    # Write meta information to a .txt file if desired.
    if metaVal is True:
      metaText = open('%s_kmerFilter_meta.txt' % outName, 'w')
      metaText.write('%s\t%f\t%s\t%d\t%d\t%d\t%d'
                     % (inputFile, timeit.default_timer() - startTime,
                        Version, kVal, merLengthVal, cleanNum, candsNum))
      metaText.close()

    # If desired, create report file and tabulate stats.
    if reportVal is True:
      reportList.sort(key=lambda x: x.split(':')[1])
      a = np.asarray(failVals)
      bin1 = ((1 < a) & (a < 10)).sum()
      bin2 = ((11 < a) & (a < 20)).sum()
      bin3 = ((21 < a) & (a < 50)).sum()
      bin4 = ((51 < a) & (a < 100)).sum()
      bin5 = (101 < a).sum()
      reportOut = open('%s_kmerFilter_log.txt' % outName, 'w')
      reportList.insert(0, 'Results produced by %s %s' % (scriptName, Version))
      reportList.insert(1, '-' * 100)
      reportList.insert(2, 'kmerFilter identified %d of %d / %0.4f%% candidate '
                           'probes as containing only %dmers occurring < %d '
                           'times' % (cleanNum, candsNum,
                                      float(cleanNum) / float(candsNum) * 100,
                                      merLengthVal, kVal))
      if len(a) == 0:
          reportList.insert(3, '-' * 100)
      else:
          reportList.insert(3, '%d of %d / %0.4f%% of failures were triggered '
                               'by the presence of a %dmer occurring between 2 '
                               'and 9 times' \
                               % (bin1, len(a),
                                  100 * float(bin1) / float(len(a)),
                                  merLengthVal))
          reportList.insert(4, '%d of %d / %0.4f%% of failures were triggered '
                               'by the presence of a %dmer occurring between '
                               '10 and 19 times' \
                               % (bin2, len(a),
                                  100 * float(bin2) / float(len(a)),
                                  merLengthVal))
          reportList.insert(5, '%d of %d / %0.4f%% of failures were triggered '
                               ' by the presence of a %dmer occurring between '
                               '20 and 49 times' \
                               % (bin3, len(a),
                                  100 * float(bin3) / float(len(a)),
                                  merLengthVal))
          reportList.insert(6, '%d of %d / %0.4f%% of failures were triggered '
                               'by the presence of a %dmer occurring between  '
                               '50 and 99 times' \
                               % (bin4, len(a),
                                  100 * float(bin4) / float(len(a)),
                                  merLengthVal))
          reportList.insert(7, '%d of %d / %0.4f%% of failures were triggered '
                               'by the presence of a %dmer occurring >100 '
                               'times' \
                               % (bin5, len(a),
                                  100 * float(bin5) / float(len(a)),
                                  merLengthVal))
          reportList.insert(8, '-' * 100)
      reportOut.write('\n'.join(reportList))
      reportOut.close()


def main():
    """Given a .bed file with probe sequences, runs Jellyfish to screen these
    for high abundance k-mers."""

    startTime = timeit.default_timer()

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
      '%s version %s. Requires a .bed file containing probe sequences in the '
      'fourth column. Also requires Jellyfish to be installed and in your '
      'PATH. Returns a .bed file in the same format as the input file '
      'containing only probes passing the specified kmer filter.' \
      % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The probe file to do kmer tallying in')
    requiredNamed.add_argument('-j', '--jf', action='store', required=True,
                               default='hg38_18mer.jf', type=str,
                               help='The Jellyfish .jf file to use')
    requiredNamed.add_argument('-m', '--merLength', action='store',
                               required=True, default=18, type=int,
                               help='The length of kmer used to build the .jf '
                                    'being used, default=18')
    requiredNamed.add_argument('-k', '--kmerThreshold', action='store',
                               required=True, default=5, type=int,
                               help='Filter probes with kmers occurring => '
                                    'than this many times, default=5')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str, help='The output name prefix')
    userInput.add_argument('-I', '--ID', action='store', type=str, default=None,
                           help='Specify an ID to be associated with temporary '
                                'file names & the default output name. Can be '
                                'useful for parallelization. Null by default. '
                                'Will not be automatically included in output '
                                'file name if -o is flagged')
    userInput.add_argument('-R', '--Report', action='store_true', default=False,
                          help='Write a Report file detailing the results of '
                               'the kmer filtering. Off by default. Note, '
                               'selecting this option will slow the script.')
    userInput.add_argument('-D', '--Debug', action='store_true', default=False,
                           help='The same as -Report, but prints info to '
                                'terminal instead of writing a log file. Off '
                                'by default')
    userInput.add_argument('-M', '--Meta', action='store_true', default=False,
                           help='Write a text file containing meta '
                                'information. Off by default. Reports input '
                                'file <tab> estimated runtime <tab> kmerFilter '
                                'version <tab> kmer occurrence threshold used '
                                '<tab> kmer length used <tab> number of input '
                                'probes passing kmer filter <tab> number of '
                                'input probes')

    # Import user-specified command line values.
    args = userInput.parse_args()
    inputFile = args.file
    outNameVal = args.output
    merLengthVal = args.merLength
    jellyfishFile = args.jf
    kVal = args.kmerThreshold
    IDval = args.ID
    reportVal = args.Report
    debugVal = args.Debug
    metaVal = args.Meta

    # Run the filter logic.
    runFilter(inputFile, outNameVal, merLengthVal, jellyfishFile, kVal, IDval,
              reportVal, debugVal, metaVal, startTime)

    # Print wall-clock runtime to terminal.
    print 'Program took %f seconds' % (timeit.default_timer() - startTime)

if __name__ == '__main__':
    main()
