#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# outputClean.py
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
scriptName = 'outputClean'

# Specify script version.
Version = '1.7'

# Import module for handling input arguments.
import argparse

# Import regex library.
import re

# Import Biopython modules.
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC

# Import timeit module and record start time. This provides a rough estimate of
# the wall clock time it takes to run the script.
import timeit

# Define Tm calculation function.
def probeTm(seq1, sal, form):
    """Calculates the melting temperature of a given sequence under the
    specified salt and formamide conditions."""
    tmval = ('%0.2f' % mt.Tm_NN(seq1, Na=sal))
    fcorrected = ('%0.2f' % mt.chem_correction(float(tmval), fmd=form))
    return fcorrected


def cleanOutput(inputFile, uniqueVal, zeroVal, probVal, tempVal, sal, form,
                reportVal, debugVal, metaVal, outNameVal, startTime):
    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    # Open input file for reading.
    with open(inputFile, 'r') as f:
      file_read = [line.strip() for line in f]

    # Determine how many unique candidates are in the .sam file
    samIDs = [x.split('\t')[0].split(':')[1].split('-')[0] \
              if x[0] is not '@' else ' ' for x in file_read]
    candsSet = set()
    for x in samIDs:
      if x is not ' ':
          candsSet.add(x)

    # Make a list to hold the output.
    outList = []

    # Make lists to hold Report info if desired.
    if reportVal or debugVal is True:
      rejectList = []
      reportList = []

    if uniqueVal or zeroVal is True:
      # Process .sam file, keeping probes with only 0 or 1 unique alignment.
      for i in range(0, len(file_read), 1):
          if file_read[i][0] is not '@':
              chromField = file_read[i].split('\t')[2]
              chrom = file_read[i].split('\t')[0].split(':')[0]
              start = file_read[i].split('\t')[0].split(':')[1].split('-')[0]
              stop = file_read[i].split('\t')[0].split('-')[1].strip(' ')
              seq = file_read[i].split('\t')[9]
              Tm = probeTm(seq, sal, form)

              # For unique mode.
              if uniqueVal is True:
                  if re.match('\*', chromField) is None \
                     and re.search('XS', file_read[i]) is None:
                      outList.append('%s\t%s\t%s\t%s\t%s' \
                                     % (chrom, start, stop, seq, Tm))
                      # Report info on selected probe if desired.
                      if reportVal is True:
                          reportList.append('Candidate probe at %s:%s-%s '
                                            'aligned 1 time, added to output' \
                                            % (chrom, start, stop))
                      if debugVal is True:
                          print('Candidate probe at %s:%s-%s aligned 1 time, '
                                'added to output' % (chrom, start, stop))

                  else:
                      # Report info on rejected candidates if desired.
                      if reportVal or debugVal is True:
                          if start not in rejectList:
                              rejectList.append(start)
                              if re.match('\*', chromField) is not None:
                                  if reportVal is True:
                                      reportList.append('Candidate probe at '
                                                        '%s:%s-%s aligned 0 '
                                                        'times, was not added '
                                                        'to output' \
                                                        % (chrom, start, stop))
                                  if debugVal is True:
                                      print('Candidate probe at %s:%s-%s '
                                            'aligned 0 times, was not added to '
                                            'output' % (chrom, start, stop))
                              elif re.search('XS', file_read[i]) is not None:
                                  if reportVal is True:
                                      reportList.append('Candidate probe at '
                                                        '%s:%s-%s aligned >1 '
                                                        'time, was not added '
                                                        'to output' \
                                                        % (chrom, start, stop))
                                  if debugVal is True:
                                      print('Candidate probe at %s:%s-%s '
                                            'aligned >1 time, was not added to '
                                            'output' % (chrom, start, stop))

              # For zero mode.
              elif zeroVal is True:
                  if re.match('\*', chromField) is not None:
                      outList.append('%s\t%s\t%s\t%s\t%s' \
                                     % (chrom, start, stop, seq, Tm))
                      # Report info on selected probe if desired.
                      if reportVal is True:
                          reportList.append('Candidate probe at %s:%s-%s '
                                            'aligned 0 times, added to output '
                                            '(Zero mode active)' \
                                            % (chrom, start, stop))
                      if debugVal is True:
                          print('Candidate probe at %s:%s-%s aligned 0 times, '
                                'added to output (Zero mode active)' \
                                % (chrom, start, stop))
                  else:
                      # Report info on rejected candidates if desired.
                      if reportVal or debugVal is True:
                          if start not in rejectList:
                              rejectList.append(start)
                              if reportVal is True:
                                  reportList.append('Candidate probe at '
                                                    '%s:%s-%s aligned >0 '
                                                    'times, was not added to '
                                                    'output (Zero mode '
                                                    'active)' \
                                                    % (chrom, start, stop))
                              if debugVal is True:
                                  print('Candidate probe at %s:%s-%s aligned '
                                        '>0 times, was not added to output '
                                        '(Zero mode active)' \
                                        % (chrom, start, stop))

    # Else use LDA model.
    else:
      # Import scikit-learn LDA module.
      # Note the module name changed between sklearn versions 0.16 and 0.17
      from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

      # Import numpy module.
      import numpy as np

      # LDA model information.
      tempList = [32, 37, 42, 47, 52, 57]
      coefList = [[[-0.14494789, 0.18791679, 0.02588474]],
                  [[-0.13364364, 0.22510179, 0.05494031]],
                  [[-0.09006122, 0.25660706, 0.1078303]],
                  [[-0.01593182, 0.24498485, 0.15753649]],
                  [[0.01860365, 0.1750174, 0.17003374]],
                  [[0.03236755, 0.11624593, 0.24306498]]]
      interList = [-1.17545204, -5.40436344, -12.45549846,
                   -19.32670233, -20.11992898, -23.98652919]
      classList = [-1, 1]

      # Convert lists to ndarrays.
      coefArray = np.asarray(coefList)
      interArray = np.asarray(interList)
      classArray = np.asarray(classList)

      # Determine which index to reference for model values.
      np_index = tempList.index(tempVal)

      # Build model from encoded values.
      clf = LinearDiscriminantAnalysis()
      clf.coef_ = coefArray[np_index]
      clf.intercept_ = interArray[np_index]
      clf.classes_ = classArray

      # Determine which classifier parameters to use.
      clfT = tempList.index(tempVal)

      # Make lists to hold data about candidates.
      testList = []
      testSet = set()
      candsInfo = []

      # Process .sam file and extract information about each candidate probe.
      for i in range(0, len(file_read), 1):
          if file_read[i][0] is not '@':
              chromField = file_read[i].split('\t')[2]
              chrom = file_read[i].split('\t')[0].split(':')[0]
              start = file_read[i].split('\t')[0].split(':')[1].split('-')[0]
              stop = file_read[i].split('\t')[0].split('-')[1].strip(' ')
              seq = file_read[i].split('\t')[9]
              Tm = probeTm(seq, sal, form)

              # First look for candidate probes with only one unique alignment.
              if re.match('\*', chromField) is None \
                 and re.search('XS', file_read[i]) is None:
                  outList.append('%s\t%s\t%s\t%s\t%s' \
                                 % (chrom, start, stop, seq, Tm))
                  # Record info on selected probe if desired.
                  if reportVal is True:
                      reportList.append('Candidate probe at %s:%s-%s aligned '
                                        '1 time, added to output' \
                                        % (chrom, start, stop))
                  if debugVal is True:
                      print('Candidate probe at %s:%s-%s aligned 1 time, '
                            'added to output' % (chrom, start, stop))

              # Populate lists that will be used to make the classification
              # model input.
              else:
                  if re.match('\*', chromField) is None \
                     and start not in testSet:
                      t = [float(len(seq)),
                           float(file_read[i].split('\t')[12].split(':')[2]),
                           GC(seq)]
                      testList.append(t)
                      testSet.add(start)
                      candsInfo.append('%s\t%s\t%s\t%s\t%s' \
                                       % (chrom, start, stop, seq, Tm))
                  else:
                      # Report info on rejected candidates if desired.
                      if reportVal or debugVal is True:
                          if re.match('\*', chromField) is not None:
                              if start not in rejectList:
                                  rejectList.append(start)
                                  if reportVal is True:
                                      reportList.append('Candidate probe at '
                                                        '%s:%s-%s aligned 0 '
                                                        'times, was not added '
                                                        'to output' \
                                                        % (chrom, start, stop))
                                  if debugVal is True:
                                      print('Candidate probe at %s:%s-%s '
                                            'aligned 0 times, was not added to '
                                            'output' % (chrom, start, stop))

      # Make ndarray for input into classifier.
      testArray = np.asarray(testList)

      # Create classifier
      clf = LinearDiscriminantAnalysis()

      # Load temperature-specific model information.
      clf.coef_ = coefArray[clfT]
      clf.intercept_ = interArray[clfT]
      clf.classes_ = classArray

      # Use model to predict the probability that candidate
      # probes will have thermodynamically relevant
      # off-target binding sites unless all have just 1
      # alignment in the .sam file.
      if len(testArray) > 1:
          probs = clf.predict_proba(testArray)[:, 1]

          # Filter through tested candidates using
          # based on user-specified probability threshold.
          for i in range(0, len(probs), 1):
              if float(probs[i]) < probVal:
                  outList.append(candsInfo[i])
                  if reportVal is True:
                      reportList.append('Candidate probe at %s:%s-%s added to '
                                        'output with %0.4f < %0.4f probability of '
                                        'having off-target sites' \
                                        % (candsInfo[i].split('\t')[0],
                                           candsInfo[i].split('\t')[1],
                                           candsInfo[i].split('\t')[2],
                                           probs[i], probVal))
                  if debugVal is True:
                      print('Candidate probe at %s:%s-%s added to output with '
                            '%0.4f < %0.4f probability of having off-target sites'
                            % (candsInfo[i].split('\t')[0],
                               candsInfo[i].split('\t')[1],
                               candsInfo[i].split('\t')[2],
                               probs[i], probVal))
              else:
                  if reportVal is True:
                      reportList.append('Candidate probe at %s:%s-%s filtered with '
                                        '%0.4f => %0.4f probability of having '
                                        'off-target sites' \
                                        % (candsInfo[i].split('\t')[0],
                                           candsInfo[i].split('\t')[1],
                                           candsInfo[i].split('\t')[2],
                                           probs[i], probVal))
                  if debugVal is True:
                      print('Candidate probe at %s:%s-%s filtered with '
                            '%0.4f => %0.4f probability of having off-target sites'
                            % (candsInfo[i].split('\t')[0],
                               candsInfo[i].split('\t')[1],
                               candsInfo[i].split('\t')[2],
                               probs[i], probVal))
      # Sort output list.
      outList.sort(key=lambda x: [int(x.split('\t')[1])])

    # Determine the name of the output file.
    if outNameVal is None:
      outName = '%s_probes' % fileName
    else:
      outName = outNameVal

    # Create the output file.
    output = open('%s.bed' % outName, 'w')

    # Write the output file.
    output.write('\n'.join(outList))
    output.close()

    # Print info about the results to terminal.
    candsNum = len(candsSet)
    cleanNum = len(outList)
    if zeroVal is True:
      print('outputClean identified %d of %d / %0.4f%% candidate probes as '
            'having zero alignments' \
            % (cleanNum, candsNum, float(cleanNum) / float(candsNum) * 100))
    elif uniqueVal is True:
      print('outputClean identified %d of %d / %0.4f%% candidate probes as '
            'unique' % (cleanNum, candsNum,
                        float(cleanNum) / float(candsNum) * 100))
    else:
      print('outputClean passed %d of %d / %0.4f%% candidate probes through '
            'specificity filtering using the %dC LDA model' \
            % (cleanNum, candsNum,
               float(cleanNum) / float(candsNum) * 100, tempVal))

    # Write meta information to a .txt file if desired.
    if metaVal is True:
      metaText = open('%s_outputClean_meta.txt' % outName, 'w')
      metaText.write('%s\t%f\t%s\t%d\t%d' \
                     % (inputFile,
                        timeit.default_timer() - startTime,
                        Version, cleanNum, candsNum))
      metaText.close()

    # If desired, create report file.
    if reportVal is True:
      reportOut = open('%s_outputClean_log.txt' % outName, 'w')
      reportList.sort(key=lambda x: [int(x.split(':')[1].split('-')[0])])
      reportList.insert(0, 'Results produced by %s %s' % (scriptName, Version))
      reportList.insert(1, '-' * 100)
      if uniqueVal is True:
          reportList.insert(2, 'outputClean returned %d of %d / %0.4f%% '
                               'candidate probes as having exactly 1 ' 
                               'alignment' \
                               % (cleanNum, candsNum,
                                  float(cleanNum) / float(candsNum) * 100))
      elif zeroVal is True:
          reportList.insert(2, 'outputClean returned %d of %d / %0.4f%% '
                               'candidate probes as having 0 alignments (Zero '
                                'mode active)' \
                               % (cleanNum, candsNum,
                                  float(cleanNum) / float(candsNum) * 100))
      else:
          reportList.insert(2, 'outputClean passed %d of %d / %0.4f%% '
                               'candidate probes through specificity filtering '
                               'using the %dC LDA model' \
                               % (cleanNum, candsNum,
                                  float(cleanNum) / float(candsNum) * 100,
                                  tempVal))
      reportList.insert(3, '-' * 100)
      reportOut.write('\n'.join(reportList))
      reportOut.close()


def main():
    """Given a Sequence Alignment/Map (SAM) file, ouputs a Browser Extendable
    Data (BED) file containing probes that pass our temperature-specific Linear
    Discriminant Analysis (LDA) model."""

    startTime = timeit.default_timer()

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .sam file as input. Returns a .bed file '
        'containing only probes predicted to have one thermodynamically '
        'relevant target at the hybridization temperature provided by -T. '
        'Classification is performed using a temperature-specific linear '
        'discriminant analysis (LDA) model that requires scikit-learn 0.17+. '
        'Calculates the Tm of each probe based on -F and -s' \
        % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    mutEx = userInput.add_mutually_exclusive_group()
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .sam file to be processed')
    mutEx.add_argument('-l', '--lda', action='store_true', default=True,
                       help='Filter the SAM file using LDA model, On by '
                            'default.')
    mutEx.add_argument('-u', '--unique', action='store_true', default=False,
                       help='Only return probes aligning exactly one time. '
                             'Does not use the LDA model. Off by default.')
    mutEx.add_argument('-0', '--zero', action='store_true', default=False,
                       help='Only return probes aligning zero times. Does not '
                            'use the LDA model. Can be useful for targeting '
                            'transgenes/exogenous sequences. Off by default.')
    userInput.add_argument('-p', '--prob', action='store', default=0.5,
                           type=float,
                           help='The probability threshold for classifying a '
                                'candidate sequence as likely to have '
                                'off-target binding using the LDA model. '
                                'Default=0.5. Selecting smaller values will '
                                'improve precision (fewer false positives), but'
                                'at the expense of recall (more false '
                                'negatives). Selecting larger values will '
                                'improve recall at the expense of precision.')
    userInput.add_argument('-T', '--Temp', action='store', type=float,
                           default=42,
                           help='Specify the temperature-specific linear '
                                'discrimination model to use in LDM. Options '
                                'are 32, 37, 42, 47, 52, 57. Default=42')
    userInput.add_argument('-s', '--salt', action='store', default=390,
                           type=int,
                           help='The mM Na+ concentration to be used for Tm '
                                'calculation, default is 390')
    userInput.add_argument('-F', '--formamide', action='store', default=50,
                           type=float,
                           help='The percent formamide to be used for Tm '
                                'calculation, default is 50')
    userInput.add_argument('-R', '--Report', action='store_true', default=False,
                           help='Write a Report file detailing the results of '
                                '.sam cleaning. Off by default. Note, '
                                'selecting this option will slow the script.')
    userInput.add_argument('-D', '--Debug', action='store_true', default=False,
                           help='The same as -Report, but prints info to '
                                'terminal instead of writing a log file. Off '
                                'by default')
    userInput.add_argument('-M', '--Meta', action='store_true', default=False,
                           help='Write a text file containing meta '
                                'information. Off by default. Reports input '
                                'file <tab> estimated runtime <tab> '
                                'outputClean version <tab> unique probes '
                                'identified <tab> number of candidate probes '
                                'inputted')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the stem of the output filename')

    # Import user-specified command line values.
    args = userInput.parse_args()
    inputFile = args.file
    uniqueVal = args.unique
    zeroVal = args.zero
    probVal = args.prob
    tempVal = args.Temp
    sal = args.salt
    form = args.formamide
    reportVal = args.Report
    debugVal = args.Debug
    metaVal = args.Meta
    outNameVal = args.output

    cleanOutput(inputFile, uniqueVal, zeroVal, probVal, tempVal, sal, form,
                reportVal, debugVal, metaVal, outNameVal, startTime)

    # Print wall-clock runtime to terminal.
    print 'Program took %f seconds' % (timeit.default_timer() - startTime)

if __name__ == '__main__':
    main()
