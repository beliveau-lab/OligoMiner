# --------------------------------------------------------------------------
# OligoMiner
# bedChainer.py
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
import sys
import csv

def main():
    """Collapses overlapping probes."""

    # Retrieve input file name from stdin if not provided as a command line arg.
    if len(sys.argv) < 2:
        inputname = raw_input('Please provide input file name: ')
    else:
        inputname = sys.argv[1]
    input = open(inputname,'rb')
    reader = csv.reader(input,delimiter='\t')
    sorteddbname = inputname.split('.')[0]+'_sorted.bed'
    sorteddb = open(sorteddbname,'wb')
    sortedwriter = csv.writer(sorteddb,delimiter='\t',quoting=csv.QUOTE_NONE)

    # Sort file first based on starting coordinate, and also dedup.
    sortedseq = sorted(reader,key=lambda col:int(col[1]))
    previousseq = ''
    for seq in sortedseq:
        if seq !=previousseq:
            sortedwriter.writerow(seq)
            previousseq = seq
        else:
            pass
    input.close()
    sorteddb.close()

    # Open sorted db.
    input = open(sorteddbname,'rb')
    reader = csv.reader(input,delimiter='\t')
    outputname = inputname.split('.bed')[0]+'_chain.bed'
    output = open(outputname,'wb')
    writer = csv.writer(output,delimiter='\t',quoting=csv.QUOTE_NONE)
    log = open(inputname.split('.bed')[0]+'-log.txt','wb')
    writelog = csv.writer(log,delimiter='\t',quoting=csv.QUOTE_NONE)

    # Create dbs.
    probedb = []
    overlapdb = {}
    chaindb = []
    totallength = 0

    for row in reader:
        # Add to probe database.
        probedb.append(row)

        writelog.writerow(['Inspecting new chain: '+str(row)])
        
        # Check for overlap.
        x = 1
        # Initiate key of starting coordinate in overlap database, no overlap
        # so far.
        overlapdb[row[1]] = []

        # Iterate through probe db to update overlap database.
        while x < len(probedb)+1:
            # As long as the current item examined is not the one just entered.
            if row != probedb[-x]:
                if int(probedb[-x][2]) >= int(row[1]):
                  # Overlap: the endCoord of the item in probedb is greater
                  # than the startCoord of new item
                    overlapdb[row[1]].append(probedb[-x][1])
                    overlapdb[probedb[-x][1]].append(row[1])
                else:
                    # No overlap, can break out of loop since won't overlap
                    # with previous entries
                    break
            x += 1

        maxTm = 0.0
        maxSize = 1

        # Start creating probe chains, starting with including this element
        # in case it is the first/only element. It contains the current row
        # and Tm.
        chaindb.append([float(row[4]),row])

        newchain = [0.0]
        for x in xrange(0, len(chaindb)):
            if (row[1] != chaindb[x][-1][1]) \
               and (row[1] not in overlapdb[chaindb[x][-1][1]]):
                # If this new element does not overlap with the last element of
                # current chain in question.
                newchain = chaindb[x] + [row]
                newchain[0] = chaindb[x][0] + float(row[4]) # Sum up total Tm.
                chaindb[x] = newchain # Replace current chain with new one.

            if maxSize < len(newchain):
                maxSize = len(newchain) # Change maxsize if needed.

            if maxTm < newchain[0]:
                maxTm = newchain[0] # Change maxTm if needed.

            if maxSize < len(chaindb[x]):
                maxSize = len(chaindb[x]) # Change maxsize if needed.

            if maxTm < chaindb[x][0]:
                maxTm = chaindb[x][0] # Change maxTm if needed.


        # Mark chains to be culled that are less than the maxsize AND have most
        # recent element. Due to the functioning of the del command, we need to
        # generate an array of indices and remove later (from end to front).
        delindex= []
        copies =0
        for x in xrange(0, len(chaindb)):
            # Anything less than maxSize should be deleted. Also, if the Tm is
            # less than the maxTm and the cahin has at least two elements, it
            # should also be deleted.
            if (len(chaindb[x]) < maxSize) or ((chaindb[x][0] < maxTm) \
               and (len(chaindb[x])>2)):
                delindex.append(x)
        # Cull out the chains, backwards to preserve index values.
        for x in xrange(1,len(delindex)+1):
            writelog.writerow(['Culling following chain because not max size ' \
                               'or smaller Tm: %s' \
                               % str(chaindb[delindex[-x]])])
            del chaindb[delindex[-x]]

        # If the last element in the first chain does not overlap with any of
        # those elements that overlap with the first element of the this same
        # chain,.and they aren't the same element, we can start storing the
        # first chain
        nooverlap = False
        for element in overlapdb[chaindb[0][1][1]]:
            # Iterate through each element that overlaps with the first unit in
            # the first chain.
            try:
                if chaindb[0][-1][1] not in overlapdb[element]:
                    # If the last element (using its start coord) overlaps with
                    # these same elements.
                    nooverlap = True 
                    pass # Continue through remaining overlapdb
                else:
                    nooverlap = False
                    break
            except KeyError:
                writelog.writerow(['Overlap item not found.  Moving along...'])

        if nooverlap == True:
            writelog.writerow(['%s does not overlap with any of these ' \
                               'elements: %s' \
                               % (str(chaindb[0][-1][1]),
                                  str(overlapdb[chaindb[0][1][1]]))])
            writelog.writerow(['culling by saving first chain...'])
            chaindb = [chaindb[0]]
            nooverlap = False

        # Log the chains.
        writelog.writerow(['writing out chains'])
        for chain in chaindb:
            writelog.writerow(chain)

        if (len(chaindb) == 1) and (len(chaindb[0])> 2):
            # If the chaindb only has one chain, and that one chain is not the
            # row just entered (length >1), start saving the probes found in the
            # chain EXCEPT for the last added item
            for x in xrange(1,len(chaindb[0])-1):
                writer.writerow(chaindb[0][x])
                totallength +=1

            # Keep the last element to nucleate the next tree of probedb and
            # overlapdb.
            probedb = probedb[-1:]
            writelog.writerow(['probedb is set at'+str(probedb)])
            overlapdb = {row[1]:[]}
            chaindb = []
            chaindb.append([float(row[4]),row]) # Reinitiate chain with Tm.
            writelog.writerow(['chaindb now has '+str(chaindb)])

        writelog.writerow(['-----------------NEXT ITERATION------------------'])

    input.close()

    if len(chaindb)>=1:
        # If the chaindb still has items after iterating through the input, we
        # must add the largest chain, or last remaining chain.
        lastchain = chaindb[0]
        for x in xrange(1,len(lastchain)):
            writer.writerow(lastchain[x])

    output.close()

    return outputname

if __name__ == '__main__':
    main()
