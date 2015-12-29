import os
import sys
from pbcore.io.FastaIO import FastaReader as sfr
import bisect

rc = dict(zip("ACTGactgNn-", "TGACtgacNn-"))
def revcomp (seq):
    return "".join(map(lambda x: rc[x], seq)[::-1])


def mapSeqPair(seqtup, seqoverhang):

     seqname, seq1, seq2 = seqtup
     qfn = "%s_temp1.fasta" %(seqname)
     rfn = "%s_temp2.fasta" %(seqname)
     q_out = open(qfn, 'w')
     q_out.write(">seq1\n%s\n" %(seq1[-seqoverhang:]))
     q_out.close()
     r_out = open(rfn, 'w')
     r_out.write(">seq2\n%s\n" %(seq2[0:seqoverhang]))
     r_out.close()
     print >>sys.stderr, "nucmer running . . ."
     os.system("nucmer --maxmatch --nosimplify %s %s --prefix %s_temp 2> %s_temp.err" %(rfn, qfn, seqname, seqname))
     os.system("delta-filter -1 %s_temp.delta > %s_temp_1.delta 2> %s_temp.err" %(seqname, seqname, seqname)) 
     os.system("show-coords -l -r -L 30 %s_temp_1.delta > %s_temp_1.coords 2> %s_temp.err" %(seqname, seqname, seqname))
     f = open("%s_temp_1.coords" %(seqname))
     for i in range(5):
         f.readline()
     line = f.readline()
     ll = line.split()
     if len(ll) < 5:
         print >>sys.stderr, "Bad Line (likely no hits)"
         print >>sys.stderr, line
         print >>sys.stderr, ll
         print >>sys.stderr, "---------------"
         return False
     # get the start and end of the ref (s2, e2) and the query (s1, e1)
     s2, e2, s1, e1 = map(int, [ll[0], ll[1], ll[3], ll[4]])
     os.system("rm %s*temp*" %(seqname))
     # if the size of the overlap is greater than the start of the overlap we'll count it as valid
     if e2-s2 > s2:
        print >>sys.stderr, "Good line"
        print >>sys.stderr, s2,e2,s1,e1
        print >>sys.stderr, line
        print >>sys.stderr, ll
        return seq2[e2:], s2,e2,s1,e1
     else:
         print >>sys.stderr, "Bad Line"
         print >>sys.stderr, line
         print >>sys.stderr, ll
         print >>sys.stderr, "---------------"
         return False

def conToSeq (numstr, conNumDict, fastaDict):
    # get strand and seq
    num, strand = numstr[:-1], numstr[-1]
    con = conNumDict[num]
    seq = fastaDict[con]
    if strand == "+":
        return seq
    else:
        return revcomp(seq)

def stichseq (seq1, s1, e1, seq2, s2, e2):
    pass
    

def fastaToDict (contigFastaFile):
    fastaDict = {}
    counter = 0
    for entry in sfr(contigFastaFile):
        counter += 1
        fastaDict[entry.name] = entry.sequence
    return fastaDict

def qMapToDict (qMapFile):
    qMap = {}
    with open(qMapFile) as f:
        for l in f:
            if l[0] == "#":
                continue
            ll = l.strip().split()
            qMap[ll[0]] = ll[1]
    return qMap

def readCmapToDict (cmapfile):
    cmapdict = {}
    with open (cmapfile) as f:
        for l in f:
            if l[0] != "#":
                ll = l.strip().split()
                scaff_id = ll[0]
                pos = int(float(ll[5]))
                cmapdict.setdefault(scaff_id, []).append(pos)
    for scaff_id in cmapdict:
        cmapdict[scaff_id].sort()
    return cmapdict

def getCutSiteGapSeq (cmapPosList, start, end, nickseq="GMWSWKCN"):# positive strand only = GCTCTTC"):
    sind = bisect.bisect(cmapPosList, start)
    eind = bisect.bisect_left(cmapPosList, end)
    nickstart = cmapPosList[sind]
    if eind < sind:
        return "N"*(end-start)
    gapSeqList = []
    nickindices = range(sind, eind+1)
    prevpos = start
    for nickindex in nickindices:
        nickpos = cmapPosList[nickindex]
        # might need to adjust this by cutsize?
        gapSeqList.append("N"*(nickpos-prevpos))
        gapSeqList.append(nickseq)
        prevpos = nickpos
    gapSeqList.append("N"*(end-prevpos))
    return "".join(gapSeqList)
    
    
                    
def overlapsBetweenContigs (fn, qMap, fastaDict, fn_tag, cmapDict, exclusions = {}, scaff2contigs={}):
    fafout = open("%s_output_test.fa" %(fn_tag), 'w')
    contiguntouchfout = open("%s_output_test_contigs_untouch.txt" %(fn_tag), 'w')
    contigoverlapfout = open("%s_output_test_contigs_overlap.txt" %(fn_tag), 'w')
    contigfailedoverlapfout = open("%s_output_test_contigs_failedoverlap.txt" %(fn_tag), 'w')
    usedContigDict = {} # contig id and revised (trimmed) sequeunce
    with open(fn) as f:
        f.readline()
        prevSeq = ""
        prevAnchorId = ""
        observedscaffs = set()
        for l in f:
            #anchorId   qId1    start1  end1    qStart1 qEnd1   orientation1    qId2    start2  end2    qStart2 qEnd2   orientation2    overlapType
            ll = l.strip().split()
            anchorId, qId1, start1, end1, qStart1, qEnd1, orientation1, qId2, start2, end2, qStart2, qEnd2, orientation2, overlapType = ll
            start1, end1, qStart1, qEnd1, start2, end2, qStart2, qEnd2 = map(int, [start1, end1, qStart1, qEnd1, start2, end2, qStart2, qEnd2])
            
            observedscaffs.add(anchorId)

            fastaId1 = qMap[qId1]
            fastaId2 = qMap[qId2]
            seq1 = fastaDict[fastaId1]
            seq2 = fastaDict[fastaId2]
    
            if orientation1 == "-":
                seq1 = revcomp(seq1)
            if orientation2 == "-":
                seq2 = revcomp(seq2)
            seqname = "%s_%s" %(qId1, qId2)
                
            usedContigDict.setdefault(fastaId1, seq1) # use full sequences as default
            usedContigDict.setdefault(fastaId2, seq2) # use full sequences as default

            optmapaligngap  = start2 - end1
            # get the distance the sequence is overhanging into the alignment space
            if orientation1 == "+":
                overhang1 = len(seq1) - qEnd1
            else:
                overhang1 = qEnd1 
            if orientation2 == "+":
                overhang2 = qStart2
            else:
                overhang2 = len(seq2) - qStart2

            if anchorId != prevAnchorId:
                if not prevAnchorId == "":
                    print >>fafout, ">%s\n%s" %(prevAnchorId, prevSeq)
                prevSeq = seq1
                prevAnchorId = anchorId

            sequencegap = optmapaligngap - overhang1 - overhang2 
            print >>sys.stderr, "optical align gap: %i, overhang1: %i, overhang2: %i, sequencegap: %i" %(optmapaligngap, overhang1, overhang2, sequencegap)
            print >>sys.stderr, "seq1 len: %i, s1: %i, e1: %i:" %(len(seq1), qStart1, qEnd1)
            print >>sys.stderr, "seq2 len: %i, s2: %i, e2: %i:" %(len(seq2), qStart2, qEnd2)
            print >>sys.stderr, ll
            

            if "sequence-untouch" in overlapType:
                print >>sys.stderr, "sequence untouch - true gap is %i" %(sequencegap)
                gapSequence = getCutSiteGapSeq (cmapDict[anchorId], end1+overhang1, start2-overhang2)
                # in case the optical mapping results DO NOT suggest overlap we used the predicted gap
                # This is true EVEN If the contigs overlap in sequence space
                # (the optical maps should be a more reliable estimate in the case of TR expansion events)
                # -----------|------- sequence contig 1
                #           nick                      --------------|-------------------- sequence contig 2
                #                                                  nick
                # 
                # -----------|-------NNNNNNNNNNNNNNNNN---------------|-------------------- merged sequence
                #prevSeq = prevSeq + "N"*sequencegap + seq2
                prevSeq = prevSeq + gapSequence + seq2
                print >>contiguntouchfout, qMap[qId1]
                print >>contiguntouchfout, qMap[qId2]
                continue

            # only try to align with some tolerance of the sequencegap
            # but, at least try 20kb to give for inconsistency in alignment optical map prediction
            # otherwise require it to be within at least 3X the predicted size
            seqoverhang = max(3*sequencegap, 10000)

            seqoverlapi = mapSeqPair((seqname, seq1, seq2), seqoverhang)
            
            print >>sys.stderr,"*"*80

            if seqoverlapi:
                seq, s2, e2, s1, e1 = seqoverlapi
                # in case of successful merge
                # -----------|-----*--------*--- sequence contig 1
                #           nick  -*--------*---------------------------|-------------------- sequence contig 2
                #               (s1,s2)  (e1, e2)                     nick
                # -----------|-------------- trim
                #                          (e1) # note nucmer is 1-based
                #                           ----------------------------|---------------------
                #                          (e2) # note nucmer is 1-based
                # -----------|------------------------------------------|-------------------- merged sequence
                
                # we need to correct back the prevSeq since this already includes the end of seq1 
                # recall that we only include the seqoverhang length of the end of seq1 (and beginning of seq2)
                # for overlap mapping
                prevseqend = len(prevSeq)-(seqoverhang - e1) # end of current scaff seq
                prevSeq = prevSeq[0:prevseqend] + seq2[e2:]
                currusedseq1 = usedContigDict[fastaId1]
                # adjust sequences in used contig dict
                usedContigDict[fastaId1] = currusedseq1[0:len(currusedseq1)-(seqoverhang-e1)]
                usedContigDict[fastaId2] = usedContigDict[fastaId2][e2:]
                
                print "successful merge %s %i %i %i %s %i %i %i %s" %(qId1, s1,e1, len(seq1), qId2, s2,e2, len(seq2), l.strip())
                print >>contigoverlapfout, qMap[qId1]
                print >>contigoverlapfout, qMap[qId2]

            else:
                # print gap size for unsucessful merges
                
                print "unsuccessful merge with gap %s (seq1len: %i, seq2len: %i), %s -1 -1 %s -1 -1 %s" %(sequencegap, len(seq1), len(seq2), qId1, qId2, l.strip())
                # in case of unsuccessful merge this means either 
                #1.) the optical mapping prediction of overlap was incorrect
                #2.) There is a structural variation at the junction
                #3.) There was insufficient sequence overlapping to create a high-quality overlap alignment
                # all three cases are directly addressed by "chewing back" the region of overlap
                # -----------|-------------- sequence contig 1
                #           nick          ----------------------------------|-------------------- sequence contig 2
                #                       ?    ?                             nick
                # here we will be conservative and return
                # -----------|NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN|-------------------- merged sequence
                if optmapaligngap <= 0:
                    print >>sys.stderr, "Note these guys had overlapping mappings!"
                    sys.exit()
                    
                prevSeq = prevSeq[0:len(prevSeq)-overhang1] + "N"*(-1*sequencegap) + seq2[overhang2:]
                currusedseq1 = usedContigDict[fastaId1]
                # adjust sequences in used contig dict
                usedContigDict[fastaId1] = currusedseq1[0:len(currusedseq1)-overhang1]
                usedContigDict[fastaId2] = usedContigDict[fastaId2][overhang2:]
                print >>contigfailedoverlapfout, qMap[qId1]
                print >>contigfailedoverlapfout, qMap[qId2]
            print >>sys.stderr, "*"*80
            sys.stderr.flush()
            fafout.flush()
            contiguntouchfout.flush()
            contigoverlapfout.flush()
            contigfailedoverlapfout.flush()
            sys.stdout.flush()
        if not prevAnchorId == "":
            print >>fafout, ">%s\n%s" %(prevAnchorId, prevSeq)
    for scaff in scaff2contigs:
        if not scaff in observedscaffs:
            print >>sys.stderr, "scaff %s not seen in overlap analysis" %(scaff)
            contigs = scaff2contigs[scaff]
            if len(contigs) > 1:
                print >>sys.stderr, "ERROR:  more than 1 contig in this missed scaffold: %s; contigs = %s" %(scaff, ", ".join(contigs))
            elif len(contigs) == 1:
                usedContigDict.setdefault(qMap[contigs[0]], fastaDict[qMap[contigs[0]]])
                try:
                    print >>sys.stderr, "added singleton: %s" %(contigs[0])
                    contigseq = fastaDict[qMap[contigs[0]]]
                except:
                    print >>sys.stderr, "contig (%s) doesn't exist" %(contigs[0])
                print >>fafout, ">%s\n%s" %(scaff, contigseq)
            else:
                print >>sys.stderr, "WARNING:  more than 0 contigs in this missed scaffold: %s" %(scaff)
                


def agp_to_dict (agpfile):
    scaff2contigs = {}
    with open(agpfile) as f:
        for l in f:
            # make sure that its a scaffold line (not an empty value)
            # make sure that its not a header line
            if l[0] != "#" and "Scaffold" in l:
                ll = l.strip().split()
                try:
                    scaff, contig = ll[0], ll[5]
                except:
                    print >>sys.stderr, "couldn't find this contig:", ll
                scaff2contigs.setdefault(scaff, []).append(contig)
    return scaff2contigs

def xmap_to_dict (xmapfile, exclusions = {}):
    scaff2contigs = {}
    with open(xmapfile) as f:
        for l in f:
            if l[0] != "#":
                ll = l.strip().split()
                scaff, contig = ll[2], ll[1]
                if (scaff, contig) in exclusions:
                    print >>sys.stderr, "%s,%s excluded" %(scaff, contig)
                    continue
                scaff2contigs.setdefault(scaff, []).append(contig)
    return scaff2contigs

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-e", "--exclusion_xmap", dest="exclusion_xmap",
                  help="list of entries to exclude", metavar="FILE", default=None)
parser.add_option("-p", dest="fn_tag",
                  help="outputprefix", default = "merged_output")
parser.add_option("-a", "--anchor_file", dest="anchor_file",
                  help="name of anchor file (agp or xmap)", metavar="FILE", default = None)
parser.add_option("--produce_agp", 
                  help="name of output agp and output contig fasta (trimmed)", metavar="FILE", default = None)
(options, args) = parser.parse_args()
#contigFastaFile = sys.argv[1]
print options
print args
contigFastaFile = args[0]
fastaDict = {}
#qMapFile = sys.argv[2]
qMapFile = args[1]
qMap = qMapToDict(qMapFile)

#overlap_file = sys.argv[3]
overlap_file = args[2]
cmap_file = args[3] # this is the accompany hybrid cmap for the xmap
anchor_file = options.anchor_file # usually an xmap
fn_tag = options.fn_tag


# read cmap to dict
cmapDict = readCmapToDict (cmap_file)

# build exclusions
# note, this isn't really used anymore
exclusion_xmap = options.exclusion_xmap
if exclusion_xmap is not None:
    exclusionDict = xmap_to_dict (exclusion_xmap)
else:
    exclusionDict = {}
exclusions = []
for key,values in exclusionDict.items():
    for v in values:
        exclusions.append((key, v))

#if len(sys.argv) == 6:
if not(anchor_file is None):
    #anchor_file = sys.argv[4]
    #if anchor_file[-5:] == ".xmap":
    scaff2contigs = xmap_to_dict(anchor_file, exclusions)
    #print scaff2contigs
    #else:
    #    scaff2contigs = agp_to_dict(anchor_file)
        
    print >>sys.stderr, "starting read fasta to dict"
    fastaDict = fastaToDict (contigFastaFile)
    print >>sys.stderr, "finished reading fasta to dict"

    overlapsBetweenContigs (overlap_file, qMap, fastaDict, fn_tag, cmapDict,  scaff2contigs=scaff2contigs)
else:

    print >>sys.stderr, "starting read fasta to dict"
    fastaDict = fastaToDict (contigFastaFile)
    print >>sys.stderr, "finished reading fasta to dict"

    overlapsBetweenContigs (overlap_file, qMap, fastaDict, fn_tag, cmapDict)
#overlapsBetweenContigsThreaded (sys.argv[3], qMap, fastaDict)

