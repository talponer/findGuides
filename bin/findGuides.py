#!/home/rdreos/miniconda2/bin/python

from __future__ import print_function
import sys
import re
import argparse
from Bio import SeqIO
from Bio import Entrez

######################################################################
### There is still one major issue to solve: if a giode has a CC in
### the active site, those CC are the seeding pont for a guide on the
### reverse strand.  These two guides would heavily interfere with
### each other (they woulg target the seeding bases of each other). It
### all depend if the guides are used independently or as a mixture
######################################################################

# parse the command line arguments:
parser = argparse.ArgumentParser(description='Given a RsfSeq Nucleotide ID, '
                                 'return all possible guides for the Cas9 '
                                 'protein. STDOUT gives the list of sequences, '
                                 'position in the CDS and a score; STDERR '
                                 'gives a log of the analysis done on each '
                                 'sequence.')
parser.add_argument('-r','--region',
                    help='Transcript region (nucleotides) with high '
                    'conservation. Example: 1:1000. Default: entire transcript',
                    type=str)
parser.add_argument('-m','--mut',
                    help='Position of known mutations separated by a comma. '
                    'Example: \'P134S,Y103H,Q68E\'',
                    type=str)
parser.add_argument('-g','--genome',
                    help='Single FASTA file containing the genome sequence. '
                    'It is used to check for guide multiple matches',
                    type=str,
                    default='/home/rdreos/Projects/annotation/human/'
                    'Homo_sapiens.GRCh38.89.dna.primary_assembly.fa')
parser.add_argument('-t','--tableLog',
                    help='Write log in tabular format. '
                    'Useful for downstream guides\' analysis and selection',
                    action="store_true")
parser.add_argument('-a','--allGuides',
                    help='Consider all guides, even if overlapping',
                    action="store_true")
parser.add_argument('-i','--id',
                    help='RefSeq Nucleotide ID. Example: \'NM_005334\'',
                    nargs=1,
                    type=str,
                    required=True)
args = parser.parse_args()

#if args.id:
#   print(args.id)

# define the genetic code:
genCode = {
    "TTT":"F",
    "TCT":"S",
    "TAT":"Y",
    "TGT":"C",
    "TTC":"F",
    "TCC":"S",
    "TAC":"Y",
    "TGC":"C",   
    "TTA":"L",
    "TCA":"S",
    "TAA":"Ter",
    "TGA":"Ter",
    "TTG":"L",
    "TCG":"S",
    "TAG":"Ter",
    "TGG":"W",
    "CTT":"L",
    "CCT":"P",
    "CAT":"H",
    "CGT":"R",
    "CTC":"L",
    "CCC":"P",
    "CAC":"H",
    "CGC":"R",
    "CTA":"L",
    "CCA":"P",
    "CAA":"Q",
    "CGA":"R",
    "CTG":"L",
    "CCG":"P",
    "CAG":"Q",
    "CGG":"R", 
    "ATT":"I",
    "ACT":"T",
    "AAT":"N",
    "AGT":"S",
    "ATC":"I",
    "ACC":"T",
    "AAC":"N",
    "AGC":"S",
    "ATA":"I",
    "ACA":"T",
    "AAA":"K",
    "AGA":"R",
    "ATG":"M",
    "ACG":"T",
    "AAG":"K",
    "AGG":"R", 
    "GTT":"V",
    "GCT":"A",
    "GAT":"D",
    "GGT":"G",
    "GTC":"V",
    "GCC":"A",
    "GAC":"D",
    "GGC":"G",
    "GTA":"V",
    "GCA":"A",
    "GAA":"E",
    "GGA":"G",
    "GTG":"V",
    "GCG":"A",
    "GAG":"E",
    "GGG":"G" 
}

def FileCheck(fn):
    try:
        open(fn, "r")
        return 1
    except IOError:
        print("Error: Genome file does not appear to exist.", file=sys.stderr)
        sys.exit()
        return 0

def eprint(*args, **end):
    if end:
        print(*args, file=sys.stderr, end=end)
    else:
        print(*args, file=sys.stderr)

def search_fasta(pattern):
#    eprint("Searching pattern: ", pattern)
    m = 0
    for chrom, seq in genome.iteritems():
#        eprint("  ", chrom)
        m += len(re.findall(pattern, seq))
    return m

def permuteBase(triplet, F, T):
    tripletM = []
    for c1 in re.finditer(F, triplet):
        c1 = c1.start()
        t1 = triplet[0:c1] + T + triplet[c1+1:3]
        tripletM.append(t1)
        if len(re.findall(F, t1)) > 0:
            for c2 in re.finditer(F, t1):
                c2 = c2.start()
                t2 = t1[0:c2] + T + t1[c2+1:3]
                tripletM.append(t2)
            if len(re.findall(F, t2)) > 0:
                for c3 in re.finditer(F, t2):
                    c3 = c3.start()
                    t3 = t2[0:c3] + T + t2[c3+1:3]
                    tripletM.append(t3)
    return list(set(tripletM))


Entrez.email = "rene.dreos@unil.ch"

refSeqId = str(args.id[0])
eprint('## ID:',refSeqId)
#refSeqId = "NM_005334"

#
#eprint('## logs:', args.tableLog)

# get the mutations
mutantStart = []
if args.mut:
    for m1 in args.mut.split(","):
        mutantStart.append((int(re.sub("[A-Z]","",m1))-1)*3)
    mutations = args.mut.split(",")
    eprint('## Mutations:',mutantStart)
else:
    mutantStart.append(-100)
    mutations = ['A1000000000A']
    
# read the genome sequence
if args.genome:
    file_path = args.genome

FileCheck(file_path)

print("## Reading genome...", file=sys.stderr, end='')
genome = {}
for record in SeqIO.parse(open(file_path, "rU"), "fasta"):
    chrom = record.id
    sequence = str(record.seq)
    genome[chrom] = sequence
print("Done", file=sys.stderr)

# specify which entry to download
print("## Fetching sequence...", file=sys.stderr, end='')
handle = Entrez.efetch(db="nucleotide", id=refSeqId, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
print("Done", file=sys.stderr)


print("## Parsing sequence...", file=sys.stderr, end='')
exonStart = []
for I in record.features:
    if (I.type == 'CDS'):
        cdsSeq = record.seq[I.location.start:I.location.end]
        cdsStart = I.location.start
        cdsEnd = I.location.end
    if (I.type == 'exon'):
        exonStart.append(I.location.start)
        #print ("%i %i" % (I.location.start, I.location.end))
print("Done", file=sys.stderr)

cdslength = cdsEnd - cdsStart
print("## CDS lenght (bp):", cdslength, file=sys.stderr)

# High conservation region, inferred with blastp
if args.region:
    highConsStart,highConsEnd = args.region.split(":")
else:
    highConsStart = 1
    highConsEnd = cdslength
eprint('## Range:',highConsStart,highConsEnd)

# Check which guides to print
if args.allGuides:
    print('## Output: all guides', file=sys.stderr)
else:
    print('## Output: non-overlapping guides (first one)', file=sys.stderr)

# print table header if requested
if args.tableLog:
    print('', file=sys.stderr)
    print('Seaquence','Strand','From','To','FOR-hit','REV-hit','Substitutions','Introns','NearMut','HighCons','ActiveC','Score','Print', sep="\t", file=sys.stderr)

######################################################################
### Start with the forward strand
######################################################################
oldStart = 0;
printed = 0
overlap = 0
lastStart = 0
for m in re.finditer('GG', str(cdsSeq)): # find the seed
    newStart = m.start()
    diff = newStart - oldStart

    if args.allGuides:
        diff = newStart
    
    # take only guides if they do not overlap, or check if first hit
    # is far enough
    if (diff > 20):

        toPrint = 0    
        beginning = newStart-20
        end = newStart+2
        guideSeq = str(cdsSeq[beginning:end])
        if args.tableLog:
            rend="\t"
            print(guideSeq,"+",str(beginning),str(end), sep="\t", end=rend, file=sys.stderr)
        else:
            rend="\n"
            print('', file=sys.stderr)
            print('Sequence:         ', guideSeq, file=sys.stderr)
            print('Transcript strand: Forward', file=sys.stderr)
            print('Position in CDS:  ', str(beginning)+'..'+str(end), file=sys.stderr)
        guideSeqRC = str(cdsSeq[beginning:end].reverse_complement())
        mGenome = search_fasta(guideSeq) # find genome matches
        mrGenome = search_fasta(guideSeqRC)
        #mGenome = 1
        #mrGenome = 0
        totMatches = mGenome + mrGenome
        if args.tableLog:
            print(mGenome,mrGenome, sep="\t", end="\t", file=sys.stderr)
        else:
            print('Genome match FOR: ', mGenome, file=sys.stderr)
            print('Genome match REV: ', mrGenome, file=sys.stderr)
        cMatch = len(re.findall('C', str(cdsSeq[beginning:beginning+5])))

        # Get the residues in the 5 last bases
        if not args.tableLog:
            print('Possible substit: ', file=sys.stderr, end=' ')
        newBeg = beginning
        printMut = []
        okMut = 0
        while True:
            modBed = newBeg % 3
            aapos = (newBeg+2)/3
            if modBed == 0:
                frame = 0
            elif modBed < 0.5:
                frame = 2
            else:
                frame = 1
            triplet = str(cdsSeq[newBeg-frame:newBeg-frame+3])
            tripletM = permuteBase(triplet,'C','T')
            amino = genCode[triplet]
            for mut in tripletM:
                aminoM = genCode[mut]
                m2 = amino + str(aapos) + aminoM
                if aminoM != amino:
                    toPrint += 1
                    printMut.append(m2)
                    okMut = 1
                else:
                    print('['+m2+']', file=sys.stderr, end=' ')
                    #toPrint -= 1
            newBeg = newBeg+3
            if newBeg > beginning+6:
                break
            
        for m1 in list(set(printMut)):
            for m2 in mutations:
                if m1 == m2:
                    toPrint += 10
                    print('*', file=sys.stderr, end='')
            print(m1, file=sys.stderr, end=' ')
        
        if args.tableLog:
            print('', end="\t", file=sys.stderr)
        else:
            print('', file=sys.stderr)
            print('Introns:          ', file=sys.stderr, end=' ')
        inIntr = 0
        for k in exonStart:
            k -= cdsStart
            if beginning < k: # if it spans an intron
                if end > k:
                    toPrint = -1000
                    inIntr = 1
        print(inIntr, file=sys.stderr, end=rend)
                    
        if not args.tableLog:
            print('Near known mut:   ', file=sys.stderr, end=' ')
        nearMut = 0
        for i in mutantStart:
            if abs(beginning-i) < 100: # if it is near the mutations
                toPrint += 1
                nearMut = 1
        print(nearMut, file=sys.stderr, end=rend)

        if not args.tableLog:
            print('High cons:        ', file=sys.stderr, end=' ')
        inCons = 0
        if beginning > int(highConsStart):
            if end < int(highConsEnd):
                toPrint += 1
                inCons = 1
            else:
                toPrint -= 1
        else:
            toPrint -= 1
        print(inCons, file=sys.stderr, end=rend)

        if not args.tableLog:
            print('Active Cs:        ', file=sys.stderr, end=' ')
        if cMatch >= 1: # if there are C in the most active region (-20:-16)
            toPrint += 1
            print('1', file=sys.stderr, end=rend)
        else:
            toPrint -= 10
            print('0', file=sys.stderr, end=rend)

        if args.tableLog:
            print(toPrint, file=sys.stderr, end=rend)
        else:
            print('Score:            ', toPrint, file=sys.stderr)

        if not args.tableLog:
            print('Printed:          ', file=sys.stderr, end=' ')        
        if totMatches == 1:
            if okMut == 1:
                if inIntr == 0:
                    print(guideSeq, str(beginning), str(end), toPrint, sep="\t")
                    #print(guideSeq, str(beginning), str(end), sep="\t", end="\t")
                    print('1', file=sys.stderr)
                    oldStart = end
                    printed = 1
                else:
                    print('0', file=sys.stderr)
            else:
                print('0', file=sys.stderr)
        else:
            print('0', file=sys.stderr)

######################################################################
### going through the reverse strand
######################################################################
printed = 0
oldStart = 0
overlap = 0
revMatch = []
for m in re.finditer('CC', str(cdsSeq)): # find the seed
    revMatch.append(m.start())

for index, newStart in enumerate(revMatch): # find the seed
    toPrint = 0
    if index+1 < len(revMatch):
        nextStart = revMatch[index+1]
    else:
        nextStart = len(str(cdsSeq))
    diff = nextStart - newStart

    if args.allGuides:
        diff = nextStart
    
    # take only guides if they do not overlap, or check if first hit
    # is far enough
    if (diff > 20 and newStart+22 < cdslength):

        toPrint = 0    
        end = newStart+22
        beginning = newStart
        guideSeq = str(cdsSeq[beginning:end])
        guideSeqRC = str(cdsSeq[beginning:end].reverse_complement())
        if args.tableLog:
            print(guideSeqRC,"-",str(end),str(beginning), sep="\t", end=rend, file=sys.stderr)
        else:
            print('', file=sys.stderr)
            print('Sequence:         ', guideSeqRC, file=sys.stderr)
            print('Transcript strand: REV', file=sys.stderr)
            print('Position in CDS:  ', str(end)+'..'+str(beginning), file=sys.stderr)
        mGenome = search_fasta(guideSeq) # find genome matches
        mrGenome = search_fasta(guideSeqRC)
        #mGenome = 1
        #mrGenome = 0
        totMatches = mGenome + mrGenome
        if args.tableLog:
            print(mGenome, mrGenome, sep="\t", end=rend, file=sys.stderr)
        else:
            print('Genome match FOR: ', mGenome, file=sys.stderr)
            print('Genome match REV: ', mrGenome, file=sys.stderr)
        cMatch = len(re.findall('C', str(cdsSeq[beginning:beginning+5])))

        # Get the residues in the 5 last bases
        if not args.tableLog:
            print('Possible substit: ', file=sys.stderr, end=' ')
        newBeg = end-6
        printMut = []
        okMut = 0
        while True:
            modBed = newBeg % 3
            aapos = (newBeg+2)/3
            if modBed == 0:
                frame = 0
            elif modBed < 0.5:
                frame = 2
            else:
                frame = 1
            triplet = str(cdsSeq[newBeg-frame:newBeg-frame+3])
            tripletM = permuteBase(triplet,'G','A')
            amino = genCode[triplet]
            for mut in tripletM:
                aminoM = genCode[mut]
                m2 = amino + str(aapos) + aminoM
                #eprint(triplet, amino, mut, aminoM)
                if aminoM != amino:
                    toPrint += 1
                    okMut = 1
                    printMut.append(m2)
                else:
                    print('['+m2+']', file=sys.stderr, end=' ')
            newBeg = newBeg+3
            if newBeg >= end:
                break
            
        for m1 in list(set(printMut)):
            for m2 in mutations:
                if m1 == m2:
                    toPrint += 10
                    print('*', file=sys.stderr, end='')
            print(m1, file=sys.stderr, end=' ')

        if args.tableLog:
            print('', end="\t", file=sys.stderr)
        else:
            print('', file=sys.stderr)
            print('Introns:          ', file=sys.stderr, end=' ')
        intr = 0
        for k in exonStart:
            k -= cdsStart
            if beginning < k: # if it spans an intron
                if end > k:
                    toPrint = -1000
                    intr = 1

        print(intr, file=sys.stderr, end=rend)

        if not args.tableLog:
            print('Near known mut:   ', file=sys.stderr, end=' ')
        nearMut = 0
        for i in mutantStart:
            if abs(beginning-i) < 100: # if it is near the mutations
                toPrint += 1
                nearMut = 1
        print(nearMut, file=sys.stderr, end=rend)
                
        if not args.tableLog:
            print('High cons:        ', file=sys.stderr, end=' ')
        inHigh = 0
        if beginning > int(highConsStart):
            #eprint('first')
            if end < int(highConsEnd):
                toPrint += 1
                inHigh = 1
            else:
                toPrint -= 1
        else:
            toPrint -= 1
        print(inHigh, file=sys.stderr, end=rend)

        if not args.tableLog:
            print('Active Cs:        ', file=sys.stderr, end=' ')
        if cMatch >= 1: # if there are C in the most active region (-20:-16)
            toPrint += 1
            print('1', file=sys.stderr, end=rend)
        else:
            toPrint -= 10
            print('0', file=sys.stderr, end=rend)

        if args.tableLog:
            print(toPrint, file=sys.stderr, end=rend)
        else:
            print('Score:            ', toPrint, file=sys.stderr)
            print('Printed:          ', file=sys.stderr, end=rend)
        if totMatches == 1:
            if okMut == 1:
                if intr == 0:
                    print(guideSeqRC, str(end), str(beginning), toPrint, sep="\t")
                    #print(guideSeqRC, str(end), str(beginning), sep="\t", end="\t")
                    print('1', file=sys.stderr)
                    oldStart = end
                    printed = 1
                else:
                    print('0', file=sys.stderr)
            else:
                print('0', file=sys.stderr)
        else:
            print('0', file=sys.stderr)
