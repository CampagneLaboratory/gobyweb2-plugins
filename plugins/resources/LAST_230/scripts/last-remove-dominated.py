#! /usr/bin/env python

# Read MAF-format alignments, and write those are not "dominated" by
# any other one.  X dominates Y if they overlap on the top sequence,
# and the score of X in the overlapping region exceeds the total score
# of Y.  This script assumes the alignments have been sorted by
# maf-sort.sh.

import fileinput, string, re, itertools, optparse, signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # stop spurious error message

op = optparse.OptionParser(usage="%prog sorted-last-output.maf")
(opts, args) = op.parse_args()

def finalize_score_matrix(score_matrix, mask_lowercase):
    '''Add lowercase and non-standard letters to the score matrix.'''
    if ('a', 'a') in score_matrix: return  # it's already finalized
    min_score = min(score_matrix.values())
    if not mask_lowercase:
        for k, score in score_matrix.items():
            i, j = k
            score_matrix[i, j.lower()]         = score
            score_matrix[i.lower(), j]         = score
            score_matrix[i.lower(), j.lower()] = score
    for i in string.printable:
        for j in string.printable:
            k = i, j
            score_matrix.setdefault(k, min_score)

score_matrix = {}
gap_open = None
gap_extend = None
gap_pair = None  # extra parameter for "generalized affine gap costs"

def gap_cost(x, y):
    '''Get the cost for x and y unaligned letters in 2 sequences.'''
    affine      = gap_open * 2 + gap_extend * (x + y)
    generalized = gap_open + gap_extend * abs(x - y) + gap_pair * min(x, y)
    return min(affine, generalized)

def aln_score(seq1, seq2):
    '''Calculate the score of the alignment represented by seq1 and seq2.'''
    score = 0
    gap1 = gap2 = 0  # gap sizes
    for pair in itertools.izip(seq1, seq2):  # faster than zip
        a, b = pair
        if   a == '-': gap1 += 1
        elif b == '-': gap2 += 1
        else:
            if gap1 + gap2 > 0:
                score -= gap_cost(gap1, gap2)
                gap1 = gap2 = 0
            score += score_matrix[pair]
    if gap1 + gap2 > 0: score -= gap_cost(gap1, gap2)
    return score

def nthNonGap(s, n):
    '''Get the start position of the n-th non-gap.'''
    for i, x in enumerate(s):
        if x != '-':
            if n == 0: return i
            n -= 1
    raise ValueError('non-gap not found')

def nthLastNonGap(s, n):
    '''Get the end position of the n-th last non-gap.'''
    return len(s) - nthNonGap(s[::-1], n)

def aln_dominates(x, y):
    '''Does alignment x dominate alignment y?'''
    yscore = y[0][4]
    ybeg = y[1][4]
    yend = y[1][-1]
    xseq1 = x[1][12]
    xseq2 = x[2][12]
    xbeg = x[1][4]
    xend = x[1][-1]
    if len(xseq1) != len(xseq2): raise Exception("bad MAF data")
    lettersFromBeg = max(ybeg-xbeg, 0)
    lettersFromEnd = max(xend-yend, 0)
    alnBeg = nthNonGap(xseq1, lettersFromBeg)
    alnEnd = nthLastNonGap(xseq1, lettersFromEnd)
    xscore = aln_score(xseq1[alnBeg:alnEnd], xseq2[alnBeg:alnEnd])
    return xscore > yscore

def parse_maf_a_line(line):
    '''Parse a MAF line starting with "a".'''
    words = re.split(r'([\s=]+)', line)  # keep the delimiters
    words[4] = int(words[4])             # alignment score
    return words

def parse_maf_s_line(line):
    '''Parse a MAF line starting with "s".'''
    words = re.split(r'(\s+)', line)   # keep the spaces
    words[4] = int(words[4])           # aln start
    words[6] = int(words[6])           # aln size
    words.append(words[4] + words[6])  # aln end
    return words

def write_aln(aln):
    if aln.pop() == True: return  # the alignment is dominated
    aln[0] = ''.join(map(str, aln[0]))
    aln[1] = ''.join(map(str, aln[1][:-1]))  # remove the end that we appended
    aln[2] = ''.join(map(str, aln[2][:-1]))  # remove the end that we appended
    print ''.join(aln)  # this prints a blank line at the end

def process_aln(alns, newaln):
    if not newaln: return
    newaln[0] = parse_maf_a_line(newaln[0])
    newaln[1] = parse_maf_s_line(newaln[1])
    newaln[2] = parse_maf_s_line(newaln[2])
    newaln.append(False)  # means: this alignment is not (yet) dominated
    newchr = newaln[1][2]
    newbeg = newaln[1][4]
    kept_alns = []
    for oldaln in alns:
        oldchr = oldaln[1][2]
        oldbeg = oldaln[1][4]
        oldend = oldaln[1][-1]
        if oldchr == newchr and oldend > newbeg:
            if oldbeg > newbeg: raise Exception("MAF data isn't sorted")
            if not oldaln[-1]:  # speed-up
                if aln_dominates(newaln, oldaln): oldaln[-1] = True
            if not newaln[-1]:  # speed-up
                if aln_dominates(oldaln, newaln): newaln[-1] = True
            kept_alns.append(oldaln)
        else:
            write_aln(oldaln)
    alns[:] = kept_alns + [ newaln ]

alns = []  # alignments that might overlap the next one
aln = []  # the current alignment
columns = []  # column headings for the score matrix
mask_lowercase = False

for line in fileinput.input(args):
    if line.startswith('#'):
        body = line[1:].strip()
        words = body.split()
        m1 = re.search(r'a=(\d+) b=(\d+) c=(\d+)', body)
        m2 = re.search(r'u=(\d+)', body)
        m3 = re.match(r'\S(\s+\S)*$', body)
        m4 = re.match(r'\S(\s+-?\d+)+$', body)
        if m1:
            gap_open, gap_extend, gap_pair = map(int, m1.groups())
        elif m2:
            mask_lowercase = (m2.group(1) == "3")
        elif m3 and not columns:
            columns = words
        elif m4 and len(words) == len(columns) + 1:
            row = words[0]  # row heading for the score matrix
            scores = map(int, words[1:])
            for col, s in zip(columns, scores):
                score_matrix[row, col] = s
        elif columns:
            finalize_score_matrix(score_matrix, mask_lowercase)
    elif line.isspace():
        process_aln(alns, aln)
        aln = []
    else:
        aln.append(line)

process_aln(alns, aln)  # don't forget the last alignment

for oldaln in alns:
    write_aln(oldaln)
