import sys
from itertools import groupby
import pysam
from Bio import SeqIO

## SAM preparation, perhaps for a separate rule in Snakemake
## SAM file needs a header with @SQ, thus:
#samtools view -h -T reference.fasta -t reference.fasta.fai  L001.bam > L001.sam
# sot the SAM according to the QNAME
#samtools sort -O sam -n L001.sam > L001.sorted.sam

####################################################################################
## micro-API for getting back to cigar, sequence and quelity from the list
####################################################################################
def cigar_from_list(ll):
    cc = ""
    for l in ll:
        if (l[0] is None):
            cc = cc + "D"
        else:
            if (l[1] is None):
                cc = cc + "I"
            else:
                cc = cc + "M"
    tmp = ["".join(g) for k, g in groupby(cc)]
    #tmp = ["".join(g) for k, g in groupby(cc) if k != '-' and k != '/']
    cc = ""
    for l in tmp:
        cc = cc + str(len(l))
        cc = cc + l[0]
    return cc


def seq_from_list(ll):
    cc = ""
    for l in ll:
        if (l[0] is not None):   # not a deletion
            cc = cc + l[2]
    return cc

def qual_from_list(ll):
    cc = ""
    for l in ll:
        if (l[0] is not None):   # not a deletion
            cc = cc + l[3]
    return cc



####################################################################################
## "microAPI" to operate on the lists of tuples from get_aligned_pairs()
####################################################################################
def get_Max_APlist(ll):
    m = 0
    for l in ll:
        if l[1] is not None:
            if m < l[1]:
                m = l[1]
    return m

def get_Min_APlist(ll):
    m = 10e10
    for l in ll:
        if l[1] is not None:
            if m > l[1]:
                m = l[1]
    return m

# get first  (Nr) for a given second (pos), may be None
def get_at_Pos(ll, pos):
    out = -1
    for l in ll:
        if l[1] == pos:
            out = l[0]
    return(out)

def get_Tuple_at_Pos(ll, pos):
    out = -1
    for l in ll:
        if l[1] == pos:
            out = l
    return(out)

# get second  (pos) for a given first (nr), may be None
def get_at_Nr(ll, nr):
    out = -1
    for l in ll:
        if l[0] == nr:
            out = l[1]
    return(out)

def get_Tuple_at_Nr(ll, nr):
    out = -1
    for l in ll:
        if l[0] == nr:
            out = l
    return(out)

# get insert after a given pos
def get_Insert_After(ll, nr):
    out = []
    i = 0
    while ll[i][1] != nr and i < len(ll):
        #print(ll[i][1])
        i = i + 1
    if i < len(ll) :
        i = i + 1
    while i < len(ll) and ll[i][1] is None :
        out.append(ll[i])
        i = i + 1
    return(out)

# get a part until a given position
def get_To_Pos(ll, en):
    out = []
    i = 0
    while ll[i][1] is None or ll[i][1] < en:
        out.append(ll[i])
        i = i + 1
    return (out)

# get a part from a given position  to the end
def get_From_Pos(ll, st):
    out = []
    i = 0
    while ll[i][1] is None or ll[i][1] < st:
        i = i + 1
    out = ll[i:len(ll)]
    return (out)

####################################################################################
# utlity: query length from CIGAR
####################################################################################
def query_len(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    """
    read_consuming_ops = ("M", "I", "S", "=", "X")
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            result += length
    return result

####################################################################################
# main algorithm
####################################################################################
# this returns the full list for the fused read
def reconcile_Overlap(r1,r2):
    ll1 = get_aligned_pairs_extended(r1)
    ll2 = get_aligned_pairs_extended(r2)
    out = []
    en = get_Max_APlist(ll1)
    i = get_Min_APlist(ll2)
    p1 = get_To_Pos(ll1, get_Min_APlist(ll2) )
    p3 = get_From_Pos(ll2, get_Max_APlist(ll1) + 1)
    while i <= en:
        e1 = get_Tuple_at_Pos(ll1, i)
        e2 = get_Tuple_at_Pos(ll2, i)
        # Deletion always wins
        if e1[0] is not None and e2[0] is None:
            out.append(e2)
        else:
            out.append(e1)  # add code, which quality is higher
        i1 = get_Insert_After(ll1, i)
        i2 = get_Insert_After(ll2, i)
        # Longer insert wins
        if not (len(i1)==0 and len(i2)==0):
            if (len(i2) > len(i1)):
                out = out + i2
            else:
                out = out + i1
        i = i + 1
    out = p1 + out + p3
    return (out)


### extending the get_aligned_pairs tuples into 4-element ones, with nucleotide and quality
def get_aligned_pairs_extended(r1):
    llo = []
    ll1 = r1.get_aligned_pairs()
    ll = len(ll1)
    seq = r1.query_sequence
    qual = r1.qual
    c = 0
    for i in ll1:
        if i[0] is not None and i[1] is not None:
            llo.append( i + (seq[c], qual[c]) )
            c = c + 1
        if i[0] is None:   #DEL
            llo.append(i + (None, None))
        if i[1] is None:   #INS
            llo.append(i + (seq[c], qual[c]))
            c = c + 1
    return(llo)

def normalize_tuples_nr(ll):
    out = []
    c = 0
    for l in ll:
        if l[0] is None: # deletion
            out.append(l)
        else:
            out.append((c, l[1], l[2], l[3]))
            c = c + 1
    return out

def get_overlap(r1, r2):
    p2 = reconcile_Overlap(r1, r2)
    out = normalize_tuples_nr(p2)
    return out

####################################################################################
# fusion function on the level of a single alignment pair
####################################################################################

def read_fusion(r1, r2, refseq, header):
    outr = pysam.AlignedSegment(header=header)
    str1 = r1.get_aligned_pairs()  # with_seq = True MD tag not present
    cig1 = r1.cigartuples
    str2 = r2.get_aligned_pairs()
    cig2 = r2.cigartuples
    outr.qname = r1.qname
    outr.flag = 0  # not paired, super nice alignment ;) we ignore for now existing flags
    outr.reference_name = r1.reference_name
    outr.template_length = max(r1.template_length, r2.template_length)  # one is negative, (always?)
    # geometry
    if r1.pos > r2.pos:  # machniom! ie. swap
        tmp = r1
        r1 = r2
        r2 = tmp
    outr.pos = r1.pos
    # gap start and end
    gs = r1.pos + r1.query_alignment_length  # problem with soft clipping
    ge = r2.pos - 1
    print("gap start and end: " + str(gs) + " " + str(ge))
    if (ge - gs) <= 0:  # when no gap, but an overlap
        zz = get_overlap(r1, r2)
        outr.cigarstring = cigar_from_list(zz)
        outr.seq  = seq_from_list(zz)
        outr.qual  = qual_from_list(zz)
        restlen = r2.query_alignment_length - gs + ge
        print(len(outr.seq), len(r1.seq), len(r2.seq[(len(r2.seq) - restlen - 1):]), len(outr.seq))
    # in case of a gap
    else:
        outr.cigarstring = str(outr.template_length) + "M"  #obvious cheating
        padder = str(refseq[gs:(ge+1)].seq)
        outr.seq = r1.seq + padder + r2.seq
        #print(len(outr.seq), len(r1.seq), len(padder), len(r2.seq), r1.template_length, r2.template_length)
        padder = "F" * (ge - gs + 1)
        outr.qual = r1.qual + padder + r2.qual
    outr.mapping_quality = min(r1.mapping_quality, r2.mapping_quality)
    outr.next_reference_name = "*"
    outr.next_reference_start = 0
    outr.pnext = 0
    print(r1.qname + "  " + r2.qname)
    #print(outr.seq)
    return outr




def fuse_reads(argv):
    fname = argv[1]
    #fname = argv
    print("start processing "+fname)
    samfile = pysam.AlignmentFile(fname, "r")
    sam_out = pysam.Samfile("fused_output.sam", "w", header=samfile.header)
    sam_out_nonfused = pysam.Samfile("nonfused_output.sam", "w", header=samfile.header)

    prev = None
    i = 0
    reference = next(SeqIO.parse("reference.fasta", "fasta"))
    for read in samfile.fetch():
        if prev is not None:
            if prev.qname == read.qname:
                #if (read.pos > 0 and ( ('D' in read.cigarstring or 'D' in prev.cigarstring) and not 'S' in read.cigarstring and not 'S' in prev.cigarstring) ):
                if (read.pos > 0 and  ( not 'S' in read.cigarstring and not 'S' in prev.cigarstring) ): # for properly aligned, skipping unaligned, skipping softclipping
                    fused = read_fusion(prev, read, reference, samfile.header)
                    i = i + 1
                    cl = query_len(fused.cigarstring)
                    ql = fused.query_length
                    if (cl != ql):
                        print ("CIGAR length and query length not matching")
                    sam_out.write(fused)
                    sam_out_nonfused.write(prev)
                    sam_out_nonfused.write(read)
        prev = read
    print("Finished frusion "+str(i)+" pairs fused")

#### go:
#fuse_reads("Lara_mini.sam")


# scriptization
if __name__ == "__main__":
    fuse_reads(sys.argv)