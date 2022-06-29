#!/usr/bin/env python
import argparse
import sys

def revcomp(x:str) -> str:
    tr = str.maketrans('ATGC','TACG')
    return x.upper()[::-1].translate(tr)

def main():
    parser = argparse.ArgumentParser('merges/patches fastas based on paf alignment')
    parser.add_argument('--fastas',nargs = '+', help = 'fastas to merge')
    parser.add_argument('--paf',help = 'paf alignment file')
    parser.add_argument('--max_overhang', 
                        help = 'max bp a hit can overhang the end of the central sequence to be merged completely, default 1000',
                        default = 1000, type = int)
    args = parser.parse_args()
    
    seqdict = {}
    for fasta in args.fastas:
        print('processing fasta file ' + fasta, file = sys.stderr)
        workingseq = None
        seq = []
        with open(fasta) as f:
            for line in f:
                if line[0] == ">":
                    print('processing seq ' + line.split()[0][1:], file = sys.stderr)
                    if workingseq:
                        seqdict[workingseq] = "".join(seq)
                    workingseq = line.split()[0][1:]
                    seq = []
                else:
                    seq.append(line.strip())
        seqdict[workingseq] = "".join(seq)
    outseq = []
    matches = []
    with open(args.paf) as f:
        for line in f:
            fields = line.strip().split('\t')
            matches.append([fields[5],int(fields[7]),int(fields[8]),int(fields[6]),fields[0],int(fields[2]),int(fields[3]),fields[4]])
    matches.sort()
    centralseq = seqdict[matches[0][0]]
    ltrim = 0
    for match in matches:
        if match[1] < args.max_overhang:
            if match[7] == '+':
                outseq.append(seqdict[match[4]][:match[6]])
            elif match[7] == '-':
                outseq.append(revcomp(seqdict[match[4]][match[5]-1:]))
            ltrim = match[2]
        elif match[3] - match[2] < args.max_overhang:
            outseq.append(centralseq[ltrim:match[1] - 1])
            outseq.append(seqdict[match[4]][match[5] - 1:])
        else:
            outseq.append(centralseq[ltrim:match[1] - 1])
            ltrim = match[2]
            outseq.append(seqdict[match[4]][match[5] - 1:match[6]])
    print('>merged_' + matches[0][0] + '\n' + "".join(outseq))

if __name__ == "__main__":
    main()