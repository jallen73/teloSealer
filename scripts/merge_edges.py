#!/usr/bin/env python
import argparse
import sys
import pandas

def revcomp(x:str) -> str:
    tr = str.maketrans('ATGC','TACG')
    return x.upper()[::-1].translate(tr)

def main():
    parser = argparse.ArgumentParser('merges/patches fastas based on paf alignment')
    parser.add_argument('--left', help = 'fasta from left telomere', default = None)
    parser.add_argument('--right', help = 'fasta from right telomere', default = None)
    parser.add_argument('--middle', help = 'fasta from main chromosome')
    parser.add_argument('--paf',help = 'paf alignment file')
    parser.add_argument('--max_overhang', 
                        help = 'max bp a hit can overhang the end of the central sequence to be merged completely, default 1000',
                        default = 1000, type = int)
    args = parser.parse_args()
    
    seqdict = {}
    seqnames = {'left':None,'right':None,'middle':None}
    for seqn,fasta in zip(['left','middle','right'],[args.left,args.middle,args.right]):
        if fasta:
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
                        seqnames[seqn] = workingseq
                        seq = []
                    else:
                        seq.append(line.strip())
            seqdict[workingseq] = "".join(seq)
    outseq = []
    matches = []
    matchdf = pandas.read_csv(args.paf,usecols=list(range(10)),header = None, sep = '\t')
    matchdf.columns = ['qid','qlen','qstart','qend','strand','sid','slen','sstart','send','matches','alen']
    lmatch = matchdf.query('qid == "' + seqnames['left'] + '"').sort_values('sstart')
    rmatch = matchdf.query('qid == "' + seqnames['right'] + '"').sort_values('send')
    lseq = seqdict[seqnames['left']][:list(lmatch['qstart'])[0]]
    mseq = seqdict[seqnames['middle']][list(lmatch['sstart'])[0]:list(rmatch['send'])[-1]]
    rseq = seqdict[seqnames['right']][list(rmatch['qend'])[-1]:]

    print('>merged_' + seqnames['middle'] + '\n' + lseq + mseq + rseq)

if __name__ == "__main__":
    main()