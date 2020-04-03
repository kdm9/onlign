#!/usr/bin/env python3
import argparse
import re
from sys import stdout, stderr, exit

EPIRE = None

def getepiID(idstr):
    m = EPIRE.search(idstr)
    if m is None:
        return None
    return m.group(1)

def fadict(fafh):
    d = {}
    seq = []
    name = None
    for line in fafh:
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                d[getepiID(name)] = {"seq": "".join(seq), "name": name}
            seq = []
            name = line[1:] # remove leading '>'
        else:
            seq.append(line)
    if name is not None:
        d[getepiID(name)] = {"seq": "".join(seq), "name": name}
        seq = []
        name = line[1:] # remove leading '>'
    return d

def getpdseqids(pdfh):
    pdids = []
    # skip until list header line
    while not pdfh.readline().startswith('The optimal PD set has'):
        pass
    for line in pdfh:
        line = line.rstrip()
        if line.startswith("Corresponding sub-tree:"):
            break
        if line:
            pdids.append(getepiID(line))
    return pdids


def writefa(name, seq, file=None, ll=80):
    print(f">{name}", file=file)
    if ll is None or ll < 1:
        ll = len(seq)
    for i in range(0, len(seq), ll):
        end=min(len(seq), i+ll)
        print(seq[i:end], file=file)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", type=argparse.FileType('w'), default=stdout,
            help="Selected sequences output file (fasta)")
    ap.add_argument("-s", "--seqs", type=argparse.FileType('r'), required=True,
            help="Sequences, as fasta")
    ap.add_argument("-p", "--pda", type=argparse.FileType('r'), required=True,
            help="IQ-TREE phylo diversity annotation file (made by iqtree -te $TREE -k $NCORE)")
    ap.add_argument("-I", "--id-regex", type=str, default=r'(EPI_ISL_\d+)',
            help="Regex to extract IDs. should have one match group")
    args = ap.parse_args()
    global EPIRE
    EPIRE = re.compile(args.id_regex)
    seqs = fadict(args.seqs)
    epiIDs = getpdseqids(args.pda)
    for epiID in epiIDs:
        seq = seqs[epiID]
        writefa(seq["name"], seq["seq"], file=args.output)
    print("Successfully extracted", len(epiIDs), "sequences", file=stderr)


if __name__ == "__main__":
    main()
