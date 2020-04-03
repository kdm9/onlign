#!/usr/bin/env python3
import argparse
import re
from sys import stdout, stderr, exit
from os import path as op

EPIRE = re.compile(r'(EPI_ISL_\d+)')

def getepiID(idstr):
    m = EPIRE.search(idstr)
    if m is None:
        return None
    return m.group(1)

def faread(fafh):
    """Reads fasta, yields (epiID, data)"""
    seq = []
    name = None
    for line in fafh:
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                yield getepiID(name), {"seq": "".join(seq), "name": name}
            seq = []
            name = line[1:] # remove leading '>'
        else:
            seq.append(line)
    if name is not None:
        yield getepiID(name), {"seq": "".join(seq), "name": name}

def getpdseqids(pdfh):
    """Reads IQTREE pda file for list of most phylogenetically-diverse sequences, returning EXTRACTED ID"""
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
    """Writes fasta, with sequences lines no longer than `ll` to `file`."""
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
    ap.add_argument("-l", "--leftovers", type=str, default=None, metavar='DIR',
            help="Output each sequence NOT in pda file to its own fa file under DIR.")
    ap.add_argument("-p", "--pda", type=argparse.FileType('r'), required=True,
            help="IQ-TREE phylo diversity annotation file (made by iqtree -te $TREE -k $NCORE)")
    ap.add_argument("-I", "--id-regex", type=str, default=r'(EPI_ISL_\d+)',
            help="Regex to extract IDs. should have one match group")
    args = ap.parse_args()
    global EPIRE
    EPIRE = re.compile(args.id_regex)
    seqs = dict(faread(args.seqs))
    epiIDs = getpdseqids(args.pda)
    for epiID in epiIDs:
        seq = seqs[epiID]
        writefa(seq["name"], seq["seq"], file=args.output)
    print("Successfully extracted", len(epiIDs), "sequences", file=stderr)
    if args.leftovers:
        if not op.isdir(args.leftovers):
            os.makedirs(args.leftovers)
        epiIDs = set(epiIDs)
        n = 0
        for seqid, seq in seqs.items():
            if seqid not in epiIDs:
                n += 1
                with open(f"{args.leftovers}/{seqid}.fasta", "w") as ofh:
                    writefa(seq["name"], seq["seq"], file=ofh)
        print(f"Saved {n} non-diverse sequences to {args.leftovers}/", file=stderr)


if __name__ == "__main__":
    main()
