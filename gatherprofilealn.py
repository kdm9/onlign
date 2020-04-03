import argparse
from sys import stdin, stdout, stderr, exit
from pdgrep import getepiID, faread, writefa
import re

EPIRE = None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", type=argparse.FileType('w'), default=stdout,
            help="Selected sequences output file (fasta)")
    ap.add_argument("-I", "--id-regex", type=str, default=r'(EPI_ISL_\d+)',
            help="Regex to extract IDs. should have one match group")
    ap.add_argument("alignments",  type=str, nargs="+",
            help="Directory of `MAFFT --add-profile`-ed individual sequences, named EPI_ID.fasta")

    args = ap.parse_args()
    global EPIRE
    EPIRE = re.compile(args.id_regex)

    seqs = []
    # add the whole of the first alignment
    with open(args.alignments[0]) as fh:
        for seqid, seqdata in faread(fh):
            seqs.append(seqdata)
    # for each subsequent alignment, extract the respective sequence
    for aln in args.alignments[1:]:
        epiID = getepiID(aln)
        print(aln, epiID)
        with open(aln) as fh:
            fa = dict(faread(fh))
            seqs.append(fa[epiID])

    # check lengths match
    seqlen = None
    for seq in seqs:
        if seqlen == None:
            seqlen = len(seq["seq"])
        if len(seq["seq"]) != seqlen:
            print(f"WARNING: {seq['name']} isn't the same length as other sequences (typical {seqlen}, this one {len(seq['seq'])}", file=stderr)

    for seq in seqs:
        writefa(seq["name"], seq["seq"], file=args.output)

if __name__ == "__main__":
    main()
