#!/bin/bash
set -euo pipefail

# Default params
k=100 # number of most divergent sequences to align
datadir=./data

function usage() {
    echo "USAGE: $(basename $0) [OPTIONS] INPUT_ALIGNMENT"
    echo
    echo "OPTIONS:"
    echo "  -k INT    Use most divergent INT samples [default $k]."
    echo "  -d DIR    Working directory [default $datadir]."
}

if [ $# -eq 0 ]
then
    usage
    exit 0
fi

while getopts 'k:d:' OPTION; do
    case "$OPTION" in
        k)
            k=$OPTARG
            ;;
        d)
            datadir="$OPTARG"
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done
shift "$(($OPTIND -1))"

if [ $# -lt 1 ]
then
    usage
    exit 1
fi

input_fasta="$1"
input_base=$(basename $input_fasta .fasta)

iqtree="${IQTREE_PATH:-iqtree}"
mafft="${MAFFT_PATH:-mafft}"
lazy=yes


#######################################################################
#                          START MAIN SCRIPT                          #
#######################################################################

# how many sequences?
n=$(grep '>' $input_fasta | wc -l)

if [[ $lazy == no || $input_fasta -nt $datadir/${input_base}_guide.tree ]]
then
    if (( $n > 50000 )); then
        # use this method for maximum speed
        mafft --retree 0 --treeout --parttree --reorder $input_fasta > $datadir/${input_base}_sorted.fasta
        # replace all newlines
        # add branchlengths because this method doesn't have any.
        # assuming equal branchlengths is obviously bunk, but nevertheless
        # still gives a pragmatic way of selecting k divergent sequences
        tr -d '\n' < $input_fasta.tree | sed -e 's/),/)0.1,/g' -e  's/))/)0.1)/g'  > $datadir/${input_base}_guide.tree
    else
        # use this method if it's feasible
        mafft --retree 0 --treeout --reorder $input_fasta > $datadir/${input_base}_sorted.fasta
        # replace all newlines
        tr -d '\n' < $input_fasta.tree > $datadir/${input_base}_guide.tree
    fi

    # now continue for all methods 

    # nuke temp guide tree output from mafft
    rm -f $input_fasta.tree
    # add a semicolon to the tree so iq-tree can read it
    echo ";" >>  $datadir/${input_base}_guide.tree
else
    echo Skipping inital guide tree generation
fi


if [[ $lazy == no || $datadir/${input_base}_guide.tree -nt $datadir/${input_base}_A_k.fasta ]]
then

    # get the k most dissimilar sequences using iq-tree
    $iqtree -pre $datadir/${input_base}_kselect -te $datadir/${input_base}_guide.tree -k $k

    # use custom python script to extract sequences by EPIid
    mkdir -p $datadir/leftovers/ 
    python3 pdgrep.py --leftovers $datadir/leftovers/ --seqs $input_fasta --pda $datadir/${input_base}_kselect.pda  --output $datadir/${input_base}_kselect.fasta

    # align the k most dissimilar sequences in MAFFT
    mafft --thread -1 $datadir/${input_base}_kselect.fasta > $datadir/${input_base}_A_k_untrimmed.fasta

    # remove columns with >=95% gap -- fixes weird alignments
    trimal -keepheader -gt 0.95 -in $datadir/${input_base}_A_k_untrimmed.fasta -out data/${input_base}_A_k.fasta

else
    echo "Skipping core alignment matrix generation (A_k)"
fi

if [[ $lazy == no || $datadir/${input_base}_A_k.fasta -nt $datadir/${input_base}_A_g.fasta ]]
then

    mkdir -p $datadir/profilealn/ 
    parallel --bar mafft --thread 1 --keeplength --addprofile {} $datadir/${input_base}_A_k.fasta \>  $datadir/profilealn/{/} 2\> $datadir/profilealn/{/.}.log :::  $datadir/leftovers/*.fasta

    python3 gatherprofilealn.py -o $datadir/${input_base}_A_g.fasta $datadir/profilealn/*.fasta
else
    echo "Skipping full alignment matrix generation (A_g)"
fi

################
#  Nuke temps  #
################

rm -rf $datadir/${input_base}_kselect*  $datadir/profilealn  $datadir/leftovers $datadir/${input_base}_guide.tree

