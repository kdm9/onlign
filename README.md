# onlign

Online alignment prototypes for ANU improvements to AUGUR

### Install deps

`conda env create -f environment.yml && conda activate onlign`

### run GISAID ncov

```
mkdir data/
wget -O data/gisaid_cov2020_sequences.fasta  $GISAID_DATA_URL

# see `bash ./alignment.sh` for advanced options
bash alignment.sh data/gisaid_cov2020_sequences.fasta 
```


## TODOs

- [ ] A more robust way of detecting the N most diverse samples that doesn't pick long tips or otherwise strange sequences
    - By which I mean prefiltering the alignments somehow so that the guide tree doesn't include strange samples
- [ ] remove known-dodgy sites and samples from alignment
- [ ] smarter handling of alignment funkyness that maintain compatibility with the recognised coordinate space
	- Alignment funkyness e.g. regions gap-or-n-only columns due to funky samples
- [ ] Verify that the "core" alignment matrix doesn't change between new sequences before just concatenating the new seqs together (in `gatherprofilealn.py`)
- [ ] Integrate treebuilding logic *a la* Rob's state machine diagram
- [ ] run with bits of Sebastian's 100k seq simulation

- 
