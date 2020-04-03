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

