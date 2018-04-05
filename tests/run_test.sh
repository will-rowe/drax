#!/usr/bin/env bash
git clone https://github.com/will-rowe/drax drax-repo && cd drax-repo/tests

# set up dummy refData as real data too large for build tests
mkdir DRAX-data && cd DRAX-data
groot get -d arg-annot
mv arg-annot.90 grootDB
cd ..
wget -q -O arg-db.fna https://github.com/will-rowe/groot/blob/master/db/full-ARG-databases/arg-annot-db/arg-annot-db.fna
samtools faidx  arg-db.fna

# run tests
run_name="Test drax Run: "$(date +%s)
cmd="nextflow run ../main.nf -name \"${run_name}\" --max_cpus 1 --max_memory '7300.MB' --reads \"./*R{1,2}.fastq.gz\" --refData ./DRAX-data --decontaminate false  --noTaxa true"
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
