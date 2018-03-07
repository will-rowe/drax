#!/usr/bin/env bash
git clone https://github.com/will-rowe/drax && cd drax/tests

# get a small dummy genome to use as a reference for decontamination step of QC
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

# index it for BBMap
wget https://sourceforge.net/projects/bbmap/files/BBMap_37.90.tar.gz/download -O bbmap.tar
tar -xvf bbmap.tar
sh bbmap/bbmap.sh ref=GCF_000005845.2_ASM584v2_genomic.fna.gz

# get NextFlow
curl -s https://get.nextflow.io | bash

# run the test
run_name="Test drax Run: "$(date +%s)
cmd="./nextflow run ../main.nf -name \"${run_name}\" -profile docker --max_cpus 1 --max_memory '2.GB' --max_time '1.h' --reads \"*R{1,2}.fq.gz\" --path2ref $PWD"
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
