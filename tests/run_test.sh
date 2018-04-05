#!/usr/bin/env bash
git clone https://github.com/will-rowe/drax drax-repo && cd drax-repo/tests


# run tests
run_name="Test drax Run: "$(date +%s)
cmd="nextflow run ../main.nf -name \"${run_name}\" --max_cpus 1 --max_memory '7300.MB' --reads \"./*R{1,2}.fastq.gz\" --refData ./ --runTest true"
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
