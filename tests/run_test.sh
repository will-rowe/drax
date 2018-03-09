#!/usr/bin/env bash
git clone https://github.com/will-rowe/drax && cd drax/tests
run_name="Test drax Run: "$(date +%s)
cmd="nextflow run ../main.nf -name \"${run_name}\" --max_cpus 1 --max_memory '7.GB' --reads \"*R{1,2}.fq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
