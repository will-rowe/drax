#!/usr/bin/env bash
git clone https://github.com/will-rowe/drax && cd drax/tests

# get NextFlow
curl -s https://get.nextflow.io | bash

# run the test
run_name="Test drax Run: "$(date +%s)
cmd="./nextflow run ../main.nf -name \"${run_name}\" -profile docker --reads \"*R{1,2}.fq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
