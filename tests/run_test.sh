#!/usr/bin/env bash
curl -s https://get.nextflow.io | bash

run_name="Test drax Run: "$(date +%s)
cmd="./nextflow run ../main.nf -name \"${run_name}\" -profile docker --max_cpus 1 --max_memory '2.GB' --max_time '1.h' --reads \"*R{1,2}.fq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
