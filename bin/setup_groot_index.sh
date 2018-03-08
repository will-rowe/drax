#!/usr/bin/env bash
# $1 --- input file from seqkit stats
# $2 --- cpus
# $3 --- output directory for index

# delete the header line of seqkit stats output and get the average length column
tail -n +2 $1 | cut -f 7 >> averageReadLength.txt
# get the mean and stdev
meanRL=$(awk '{ sum +=$1; n++ } END { if (n > 0) printf "%3.0f", sum /n}' averageReadLength.txt)
stdev=$(awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)**2)}' averageReadLength.txt)

# check that we can generate a GROOT index suitable for all samples
if [ $stdev -gt 10 ]; then
    echo "The mean read length of the cleaned data is too variable to run GROOT on all samples!"; exit $ERRCODE;
fi

# set up the commands
#download an ARG database for groot
grootGetCMD="groot get -d resfinder"
# index the database
grootIndexCMD="groot index -i resfinder.90 -o groot-index -l $meanRL -p $2 -o $3"

# run the commands
exec $grootGetCMD 2>&1 | tee setup_groot_index.log
exec $grootIndexCMD 2>&1 | tee setup_groot_index.log

rm averageReadLength.txt
