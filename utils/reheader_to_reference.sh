#!/bin/bash -

# reheader_to_reference --- reheader a bam file and keep only reference sequences that match regex

USAGE="USAGE reheader_to_reference -i input.bam -o output.bam -r regex_to_match"

if [ $# == 0 ]; then
    echo $USAGE
    exit 1;
fi

while getopts "i:o:r:" opt
do
    case $opt in
	i)
	    inputFile=$OPTARG ;;
	o)
	    outputFile=$OPTARG ;;
	r)
	    refKeep=$OPTARG ;;
    esac
done


#    inputFile="EV0004.H004.fresh.a1.std.EX2.media_88_Novaseq_20190829.bam_EBOV.sorted.bam"
#    refKeep="KU"
#    outputFile='output.bam'


tmpfile=$(mktemp /tmp/temp_header.XXXXXX)
trap "rm -f $tmpfile" EXIT

# Make the header
samtools view -H $inputFile | awk -v refKeep="$refKeep" '!/^@SQ/ { print } /@SQ/ && $2~refKeep { print }' > ${tmpfile}

# Reheader File
cat ${tmpfile} <( samtools view $inputFile ) | samtools view -b > ${outputFile}
