#!usr/bin/sh

# adaptor and quality trimming by trim_galore

FASTQ_PATH=path/to/fastq
OUTPUT_PATH=path/to/output

cd ${FASTQ_PATH}
trim_galore --paired fastq_1.fq.gz fastq_2.fq.gz -o ${OUTPUT_PATH}