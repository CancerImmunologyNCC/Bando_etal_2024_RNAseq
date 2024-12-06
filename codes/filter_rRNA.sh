#!usr/bin/sh

# filtering rRNA reads

FASTQ_PATH=path/to/fastq
rRNA_bowtie2_index_PATH=path/to/index 
OUTPUT_PATH=path/to/output

cd ${FASTQ_PATH}
bowtie2 -x ${rRNA_bowtie_index_PATH} -1 trimmedfastq_1_val_1.fq.gz -2 trimmedfastq_2_val_2.fq.gz \
-S ${OUTPUT_PATH}/output_rRNA.mapped.sam --un-conc ${OUTPUT_PATH}/output_rRNA_unmapped.fastq

cd ${OUTPUT_PATH}
rm output_rRNA.mapped.sam
gzip output_rRNA_unmapped.1.fastq
gzip output_rRNA_unmapped.2.fastq