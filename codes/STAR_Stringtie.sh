#!usr/bin/sh

# Calculation of RNAseq counts by SATR and Stringtie

FATSQ_PATH=path/to/fastq
STAR_INDEX_PATH=path/to/STARindex
GTF_PATH=path/to/gtf
STAR_OUTPUT_PATH=path/to/star_output
ballgown_OUTPUT_PATH=path/to/ballgown_output


cd ${STAR_OUTPUT_PATH}
mkdir -p SampleName
cd SampleName

STAR --runThreadN 6 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts \
     --genomeDir ${STAR_INDEX_PATH} \
     --outSAMstrandField intronMotif \
     --outFilterMultimapNmax 20 \
     --sjdbGTFfile ${GTF_PATH} \
     --readFilesIn ${FASTQ_PATH}/filteredfastq_rRNA_unmapped.1.fastq.gz ${FASTQ_PATH}/filteredfastq_rRNA_unmapped.2.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix SampleName_

mkdir -p ${ballgown_OUTPUT_PATH}/SampleName

stringtie -p 6 \
SampleName_Aligned.sortedByCoord.out.bam -e -B -G ${GTF_PATH} \
-o ${ballgown_OUTPUT_PATH}/SampleName/SampleName.gtf \
-A ${ballgown_OUTPUT_PATH}/SampleName/SampleName_gene_abund.tab

##  After all samples have been computed

cd ${ballgown_OUTPUT_PATH}

python3 prepDE.py3


