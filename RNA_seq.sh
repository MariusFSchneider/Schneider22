#!/bin/bash 


files =`ls ./*.fastq`

STAR_genome="../rsem_star_index"
STAR_GTF="Mus_musculus.GRCm38.99.gtf"
STAR_opts="--runThreadN 34 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand gunzip -c  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000"


for SRUN in $files

do
STAR $STAR_opts --genomeDir $STAR_genome --sjdbGTFfile ${STAR_GTF} --readFilesIn ../${SFILE} --outFileNamePrefix ${SRUN}.

rsem-calculate-expression --bam -p 34 \
                        ${SRUN}.Aligned.toTranscriptome.out.bam \
                        $STAR_genome/genome \
                         ${SRUN}_rsem

done
