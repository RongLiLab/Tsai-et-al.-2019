#!/usr/bin/env bash


#Load FastQC
module load fastqc

#RNASeq directory
cd ~/SeqHome/RNASeq

#Pre-alignment quality check
mkdir RNASeq_QC
for file in *.fastq ; do  fastqc -f ${file} -o /RNASeq_QC/${file} ; done

#Build HISAT2 index
hisat2-build –p 4 –f referencegenome.fa 

#Align fastq files with HISAT2
mkdir SAMfiles
for f in *.fastq ; do  hisat2 –p 4 -x referencegenome –q –U $f –S /SAMfiles/${f}.sam ; done

#Index for samtools
samtools faidx referencegenome.fa

#Convert .sam files to indexed and sorted .bam
cd SAMfiles
for i in *.sam ; do samtools view -@ 4 –t referencegenome.fa.fai  –bS $i | samtools sort -$i.sorted && echo $i "bam sorted" && samtools index $i.sorted.bam; done

#Use HTSeq to get read count
mkdir RNASeq_Reads
for j in *.sorted.bam; do htseq-count –f = bam -s = reverse --type= mRNA --idattr = Parent $j reference.gff  > /RNASeq_Reads/${j}.txt ; done

