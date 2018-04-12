#!/usr/bin/env bash


#Load FastQC
module load fastqc

#DNASeq directory
cd ~/SeqHome/DNASeq

#Pre-alignment quality check
mkdir DNASeq_QC
for file in *.fastq ; do  fastqc $file -o /DNASeq_QC ; done

#Index the reference genome
bowtie2-build referencegenome.fa

#Align fastq files with bowtie2
mkdir SAMfiles
for f in *.fastq ; do  bowtie2 –x –p 4 referencegenome –q –U $f –S /SAMfiles/${f}.sam ; done

#Index for samtools
samtools faidx referencegenome.fa

#Convert .sam files to indexed and sorted .bam
cd SAMfiles
for i in *.sam ; do samtools view -@ 4 –t referencegenome.fa.fai  –bS $i | samtools sort -$i.sorted && echo $i "bam sorted" && samtools index $i.sorted.bam; done

#Extract uniquely mappable reads
for j in *.sorted.bam; do 
samtools view –H $j  > header.sam
samtools view –F 4 $j | grep –v “XS:” | cat header.sam - | \
samtools view –b - > ${j}_unique.bam
rm header.sam
done

#Get read depth
mkdir DNASeq_UniqueReads
for file in *unique.bam; do samtools bedcov BAMSegments.bed  $file > /DNASeq_UniqueReads/${file}.txt ; done

