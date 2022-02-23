#!/bin/sh

## downloaded files from SRA

bowtie2 --maxins 2000 -x /Volumes/tbDrive_kvdb/atacseq/hg19_botwie2Index/hg19  -1 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650766_1.fastq  -2 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650766_2.fastq -S /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650766.sam
# high alignment rate
/Applications/samtools/bin/samtools view -bo SRR7650766.bam SRR7650766.sam

bowtie2 --maxins 2000 -x /Volumes/tbDrive_kvdb/atacseq/hg19_botwie2Index/hg19  -1 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650808_1.fastq  -2 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650808_2.fastq -S /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650808.sam
# 97.32% overall alignment rate
/Applications/samtools/bin/samtools view -bo SRR7650808.bam SRR7650808.sam

bowtie2 --maxins 2000 -x /Volumes/tbDrive_kvdb/atacseq/hg19_botwie2Index/hg19  -1 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650848_1.fastq  -2 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650848_2.fastq -S /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650848.sam
# 55.95% overall alignment rate
/Applications/samtools/bin/samtools view -bo SRR7650848.bam SRR7650848.sam

bowtie2 --maxins 2000 -x /Volumes/tbDrive_kvdb/atacseq/hg19_botwie2Index/hg19  -1 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650885_1.fastq  -2 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650885_2.fastq -S /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650885.sam
# 95.51% overall alignment rate
/Applications/samtools/bin/samtools view -bo SRR7650885.bam SRR7650885.sam

bowtie2 --maxins 2000 -x /Volumes/tbDrive_kvdb/atacseq/hg19_botwie2Index/hg19  -1 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650913_1.fastq  -2 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650913_2.fastq -S /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650913.sam
# 97.78% overall alignment rate
/Applications/samtools/bin/samtools view -bo SRR7650913.bam SRR7650913.sam

bowtie2 --maxins 2000 -x /Volumes/tbDrive_kvdb/atacseq/hg19_botwie2Index/hg19  -1 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650922_1.fastq  -2 /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650922_2.fastq -S /Volumes/tbDrive_kvdb/atacseq/calderon/SRR7650922.sam
# 97.88% overall alignment rate
/Applications/samtools/bin/samtools view -bo SRR7650922.bam SRR7650922.sam


## sort BAM files
/Applications/samtools/bin/samtools sort -o SRR7650766.sorted.bam SRR7650766.bam
/Applications/samtools/bin/samtools sort -o SRR7650808.sorted.bam SRR7650808.bam
/Applications/samtools/bin/samtools sort -o SRR7650848.sorted.bam SRR7650848.bam
/Applications/samtools/bin/samtools sort -o SRR7650885.sorted.bam SRR7650885.bam
/Applications/samtools/bin/samtools sort -o SRR7650913.sorted.bam SRR7650913.bam
/Applications/samtools/bin/samtools sort -o SRR7650922.sorted.bam SRR7650922.bam

## filter mapping qualities <30 and use flags reproted in paper
/Applications/samtools/bin/samtools view -bq 30 -F 1804 -f 2 SRR7650766.sorted.bam > SRR7650766.sorted.filtered.bam
/Applications/samtools/bin/samtools view -bq 30 -F 1804 -f 2 SRR7650808.sorted.bam > SRR7650808.sorted.filtered.bam
/Applications/samtools/bin/samtools view -bq 30 -F 1804 -f 2 SRR7650848.sorted.bam > SRR7650848.sorted.filtered.bam
/Applications/samtools/bin/samtools view -bq 30 -F 1804 -f 2 SRR7650885.sorted.bam > SRR7650885.sorted.filtered.bam
/Applications/samtools/bin/samtools view -bq 30 -F 1804 -f 2 SRR7650913.sorted.bam > SRR7650913.sorted.filtered.bam
/Applications/samtools/bin/samtools view -bq 30 -F 1804 -f 2 SRR7650922.sorted.bam > SRR7650922.sorted.filtered.bam

## remove duplicate reads using Picard tool
java -jar /Applications/picardtools/picard.jar MarkDuplicates --INPUT SRR7650766.sorted.filtered.bam \
--OUTPUT SRR7650766.sorted.filtered.rmDup.bam \
--METRICS_FILE SRR7650766.dupMetrics \
--REMOVE_DUPLICATES true
java -jar /Applications/picardtools/picard.jar MarkDuplicates --INPUT SRR7650808.sorted.filtered.bam \
--OUTPUT SRR7650808.sorted.filtered.rmDup.bam \
--METRICS_FILE SRR7650808.dupMetrics \
--REMOVE_DUPLICATES true
java -jar /Applications/picardtools/picard.jar MarkDuplicates --INPUT SRR7650848.sorted.filtered.bam \
--OUTPUT SRR7650848.sorted.filtered.rmDup.bam \
--METRICS_FILE SRR7650848.dupMetrics \
--REMOVE_DUPLICATES true
java -jar /Applications/picardtools/picard.jar MarkDuplicates --INPUT SRR7650885.sorted.filtered.bam \
--OUTPUT SRR7650885.sorted.filtered.rmDup.bam \
--METRICS_FILE SRR7650885.dupMetrics \
--REMOVE_DUPLICATES true
java -jar /Applications/picardtools/picard.jar MarkDuplicates --INPUT SRR7650913.sorted.filtered.bam \
--OUTPUT SRR7650913.sorted.filtered.rmDup.bam \
--METRICS_FILE SRR7650913.dupMetrics \
--REMOVE_DUPLICATES true
java -jar /Applications/picardtools/picard.jar MarkDuplicates --INPUT SRR7650922.sorted.filtered.bam \
--OUTPUT SRR7650922.sorted.filtered.rmDup.bam \
--METRICS_FILE SRR7650922.dupMetrics \
--REMOVE_DUPLICATES true

## index BAM files
/Applications/samtools/bin/samtools index SRR7650766.sorted.filtered.rmDup.bam



# on SCF
# ssh koenvdberge@gandalf.berkeley.edu
# module load R/4.1.0

## in R
#BiocManager::install("gcapc", lib = "/accounts/campus/koenvdberge/rpkgs/")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", lib = "/accounts/campus/koenvdberge/rpkgs/")


#  preprocessCalderon.sh
#  
#
#  Created by Koen Van den Berge on 10/3/21.
#  
