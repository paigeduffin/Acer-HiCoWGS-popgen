#!/bin/bash

ID=$1

source /home1/pduffin/miniconda3/etc/profile.d/conda.sh
conda activate samtools_env

samtools merge ${ID}_host_merged.bam ./unpaired_host/${ID}.unpaired_sorted_host.bam ./paired_host/${ID}.paired_sorted_host.bam
samtools sort ${ID}_host_merged.bam -@ 32 -O BAM -o ${ID}_host_merged.sorted.bam

# SKIPPING FOR NOW
#ml gcc/11.3.0
# ml bamtools/2.5.2
# echo $ID
# bamtools stats -in ${ID}_host_merged.sorted.bam

echo $ID
samtools depth ${ID}_host_merged.sorted.bam | awk '{sum+=$3} END {print "Average = ",sum/NR}'

###########################
#########################################################################################

### STEP 1: ADD READ GROUPS ###
# ml gcc/13.3.0 picard/2.26.2 jdk/17.0.5

ml legacy/CentOS7 gcc/11.3.0 jdk/17.0.5

export PICARD=/spack/2206/apps/linux-centos7-x86_64_v3/gcc-11.3.0/picard-2.26.2-bk5qcm2/bin/picard.jar

java -Xmx4G -jar $PICARD AddOrReplaceReadGroups        I=${ID}_host_merged.sorted.bam        O=${ID}_host_RG.bam        RGID=${ID}        RGLB=lib1        RGPL=ILLUMINA        RGPU=unit1        RGSM=${ID}

#########################################################################################

### STEP 2: SORT ###
samtools sort ${ID}_host_RG.bam > ${ID}_host_RG.sort.bam

rm ${ID}_host_RG.bam
#########################################################################################
### STEP 3: FIX MATES ###
export PICARD=/spack/2206/apps/linux-centos7-x86_64_v3/gcc-11.3.0/picard-2.26.2-bk5qcm2/bin/picard.jar

java -Xmx4G -jar $PICARD  FixMateInformation I=${ID}_host_RG.sort.bam O=${ID}_host_RG.sort.fxm8s.bam

rm ${ID}_host_RG.sort.bam
#########################################################################################
### STEP 4: MARK DUPLICATES ###
export PICARD=/spack/2206/apps/linux-centos7-x86_64_v3/gcc-11.3.0/picard-2.26.2-bk5qcm2/bin/picard.jar

java -Xmx4G -jar $PICARD  MarkDuplicates I=${ID}_host_RG.sort.fxm8s.bam O=${ID}_host_GATKd.bam M=${ID}_host_metrics.txt

rm ${ID}_host_RG.sort.fxm8s.bam

#########################################################################################
### STEP 5: FINAL SORT AND INDEX ###
samtools sort -@ 32 -O BAM -o ${ID}_GATKd_sort.bam ${ID}_host_GATKd.bam
samtools index ${ID}_GATKd_sort.bam
