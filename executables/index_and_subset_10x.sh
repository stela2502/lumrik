samtools index /tmp/pbmc_1k.bam
samtools view -b /tmp/pbmc_1k.bam GL000220.1 > /tmp/GL000220.1.bam
samtools index /tmp/GL000220.1.bam
