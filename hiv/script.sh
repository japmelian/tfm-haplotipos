abyss-pe name=study k=90 in='forward.fastq reverse.fastq'
cp study-contigs.fa contigs.fasta

bwa index reference.fasta
bwa mem reference.fasta contigs.fasta > contigs.sam

#g++ fusion-overlapp.cpp
./a.out

bwa index final-contigs.fasta
bwa mem final-contigs.fasta forward.fastq reverse.fastq > contigs-mapped.sam
samtools view -bS contigs-mapped.sam > contigs-mapped.bam
samtools sort contigs-mapped.bam -o contigs-mapped-sorted.bam
samtools index contigs-mapped-sorted.bam
samtools idxstats contigs-mapped-sorted.bam