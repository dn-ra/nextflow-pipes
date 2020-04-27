#!/bin/bash -e


#basecalled fastq file
READS='/DataOnline/Data/raw_external/Coronavirus/Direct_RNA_Sequence_Cellular/20200211_1144_GA50000_FAL90649_37c13261/fastq_pass/Cells.fastq'
HOST=''
VIRUS=''

#minimap2 location
MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'


echo '>Leader
ACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAAC' > leader.fa


#map to whole genome
#$MM2 -ax splice -un /home/daniel/wuhan_coronavirus_australia.fasta $READS | samtools view -hF4 -b | samtools sort - > sorted.whole_genome_mapped.bam

#map to leader. Modify -k paramater as fits
#$MM2 -ax map-ont -un leader.fa $READS | samtools view -hF4 -b | samtools sort - > sorted.leader_mapped.bam

#map to transcriptome. Primary alignments only.
#$MM2 -ax map-ont -un -N 0 $HOST $READS | samtools view -hF4 -b | samtools sort - > sorted.host_mapped.bam

#extract lengths of reads with leader
seqtk subseq $READS <(samtools view leader_mapped.bam | cut -f 1 -| sort | uniq) > leader_reads.fq
seqtk comp leader_reads.fq | tee read_comp.txt | cut -f 2 > leader_read_lengths.txt

#R script to find main subgenome mRNAs
./find_subgenomes.R leader_read_lengths.txt > bin_bounds.txt

#get reads for each length bin
#TODO - change temp_file to pipe?

i=1
while read line; do bounds=( $line ); file_name="bin${i}_${bounds[0]}_${bounds[1]}.fq" ; \
temp_file=$(mktemp); awk -v min=${bounds[0]} -v max=${bounds[1]} '{if ($2 >= min && $2 <= max) print $1}' read_comp.txt > $temp_file; \
seqtk subseq leader_reads.fq $temp_file > $file_name; rm $temp_file; i=$(expr $i + 1); \
done < bin_bounds.txt

#minimap for each separate bin here

