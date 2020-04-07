#!/usr/bin/env nextflow

/*	Author: Daniel Rawlinson, Coin Group, The Peter Doherty Institute for Infection & Immunity
*	Email:	daniel.rawlinson@unimelb.edu.au
*		This nextflow pipeline automates the processing of Direct-RNA-Sequenced SARS-CoV-2 samples for mRNA detection, 5mC methylation,
*		and host RNA identification, as employed in:
			Taiaroa G, et al. Direct RNA sequencing and early evolution of SARS-CoV-2. 2020. BioRxiv.
*/


//--------inputs------------
params.virus_reference = '/home/daniel/wuhan_coronavirus_australia.fasta'
params.host_reference = '//DataOnline/Data/raw_external/Coronavirus/monkey/newdb/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz'
params.leader = '/DataOnline/Data/raw_external/Coronavirus/Direct_RNA_Sequence_Cellular/doherty_analysis/leader.fa'
params.reads = ''


//--------check inputs--------

host = file(params.host_reference)
virus = file(params.virus_reference)
leader = file(params.leader)

if (!host.exists()) println 'Warning: Host genome not provided. Analysis will proceed without host genome mapping.'
if (!virus.exists()) exit 1, 'SARS-CoV-2 reference file does not exist: ${virus}'
if (!leader.exists()) exit 1, 'Leader sequence reference file does not exist: ${leader}'

//-------launch reads channel------

input_channel = Channel.fromPath(params.reads)

//--------minimap for virus----------///

process whole_genome_map {
	input:
	file read_file from input_channel
	file(virus)
	
	output:
	file "virus_map.bam" into minimap_out

	script:
	'''
	/sw/minimap2/minimap2-2.11_x64-linux/minimap2 -ax splice -un -k 14 $reference_file $read_file | samtools view -hb > minimap_out.bam
	
	'''
	}
process leader_map {
	input:
	file read_file from input channel
	file(leader)

	output:
	file "leader_map.bam" into detect_subgenomes_from_leader
	
	script:
	'''
	/sw/minimap2/minimap2-2.11_x64-linux/minimap2 -ax map-ont -k 14 -un $leader_reference $read_file | samtools view -hb > minimap_out.bam
	'''
	
//-------minimap for host---------------///

process host_map {
	input:
	file read_file from input_channel
	file(host)

	output:
	
	script:
	
}


//-------find mRNAs--------------------///



//--Breakpoint analysis---------------///





/*##map to genome first
*
*##extract reads to leader
*##read-length histogram
*##smoothed z-score 
*
*
*##breakpoint analysis
*
*
*##meth analysis
*
*
*##will need:
*#	japsa, smoothed-z score function, 
*/