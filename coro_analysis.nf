#!/usr/bin/env nextflow

/*	Author: Daniel Rawlinson, Coin Group, The Peter Doherty Institute for Infection & Immunity
*	Email:	daniel.rawlinson@unimelb.edu.au
*		This nextflow pipeline automates the processing of Direct-RNA-Sequenced SARS-CoV-2 samples for mRNA detection, 5mC methylation,
*		and host RNA identification, as employed in:
			Taiaroa G, et al. Direct RNA sequencing and early evolution of SARS-CoV-2. 2020. BioRxiv.
*/


//--------inputs------------///
params.virus_reference = '/home/daniel/wuhan_coronavirus_australia.fasta'
//params.host_reference = '/DataOnline/Data/raw_external/Coronavirus/monkey/newdb/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz'
params.host_reference = ''
params.leader = '/DataOnline/Data/raw_external/Coronavirus/Direct_RNA_Sequence_Cellular/doherty_analysis/leader.fa'
params.reads = ''
println(params.reads)

host_map="ON"

//--------evaluate inputs--------///

if (params.reads == '') exit 1, 'Please provide a one or more fastq read files with --reads []'


if (params.host_reference == '') { log.warn 'Host genome not provided. Analysis will proceed without host genome mapping.'
			host_map = "OFF" } 
/* DO THIS IN THE HOST MAP SECTION?
else if (host_map != "OFF") host = file(params.host_reference) {
if (!host.exists()) exit 1, 'Expected host mapping file does not exist: ${host}'exit 1, 'Expected host mapping file does not exist: ${host}'
}
*/

virus = file(params.virus_reference)
leader = file(params.leader)
if (!virus.exists()) exit 1, 'SARS-CoV-2 reference file does not exist: ${virus}'
if (!leader.exists()) exit 1, 'Leader sequence reference file does not exist: ${params.leader}'


//-------launch reads channel------///
/* will need to set this up to handle folders of multiple fast5s */

Channel 
	.fromPath(params.reads, checkIfExists: true) 
	.collect()
	.subscribe { println it }
	.into {fastq_reads} 

//--------output directories---------///

//-------demux???-------------------///

//--------minimap for virus----------///

process whole_genome_map {
	input:
	file file(read_file) from fastq_reads
	file reference_file
	
	output:
	file "virus_map.bam" into minimap_out

	script:
	'''
	/sw/minimap2/minimap2-2.11_x64-linux/minimap2 -ax splice -un -k 14 $reference_file $read_file | samtools view -hb > minimap_out.bam
	
	'''
	}
process leader_map {
	input:
	file read_file from ___
	file(leader)

	output:
	file "leader_map.bam" into detect_subgenomes_from_leader
	
	script:
	'''
	/sw/minimap2/minimap2-2.11_x64-linux/minimap2 -ax map-ont -k 14 -un $leader_reference $read_file | samtools view -hb > minimap_out.bam
	'''
	}
//-------minimap for host---------------///

process host_map {
	input:
	file read_file from input_channel
	file(host)

	output:
	
	script:
	'''
	'''
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
