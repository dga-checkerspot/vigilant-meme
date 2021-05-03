#!/usr/bin/env nextflow

params.reads='s3://algaetranscriptomics/CHK17*_R{1,2}_001.fastq.gz'
pairInt='s3://transcriptomepipeline/PairInterleaves.sh'


Channel
	.fromFilePairs(params.reads)
	.ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
	.set { read_pairs_ch }

process pairs {

	memory '8G'
	
	input:
	tuple val(pair_id), path(reads) from read_pairs_ch
	
	output:
	set pair_id, "mid.fq" into norm_ch
	
	script:
	"""
	gunzip $reads
	bbnorm.sh in="${pair_id}_R1_001.fastq" in2="${pair_id}_R2_001.fastq" outlow=low.fq outmid=mid.fq outhigh=high.fq passes=1 lowbindepth=6 highbindepth=60 -Xmx3g
	
	"""

}



process pairInt {

	memory '1G'

	input:
	path 'pairInt' from pairInt
	tuple val(pair_id), path('midfq') from norm_ch

	output:
	tuple val(pair_id), "R1reads.fastq" into reads_ch
	tuple val(pair_id), "R2reads.fastq" into reads_ch2

	"""
	chmod 744 $pairInt
	./$pairInt < $midfq R1reads.fastq R2reads.fastq
	"""

}


process join {

	memory '4G'

	input:
	path 'seq1' from reads_ch
	path 'seq2' from reads_ch2
	
	output:
	file 'R1reads.fa' into R1reads
	file 'R2reads.fa' into R2reads
	
	"""
	cat *R1* > R1reads.fa
	cat *R2* > R2reads.fa
	"""
}





