#!/usr/bin/env nextflow

params.reads='s3://algae.chk14.transcriptome/CHK14*_R{1,2}_001.fastq.gz'
pairInt='s3://transcriptomepipeline/PairInterleaves.sh'
sequences2='s3://transcriptomepipeline/ContaminantsForRemove.fasta'
sequences22='s3://transcriptomepipeline/ContaminantsForRemove.fasta'
adapters='s3://transcriptomepipeline/TruSeq3-PE.fa'
pairInt2='s3://transcriptomepipeline/PairInterleaves.sh'



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

	memory '4G'

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




process join1 {

	memory '4G'

	input:
	tuple val(pair_id), path(fileList) from reads_ch.collect()
	
	output:
	file 'R1reads.fa' into R1reads
	
	"""
	cat `ls * ` > R1reads.fa
	"""
}


process join2 {

	memory '4G'

	input:
	tuple val(pair_id), path(fileList) from reads_ch2.collect()
	
	output:
	file 'R2reads.fa' into R2reads
	
	"""
	cat `ls *` > R2reads.fa
	"""
}


process minimapS31 {
	
	memory '16G'
	
	input:
	path fastq from R1reads
	path contam from sequences2
	
	output:
	file 'aln.sam.gz' into align
	

    """
	minimap2 -a $contam $fastq > aln.sam; gzip aln.sam
    """

}


process samtools1 {

	memory '12G'
	
	input:
	path samalign from align
	
	output:
	file 'clean1.tst.fastq' into cleanReads1
	
	"""
	gunzip -f $samalign; samtools fastq -n -f 4 aln.sam > clean1.tst.fastq
	
	"""

}

process cutadapt1 {

	memory '10G'
	
	input:
	path cleanfas from cleanReads1
	
	output:
	file 'R1.fastq' into reads1
	
	"""
	cutadapt --rename='{id}/1' $cleanfas -j 7 -o R1.fastq
	"""

}




//Second  pair

process minimapS32 {

	memory '16G'
	
	input:
	path fastq from R2reads
	path contam from sequences22
	
	output:
	file 'aln.sam.gz' into align2
	

    """
	minimap2 -a $contam $fastq > aln.sam; gzip aln.sam
    """

}


process samtools12 {

	memory '8G'
	
	input:
	path samalign from align2
	
	output:
	file 'clean1.tst.fastq' into cleanReads12
	
	"""
	gunzip -f $samalign ; samtools fastq -n -f 4 aln.sam > clean1.tst.fastq
	
	"""

}



process cutadapt12 {

	memory '10G'
	
	input:
	path 'cleanfas' from cleanReads12
	
	output:
	file 'R2.fastq' into reads12
	
	"""
	cutadapt --rename='{id}/2' $cleanfas -j 7 -o R2.fastq
	"""

}

process fastqpair {

	memory '16G'

	input:
	path 'R1fastq' from reads1
	path 'R2fastq' from reads12

	output:
	file 'R1fastq.paired.fq' into pairR1
	file 'R2fastq.paired.fq' into pairR2
	file 'R1fastq.single.fq' into pairR3
	file 'R2fastq.single.fq' into pairR4

	"""
	fastq_pair -t 10000000 $R1fastq $R2fastq
	"""

}

process Trimmomatic {

	memory '1G'

	input:
	path 'R1pair' from pairR1
	path 'R2pair' from pairR2
	path 'adapt' from adapters


	output:
	file 'R1p.fq' into readTrim1
	file 'R2p.fq' into readTrim2

	"""
	trimmomatic PE -threads 12 $R1pair $R2pair R1p.fq R1up.fq R2p.fq R2up.fq ILLUMINACLIP:$adapt:2:30:10 SLIDINGWINDOW:4:20
	"""

}



process bbnorm {

	memory '16G'
	
        input:
        path seq1 from readTrim1
        path seq2 from readTrim2
        
        output:
        file 'mid.fq' into ReadTrimNorm1

	"""
	bbnorm.sh in=$seq1 in2=$seq2 outlow=low.fq outmid=mid.fq outhigh=high.fq passes=1 lowbindepth=6 highbindepth=60 -Xmx15g
	"""
}


process pairInt2 {

	memory '1G'

	input:
	path 'pairInt' from pairInt2
	path 'Intpair' from ReadTrimNorm1

	output:
	file 'R1reads.fastq' into R1Tofastq
	file 'R2reads.fastq' into R2Tofastq

	"""
	chmod 744 $pairInt
	./$pairInt < $Intpair R1reads.fastq R2reads.fastq
	"""

}


process fastqpair2 {

	memory '2G'

	input:
	path R1p from R1Tofastq
	path R2p from R2Tofastq

	output:
	file 'R1reads.fastq.paired.fq' into pairR1T
	file 'R2reads.fastq.paired.fq' into pairR2T
	//For now not even bothering with unpaired

	"""
	fastq_pair -t 10000000 $R1p $R2p
	"""
}

pairR1T.into{P1NormSpades; P1NormTrinity}
pairR2T.into{P2NormSpades; P2NormTrinity}

process SpadeAssemble {

	memory '24G'

        input:
        path R1Norm from P1NormSpades
	path R2Norm from P2NormSpades

        output:
        file 'hard_filtered_transcripts.fasta' into Spades

        """
        rnaspades.py  --pe1-1 $R1Norm --pe1-2 $R2Norm  -o spades_output
	cp ./spades_output/hard_filtered_transcripts.fasta .
        """
}






process TrinityAssemble {

	memory '64G'

	input:
	path R1pair from P1NormTrinity
	path R2pair from P2NormTrinity

	output:
	file 'Trinity.fasta' into Trinity


	"""
	conda install tbb=2020.2
	Trinity --seqType fq --left $R1pair --right $R2pair --max_memory 56G --output trinity_output
	cp ./trinity_output/Trinity.fasta .
	"""

}












