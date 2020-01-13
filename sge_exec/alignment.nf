
// Script Parameters

ref = file(params.reference)
known_sites = file(params.known_sites)
interval_list = file(params.interval_list)


clusterOptions_single = "-S /bin/bash -cwd -q long.q -j y -l scratch=75G -l mem_free=4G"
clusterOptions_multi = clusterOptions_single.concat(" -pe smp ${params.threads}")

/*
Create 'read_pairs' channel that emits for each read pair a 
tuple containing 3 elements: pair_id, R1, R2
*/
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs } 


// fastp based quality filtering and QC
process fastp {

	publishDir "${params.outdir}/$pair_id", 
		saveAs: {filename ->
			if (filename.indexOf(".fastq.gz") > 0) "filtered/$filename"
			else if (filename.indexOf("_fastp") > 0) "metrics/$filename"
			else filename
		}

	clusterOptions = clusterOptions_multi.concat(" -o qc.log")

	input:
	set pair_id, file(reads) from read_pairs

	output:
	set pair_id, "filt_${pair_id}_R{1,2}.fastq.gz" into filtered_reads
	file("${pair_id}_fastp.{json,html}")

	script:
	"""
	#!/usr/bin/env bash

	fastp \
		-i ${reads[0]} \
		-I ${reads[1]} \
		-o "filt_${pair_id}_R1.fastq.gz" \
		-O "filt_${pair_id}_R2.fastq.gz" \
		-w ${params.threads} \
		-j "${pair_id}_fastp.json" \
		-h "${pair_id}_fastp.html" \
		--trim_poly_g
	
	"""
}


// bwa based alignment
process bwa_align {

	publishDir "${params.outdir}/$pair_id"

	clusterOptions = clusterOptions_multi.concat(" -o align.log")
	
	input:
	set pair_id, file(reads) from filtered_reads

	output:
	set pair_id, "${pair_id}.bam" into rg_bam_list

	script:
	
	"""
	#!/usr/bin/env bash

	# alignment
	bwa mem -t ${params.threads} $ref ${reads} | \
		samtools view -b -o temp.bam

	# sorting reads
	samtools sort -@ ${params.threads} -o sorted.bam temp.bam

	# adding read groups
	gatk AddOrReplaceReadGroups \
		-I sorted.bam \
		-O ${pair_id}.bam \
	    --RGLB 1 --RGPL illumina \
	    --RGPU 1 --RGSM ${pair_id}

	# remove intermediary bams
	rm temp.bam
	rm sorted.bam

	"""
}


// gatk base recalibrator
process base_recal {

	publishDir "${params.outdir}/$pair_id",
		saveAs: {filename ->
			if (filename.indexOf("_recal.bam") > 0) filename
			else if (filename.indexOf(".recal") > 0) "metrics/$filename"
			else null
		}

	clusterOptions = clusterOptions_multi.concat(" -o recal.log")

	input:
	set pair_id, file(rg_bam) from rg_bam_list

	output:
	set pair_id, "${pair_id}_recal.bam" into recal_bam_list
	file("${pair_id}.recal")
	file("${pair_id}_recal.bam.bai")

	script:
	"""
	#!/usr/bin/env bash

    gatk BaseRecalibrator \
        -I ${rg_bam} \
        -O ${pair_id}.recal \
        -R ${ref} \
        --known-sites ${known_sites}

    gatk ApplyBQSR \
        -bqsr ${pair_id}.recal \
        -I ${rg_bam} \
        -O ${pair_id}_recal.bam

    samtools index -@ {threads} ${pair_id}_recal.bam

	"""
}


// flagstat metrics and coverage statistics
process metrics {

	publishDir "${params.outdir}/$pair_id/metrics"

	clusterOptions = clusterOptions_multi.concat(" -o metrics.log")

	input:
	set pair_id, file(recal_bam) from recal_bam_list

	output:
	file("${pair_id}.flagstat")
	file("${pair_id}.coverage")

	script:
	"""

	# gather flagstat statistics
	samtools flagstat -@ ${params.threads} ${recal_bam} | \
        ~/bin/binftools/flagstat_parser/flagstat_parser.py \
            -n ${pair_id} > ${pair_id}.flagstat
	
    # gather coverage statistics
    gatk CollectWgsMetrics \
        -I ${recal_bam} \
        -R ${ref} \
        -O ${pair_id}.coverage \
        --INTERVALS ${interval_list}

	"""
}
