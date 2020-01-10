
// Script Parameters

ref = file(params.reference)
known_sites = file(params.known_sites)
interval_list = file(params.interval_list)

Channel
	.fromPath( "bams/*.bam" )
	.map {x -> [x.name.replace('.bam', ''), x]}
	.set { input_bams }


process sort_and_addRG {

	publishDir "${params.outdir}/${sample_id}"

	input:
	set sample_id, file(bam) from input_bams

	output:
	set sample_id, "${sample_id}.rg{.bam,.bam.bai}" into rg_added
	
	script:
	"""
	#!/usr/bin/env bash

	samtools sort \
		-@ ${params.threads} \
		-o temp.sorted.bam \
		${bam}

	# adding read groups
	gatk AddOrReplaceReadGroups \
		-I temp.sorted.bam \
		-O ${sample_id}.rg.bam \
	    --RGLB 1 --RGPL illumina \
	    --RGPU 1 --RGSM ${sample_id}

	samtools index ${sample_id}.rg.bam

	"""

}


// gatk base recalibrator
process base_recal {

	publishDir "${params.outdir}/$sample_id"

	input:
	set sample_id, file(rg_bam) from rg_added

	output:
	set sample_id, "${sample_id}_recal.bam" into recal_bam_list
	file("${sample_id}.recal")
	file("${sample_id}_recal.bam.bai")

	script:
	"""
	#!/usr/bin/env bash

    gatk BaseRecalibrator \
        -I ${rg_bam[0]} \
        -O ${sample_id}.recal \
        -R ${ref} \
        --known-sites ${known_sites}

    gatk ApplyBQSR \
        -bqsr ${sample_id}.recal \
        -I ${rg_bam[0]} \
        -O ${sample_id}_recal.bam

    samtools index -@ {threads} ${sample_id}_recal.bam

	"""
}