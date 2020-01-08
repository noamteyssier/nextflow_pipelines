
// Script Parameters


ref = file(params.reference)
known_sites = file(params.known_sites)
interval_list = file(params.intervals)
regions = file(params.regions)

clusterOptions_single = "-S /bin/bash -cwd -q long.q -j y"
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

	publishDir "$params.outdir/$pair_id", 
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

	publishDir "$params.outdir/$pair_id"
	
	clusterOptions = clusterOptions_multi.concat(" -o align.log")

	input:
	set pair_id, file(reads) from filtered_reads


	output:
	set pair_id, "${pair_id}.bam" into rg_bam_list

	script:
	"""
	#!/usr/bin/env bash

	# alignment
	bwa mem -t ${params.threads} $ref ${reads} | samtools view -b -o temp.bam

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

	publishDir "$params.outdir/$pair_id",
		saveAs: {filename ->
			if (filename.indexOf("_recal.bam") > 0) filename
			else if (filename.indexOf(".recal") > 0) "metrics/$filename"
			else null
		}

	clusterOptions = clusterOptions_multi.concat(" -o recal.log")

	input:
	set pair_id, file(rg_bam) from rg_bam_list

	output:
	set pair_id, "${pair_id}_recal.bam", "${pair_id}_recal.bam.bai" into recal_for_gvcf, recal_for_metrics
	file("${pair_id}.recal")
	// file("${pair_id}_recal.bam.bai")

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

    samtools index -@ ${params.threads} ${pair_id}_recal.bam

	"""
}


// flagstat metrics and coverage statistics
process metrics {

	publishDir "$params.outdir/$pair_id/metrics"

	clusterOptions = clusterOptions_multi.concat(" -o metrics.log")

	input:
	set pair_id, file(recal_bam), file(recal_bai) from recal_for_metrics

	output:
	file("${pair_id}.flagstat")
	file("${pair_id}.coverage")

	script:
	"""
	#!/usr/bin/env bash

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

process call_gvcf {
	
	publishDir "$params.outdir/$pair_id"

	clusterOptions = clusterOptions_multi.concat(" -o call_gvcf.log")

	input:
	set pair_id, file(recal_bam), file(recal_bai) from recal_for_gvcf

	output:
	file "${pair_id}.g.vcf" into gvcf_list
	file "${pair_id}.g.vcf.idx" into gvcf_idx_list

	script:
	"""
	#!/usr/bin/env bash

	# call gvcf for sample
	gatk HaplotypeCaller \
		-R ${ref} \
		-I ${recal_bam} \
		--ERC GVCF \
		-L ${interval_list} \
		-O ${pair_id}.g.vcf

	# index genomic vcf
	gatk IndexFeatureFile -F ${pair_id}.g.vcf

	"""
}

process joint_genotype {

	publishDir "$params.outdir/variants"

	clusterOptions = clusterOptions_multi.concat(" -o joint_genotyping.log")

	input:
    file sample_gvcf from gvcf_list.collect()
    file sample_index from gvcf_idx_list.collect()

	output:
	file "variant_set.vcf"

	script:
	"""
	#!/usr/bin/env bash

	# Load gvcfs into a gendb 
	gatk GenomicsDBImport \
		--genomicsdb-workspace-path gendb \
		\$(while read interval; do echo "-L \$interval"; done < ${regions}) \
		\$(for v in \$(ls *.vcf); do echo "-V \$v"; done)	

	# genotype
	gatk GenotypeGVCFs \
	    -R ${ref} \
	    -V gendb://gendb \
	    -G StandardAnnotation \
	    --new-qual true \
	    -O variant_set.vcf

	"""
}



// 	"""
// }
