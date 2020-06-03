
// Script Parameters


params.indir = "recal_bams"

params.outdir = "results"

ref = file(
	params.reference
	)

known_sites = file(
	params.known_sites
	)

interval_list = file(
	params.interval_list
	)

regions = file(
	params.regions
	)

clusterOptions_single = "-S /bin/bash -cwd -q long.q -j y -l scratch=75G"
clusterOptions_multi = clusterOptions_single.concat(" -pe smp ${params.threads}")


Channel
    .fromFilePairs( params.indir + "/**{.bam,.bam.bai}" )                                             
    .set { recal_for_gvcf } 


// call gvcf on intervals independently
process call_gvcf {
	
	publishDir "$params.outdir/variants/intervals/${region}"

	clusterOptions = clusterOptions_multi.concat(" -o call_gvcf.log")


	input:
	set sample_id, file(recal_bam) from recal_for_gvcf
	each region from regions.readLines()

	output:
	file "${region.replace(':', '-')}.${sample_id}.g.vcf" into gvcf_list
	file "${region.replace(':', '-')}.${sample_id}.g.vcf.idx" into gvcf_index_list

	script:
	"""
	#!/usr/bin/env bash

	# call gvcf for sample
	gatk HaplotypeCaller \
		-R ${ref} \
		-I ${recal_bam[0]} \
		--ERC GVCF \
		-L ${region} \
		-O ${region.replace(':', '-')}.${sample_id}.g.vcf

	# index gvcf
	gatk IndexFeatureFile -F ${region.replace(':', '-')}.${sample_id}.g.vcf

	"""
}

// // Create channel for interval and all gvcfs related
// gvcf_list
//   .map { 
//   	gvcf -> tuple(gvcf.getSimpleName(), gvcf) 
//   }
//   .groupTuple()
//   .set { grouped_vcfs }

// // Create channel for interval and all indices related
// gvcf_index_list
//   .map { 
//   	gvcf -> tuple(gvcf.getSimpleName(), gvcf) 
//   }
//   .groupTuple()
//   .set { grouped_idx }

// // genotype each interval independently
// process joint_genotype {

// 	publishDir "$params.outdir/variants/merged/${region}"

// 	clusterOptions = clusterOptions_multi.concat(" -o joint_geno.log")


// 	input:
// 	set region, file(gvcf) from grouped_vcfs
// 	set region2, file(idx) from grouped_idx

// 	output:
// 	file "${region}.merged.vcf" into merged_regions
// 	file "${region}.merged.vcf.idx" into merged_regions_idx
// 	file "${region}_gendb"

// 	script:
// 	"""
// 	#!/usr/bin/env bash

// 	# Load gvcfs into a gendb 
// 	gatk GenomicsDBImport \
// 		--genomicsdb-workspace-path ${region}_gendb \
// 		-L ${region.replaceFirst('-' , ':')} \
// 		\$(for v in ${gvcf}; do echo "-V \$v"; done)


// 	# genotype
// 	gatk GenotypeGVCFs \
// 	    -R ${ref} \
// 	    -V gendb://${region}_gendb \
// 	    -G StandardAnnotation \
// 	    --new-qual true \
// 	    -L ${region.replaceFirst('-' , ':')} \
// 	    -O ${region}.merged.vcf 

// 	"""
// }

// // Channel of sorted regions
// merged_regions
// 	.toSortedList()
// 	.set {sorted_merged_regions}

// // Concatenated merged regions and prepare final output file
// process concatenate_vcfs {
	
// 	publishDir "$params.outdir/variants"

// 	clusterOptions = clusterOptions_multi.concat(" -o concat.log")


// 	input:
// 	file merged_region from sorted_merged_regions.collect()
// 	file merged_region_idx from merged_regions_idx.collect()

// 	output:
// 	file "sampleset.vcf.gz"
// 	file "sampleset.vcf.gz.csi"

// 	script:
// 	"""
// 	#!/usr/bin/env bash

// 	# collect and sort
// 	bcftools concat *.merged.vcf | bcftools sort -Oz -o "sampleset.vcf.gz"

// 	# index
// 	bcftools index "sampleset.vcf.gz"

// 	"""
// }