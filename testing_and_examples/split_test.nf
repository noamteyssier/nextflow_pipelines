
params.reads = "samples/sub1*_R1.fastq.gz"
params.regions = "/home/noam/data/intervals/pf3d7_misc/core.regions"
params.outdir = "out_split"

regions = file(params.regions)

Channel
    .fromPath( params.reads )                                             
    .set { sample_ids } 


process split_interval {
	publishDir "$params.outdir/variants/${sample}"

	input:
	file sample from sample_ids
	each region from regions.readLines()

	output:
	file "${region.replace(':', '-')}.${sample.getSimpleName()}.txt" into sample_merge

	script:
	"""
	#!/usr/bin/env bash

	echo ${sample} > ${region.replace(':', '-')}.${sample.getSimpleName()}.txt
	"""
}

sample_merge
  .map { 
  	file -> tuple(file.getSimpleName(), file) 
  }
  .groupTuple()
  .set { grouped_files }

process merge_samples {
	publishDir "$params.outdir/variants/merged"

	input:
	set region, file(id) from grouped_files

	output:
	file "${region}.output" into merged_regions

	script:
	"""
	#!/usr/bin/env bash

	echo ${region} > ${region}.output
	cat ${id} >> ${region}.output
	"""
}

process merge_regions {
	publishDir "$params.outdir/"

	input:
	file all_regions from merged_regions.collect()

	output:
	file "concat.txt"

	script:
	"""
	#!/usr/bin/env bash

	cat ${all_regions} > concat.txt
	"""


}

