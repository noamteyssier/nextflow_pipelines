# Nextflow Pipelines for Data Analysis

## Requirements
```
# Must Be in Path

fastp
bwa
samtools
GATK (4.0)

```

# Recommended
```bash
# Create an analysis directory for the sample set
mkdir my_analysis
cd my_analysis
mkdir samples
mv mysamples_{1..100}_R[12].fastq.gz samples

# Identify which pipeline is necessary for the analysis
# Symlink the appropriate pipeline and config to the analysis directory
ln -s ~/nextflow_pipelines/sge_exec/alignment.nf .
ln -s ~/nextflow_pipelines/sge_exec/nextflow.config .

# run analysis
nextflow run alignment.nf

########
# Note #
########

# Make sure to adjust nextflow.config 
# to point to correct meta data directory

```

## Alignment Pipeline
QC Filtering, Alignment, Read Group addition, Base Recalibration, and alignment and sequencing Metrics calculated.

## GATK Best Practice
QC Filtering, Alignment, Read Group addition, Base Recalibration,
Alignment and sequencing Metrics calculated, Variants Called

### Example Data
Data collected from subsetting a P. falciparum WGS fastq file.