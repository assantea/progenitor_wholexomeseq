#!/bin/bash
#SBATCH -J sarek
#SBATCH -o log/sarek_%j.out  # Output file with the job ID
#SBATCH -e log/sarek_%j.err  # Error file with the job ID
#SBATCH -t 06:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH --qos=normal
#SBATCH -n 1 -c 2  # Ask for number of nodes/cores
#SBATCH --mem=8GB  # Specify memory allocation

set -eu  

# Make software accessible:
module load nextflow/25.10.2
module load singularity/3.7.4
sep=----------------------------------------
printf -- "$sep\n%s\n%s\n$sep\n%s\n%s\n$sep\n\n" \
        "$(which nextflow)" "$(nextflow -version)" \
        "$(which singularity)" "$(singularity --version)"


## Override the settings in the Nextflow module to organize by project:
export NXF_WORK="/path/to/nf_core_dir/work"
export NXF_TEMP="/path/to/nf_core_dir/tmp"
export NXF_HOME="/path/to/nf_core_dir/nextflow"

## Sample file
samplefile="/path/to/samplefilelist.csv"

## Create output dir for nextflow
mkdir -p "$pipeoutdir"
cd "$pipeoutdir"

ptargetsbed="/path/to/probes_targetregions/S33699751_Padded.bed"

## Make temporary headerless 3-column copy of BED file to play nice with tools:
ln1=$(grep -nPm1 '^chr' "$ptargetsbed" | cut -d : -f 1)
tail -n +"${ln1-1}" "$ptargetsbed" | cut -f 1-3 >"/path/to/ptargets.bed"

sarek381="/path/to/nf-core/nf-core-sarek-3.8.1/3_8_1"

nextflow run "$sarek381" -profile <institute_name_w_preconfig_settings> -ansi-log false \
	--wes \
	--input "$samplefile" \
	--genome GATK.GRCh38 \
	--igenomes_base "/path/to/igenomes/" \
	--fasta "/path/to/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
	--step variant_calling \
	--tools cnvkit,mutect2,snpeff \
    --outdir "$pipeoutdir" \
	--intervals "/path/to/ptargets.bed" \
	-resume

