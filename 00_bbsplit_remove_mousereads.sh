#!/bin/bash
#SBATCH -J bbsplit
#SBATCH -o log/bbsplit_%j.out  # Output file with the job ID
#SBATCH -e log/bbsplit_%j.err  # Error file with the job ID
#SBATCH -t 8:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH --qos=normal
#SBATCH --nodes=1  # ask for number of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 #threads
#SBATCH --mem=64GB  # Specify memory allocation/RAM

set -euo pipefail

module load bbtools

# Paths to human and mouse igenomes reference files
HUMAN_REF=/path/to/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
MOUSE_REF=/path/to/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa

# Loop through all Read 1 files
for r1 in *_R1_001.fastq.gz; do
    # Define Read 2
    r2="${r1/_L002_R1_001.fastq.gz/_L002_R2_001.fastq.gz}"
    # Base output name
    base_name="${r1%_L002_R1_001.fastq.gz}"
    echo "Processing sample: $base_name"
    bbsplit.sh \
        in1="$r1" \
        in2="$r2" \
        ref="$HUMAN_REF,$MOUSE_REF" \
        out1="${base_name}_human_R1.fastq.gz" \
        out2="${base_name}_human_R2.fastq.gz" \
        ambig2=toss \
        overwrite=true
done

