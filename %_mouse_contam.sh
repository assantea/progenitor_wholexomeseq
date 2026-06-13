#!/bin/bash
#SBATCH --partition amilan
#SBATCH -t 02:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH --qos=normal
#SBATCH --nodes=1  # ask for number of nodes
#SBATCH --ntasks=1
#SBATCH --mem=8GB  # Specify memory allocation/RAM


# Define output file
OUTPUT="summary.tsv"

# Write the header to the file (overwrites existing file)
echo -e "sample\traw file\thuman only\tdiscarded reads\t% discarded" > "$OUTPUT"

# Loop through human-only R1 files
for human_file in *_human_R1.fastq.gz; do
    
    # 1. Identify the matching raw R1 file
    # Pattern: swaps '_human_R1' for raw suffix '_L002_R1_001'
    raw_file=$(echo "$human_file" | sed 's/_human_R1/_L002_R1_001/')
    
    # 2. Extract sample name (e.g., 220754_RB6_LowSR)
    # This removes 'p01_', the '_S##' part, and the file extension
    sample_base=$(echo "$human_file" | sed -E 's/^p01_//; s/_S[0-9]+_human_R1.fastq.gz//')

    # 3. Process if raw file exists
    if [ -f "$raw_file" ]; then
        echo "Processing $sample_base..."

        # Count reads (total lines / 4)
        raw_count=$(zcat "$raw_file" | wc -l | awk '{print $1/4}')
        clean_count=$(zcat "$human_file" | wc -l | awk '{print $1/4}')

        # 4. Perform calculations
        discarded=$((raw_count - clean_count))
        percent=$(awk -v r="$raw_count" -v d="$discarded" 'BEGIN {printf "%.8f", (d/r)*100}')

        # 5. Append results to the summary.tsv
        echo -e "${sample_base}_R1\t$raw_count\t$clean_count\t$discarded\t$percent" >> "$OUTPUT"
    else
        echo "Warning: Could not find raw file for $human_file"
    fi
done

echo "Done"
