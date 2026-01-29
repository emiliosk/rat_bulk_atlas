#!/bin/bash
#SBATCH -A naiss2025-22-779
#SBATCH -p main
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH -J align_to_ref
#SBATCH -e logs/align_to_ref_%j.err
#SBATCH -o logs/align_to_ref_%j.out
#SBATCH --mail-user=emilio.skarwan@scilifelab.se
#SBATCH --mail-type=ALL

# Load modules
module load PDC/24.11
module load kallisto/0.51.1-cpeCray-24.11

# --- Configuration ---
WORK_DIR=/cfs/klemming/home/s/skarwan/EmilioTemp/rat/rat_bulk
output_dir=$WORK_DIR/local/kallisto_out
metadata="$WORK_DIR/metadata/filename_metadata_rat.tsv"
reference_file=$WORK_DIR/local/reference/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.index
data_dir=/cfs/klemming/projects/snic/rnaatlas/private

export WORK_DIR output_dir reference_file data_dir


process_sample_group() {
    local sample_name="$1"
    shift 
    local raw_files=("$@")

    if [ -f "$output_dir/$sample_name/abundance.h5" ]; then
        echo "   [SKIPPING] Results already exist for $sample_name"
        return
    fi
    
    if [[ -z "$sample_name" ]]; then return; fi
    echo "Processing Sample: $sample_name"
    echo "   Found ${#raw_files[@]} file(s) for this sample."

    local full_paths=()
    for f in "${raw_files[@]}"; do
        local f_path="$data_dir/$f"
        
        if [ -f "$f_path" ]; then
            echo "  [FOUND]   $f"
            full_paths+=("$f_path")
        else
            echo "   [ERROR] File missing: $f_path"
        fi
    done

    if [ ${#full_paths[@]} -gt 0 ]; then
        mkdir -p "$output_dir/$sample_name"
        
        # Note: Kallisto expects paired files in order (fwd rev fwd rev). 
        # Metadata is ordered (1 then 2), so this should work.
        kallisto quant -i "$reference_file" -t 8 -o "$output_dir/$sample_name" "${full_paths[@]}"
    else
        echo "   [SKIP] No valid files found for $sample_name"
    fi
}

export -f process_sample_group

# --- Main Execution ---

# 1. Parse Metadata with AWK
#    NR>1      : Skip the header row
#    $(NF-1)   : The second to last column (Sample Name)
#    $1        : The first column (File Path)
#    Logic     : Store files in a map: map[sample_name] = "file1 file2 file3..."
#    END       : Print "sample_name file1 file2..." for every sample

awk 'NR>1 { 
    sample=$(NF-1); 
    file=$1; 
    
    if (sample in map) {
        map[sample] = map[sample] " " file
    } else {
        map[sample] = file
    }
} 
END { 
    for (s in map) {
        print s, map[s]
    }
}' "$metadata" | \
# 2. Pass to xargs
#    -L 1 : Take 1 line at a time (One sample + all its files)
#    -P 16: Run 16 samples in parallel
xargs -L 1 -P 16 bash -c 'process_sample_group "$@"' _


