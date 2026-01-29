#!/bin/bash

# --- Configuration ---
WORK_DIR=/cfs/klemming/home/s/skarwan/EmilioTemp/rat/rat_bulk
output_dir=$WORK_DIR/local/kallisto_out
metadata="$WORK_DIR/metadata/filename_metadata_rat.tsv"
reference_file=$WORK_DIR/local/reference/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz
data_dir=/cfs/klemming/projects/snic/rnaatlas/private

export WORK_DIR output_dir reference_file data_dir

# --- Dry Run Function ---
process_sample_group() {
    local sample_name="$1"
    shift
    local raw_files=("$@")
    
    if [[ -z "$sample_name" ]]; then return; fi

    echo "---------------------------------------------------"
    echo "SAMPLE: $sample_name"

    # 1. Verify Files
    local full_paths=()
    local missing_count=0
    
    for f in "${raw_files[@]}"; do
        local f_path="$data_dir/$f"
        
        if [ -f "$f_path" ]; then
            echo "  [FOUND]   $f"
            full_paths+=("$f_path")
        else
            echo "  [MISSING] $f"  # <--- Watch out for these in the output
            ((missing_count++))
        fi
    done

    # 2. Print the Command (Do not execute)
    if [ ${#full_paths[@]} -gt 0 ]; then
        # Construct the command string for display
        local cmd="kallisto quant -i $reference_file -t 8 -o $output_dir/$sample_name ${full_paths[*]}"
        
        echo ""
        echo "  [COMMAND WOULD BE]:"
        echo "  $cmd"
        echo ""
    else
        echo "  [SKIP] No valid files found. Command would not run."
    fi

    if [ $missing_count -gt 0 ]; then
        echo "  [WARNING] $missing_count files were missing from disk!"
    fi
    echo "---------------------------------------------------"
}

export -f process_sample_group

# --- Execution ---
# changed -P 16 to -P 1 so the output prints nicely in order for you to read
awk 'NR>1 { 
    sample=$(NF-1); 
    file=$1; 
    if (sample in map) { map[sample] = map[sample] " " file } 
    else { map[sample] = file }
} 
END { 
    for (s in map) { print s, map[s] }
}' "$metadata" | \
xargs -L 1 -P 1 bash -c 'process_sample_group "$@"' _

