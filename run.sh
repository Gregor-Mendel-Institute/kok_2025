#!/bin/bash

clean=TRUE

# stage files if they are not already staged
if [ ! -d "/scratch-cbe/users/elin.axelsson/IsoSeq_Results" ]; then
    tar -xzf /groups/berger/lab/Raw/20250313_skku_YoonHwanSu_IsoSeq_1.tar.gz \
        -C /scratch-cbe/users/elin.axelsson \

else
    echo "IsoSeq results already staged."
fi

if [ clean ] ; then
    echo "Cleaning up old output files..."
    rm -rf output
    rm -rf modified
    echo "Cleanup complete."
else
    echo "Skipping cleanup of old output files."
fi

# always remove out files first
rm *.out
rm *.err

scripts=("scripts/minimap2_trans.sbatch" 
    "scripts/tama_collaps.sbatch" 
    "scripts/gffcompare.sh" 
    "scripts/bed_to_gff.sh" 
    "scripts/bed2gffandbeyond.sh" 
    "scripts/kallisto.sbatch"
    "scripts/as_bedtogtf.sh" 
    "scripts/compare_with_antisense.sh"
    "scripts/R.sbatch")



job_id=""


# Submit jobs in order, with dependencie

for script in "${scripts[@]}"; do
    if [ -z "$job_id" ]; then
        # First job
        echo "Submitting first job: $script"
        job_id=$(sbatch --parsable "$script")
        echo "Submitted $script with job ID: $job_id"
    else
        # Dependent jobs
        echo "Submitting dependent job: $script with dependency on job ID: $job_id"
        job_id=$(sbatch --parsable --dependency=afterok:$job_id "$script")
        echo "Submitted $script with job ID: $job_id"
    fi
    echo "Submitted $script with job ID: $job_id"
done