#!/bin/bash

#===================================================================
# SBATCH DIRECTIVES
#===================================================================
# Job name:
#SBATCH --job-name=APS_sim_array
#
# Slurm logs:
#SBATCH --output=slurm_logs/sim-%A_%a.out
#SBATCH --error=slurm_logs/sim-%A_%a.err
#
# Array setup:
# 7 bins * N realizations = 7*N tasks
#SBATCH --array=1-2000%30
#
# Resource allocation:
#SBATCH --time=100:30:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

#===================================================================
# SCRIPT 
#===================================================================

echo "------------------------------------------------"
echo "SLURM JOB ID:     $SLURM_JOB_ID"
echo "SLURM ARRAY TASK: $SLURM_ARRAY_TASK_ID"
echo "------------------------------------------------"

eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/PY10

DIR_SCRIPTS='/data/AMARINS/LSSxHI-CODES/scripts'
VERBOSE=False 
FIELD1="HI"
FIELD2="DENSITY"
BASE_PROJECT='LSSTxLOWZ_40_0001' 
PATHOUT="/data/AMARINS/LSSxHI-DATA/simulations/GCxHI"
SUFFIX="LOWZ"

NREALIZATIONS_TOTAL=2000  
SEED_BASE=9000           
LIMITED_CORRELATED_CHANNELS='False'
CHANNEL_TAX=1010
CHANNEL_MIN_CORR=1
CHANNEL_MAX_CORR=30
jTOMO=40
SIGMAZ="0001"
RUN_FILE="$DIR_SCRIPTS/generating_APS_simulations_parallel.py"

#####
# Map 1 - M*N to Bins (1-M) and Realizations (0 - N-1)
# i_BIN: 1 to 7 [galaxy bins]
i_BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) / $NREALIZATIONS_TOTAL + 1 ))
# i_REAL: 0 to N-1 (This acts as our offset)
i_REAL=$(( ($SLURM_ARRAY_TASK_ID - 1) % $NREALIZATIONS_TOTAL ))
# CURRENT_SEED: Calculates the specific seed for this single job
CURRENT_SEED=$(( $SEED_BASE + $i_REAL ))
echo "Processing Bin: $i_BIN | Realization Offset: $i_REAL | Seed: $CURRENT_SEED"

#####
# PATHs
PROJECT_NAME="${BASE_PROJECT}_bin${i_BIN}"
FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i_BIN}_${jTOMO}_${SIGMAZ}_${FIELD1}_cl_${SUFFIX}.txt"
FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i_BIN}_${jTOMO}_${SIGMAZ}_${FIELD2}_cl_${SUFFIX}.txt"
FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i_BIN}_${jTOMO}_${SIGMAZ}_${FIELD1}x${FIELD2}_cl_${SUFFIX}.txt"

####
# OPTIMIZATION 
# Force libraries to respect the single-core allocation.
# This prevents "thrashing" where multiple threads fight for the single CPU.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- 5. Execute Python ---
# We pass --nrealizations 1 because this job handles exactly one realization.
# We pass --sim_index_offset so the folder is named 'sim_15' instead of always 'sim_0'.
echo "Starting python script..."
python3 $RUN_FILE \
    --verbose $VERBOSE \
    --project $PROJECT_NAME \
    --pathout $PATHOUT \
    --filepath_field1 $FILEPATH_FIELD1 \
    --filepath_field2 $FILEPATH_FIELD2 \
    --filepath_cross $FILEPATH_CROSS \
    --limited_correlated_channels $LIMITED_CORRELATED_CHANNELS \
    --channel_min_corr $CHANNEL_MIN_CORR \
    --channel_max_corr $CHANNEL_MAX_CORR \
    --channel_tax $CHANNEL_TAX \
    --nrealizations 1 \
    --seed0 $CURRENT_SEED \
    --sim_index_offset $i_REAL

echo "Job Finished"