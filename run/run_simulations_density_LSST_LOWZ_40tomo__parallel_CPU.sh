#!/bin/bash
#####################################################
# ================= CONFIGURATION =================
# Cores to use at once
MAX_JOBS=30 
# Total number of simulations (7 bins * 2000 sims)
TOTAL_TASKS=14000
NREALIZATIONS_TOTAL=2000
SEED_BASE=9000

#####################################################
eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/PY10
#####################################################
# Performance control
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export JOBLIB_MULTIPROCESSING=0
#####################################################
DIR_SCRIPTS='/data/AMARINS/LSSxHI-CODES/scripts'
RUN_FILE="$DIR_SCRIPTS/generating_APS_simulations_parallel.py"
BASE_PROJECT='LSSTxLOWZ_40_0001' 
PATHOUT="/data/AMARINS/LSSxHI-DATA/simulations/GCxHI"
SUFFIX="LOWZ"
FIELD1="HI"
FIELD2="DENSITY"
jTOMO=40
SIGMAZ="0001"

#####################################################
echo "Starting Local Job Manager..."
echo "Total Tasks: $TOTAL_TASKS"
echo "Max Parallel Jobs: $MAX_JOBS"

for (( TASK_ID=1; TASK_ID<=TOTAL_TASKS; TASK_ID++ )); do
    # i_BIN: 1 to 7 #galaxy bins
    i_BIN=$(( ($TASK_ID - 1) / $NREALIZATIONS_TOTAL + 1 ))
    # i_REAL: 0 to N-1 (The offset) #N=number of simulations
    i_REAL=$(( ($TASK_ID - 1) % $NREALIZATIONS_TOTAL ))
    CURRENT_SEED=$(( $SEED_BASE + $i_REAL ))
    PROJECT_NAME="${BASE_PROJECT}_bin${i_BIN}"
    FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i_BIN}_${jTOMO}_${SIGMAZ}_${FIELD1}_cl_${SUFFIX}.txt"
    FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i_BIN}_${jTOMO}_${SIGMAZ}_${FIELD2}_cl_${SUFFIX}.txt"
    FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i_BIN}_${jTOMO}_${SIGMAZ}_${FIELD1}x${FIELD2}_cl_${SUFFIX}.txt"

    # We use '&' to put it in background
    (
        python3 $RUN_FILE \
            --verbose False \
            --project $PROJECT_NAME \
            --pathout $PATHOUT \
            --filepath_field1 $FILEPATH_FIELD1 \
            --filepath_field2 $FILEPATH_FIELD2 \
            --filepath_cross $FILEPATH_CROSS \
            --limited_correlated_channels False \
            --channel_min_corr 1 \
            --channel_max_corr 30 \
            --channel_tax 1010 \
            --nrealizations 1 \
            --seed0 $CURRENT_SEED \
            --sim_index_offset $i_REAL
    ) &

    # Queue Management # Check how many jobs are running
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 1 # If we hit the limit, wait 1 second and check again
    done

    # Print progress every 100 jobs
    if (( $TASK_ID % 100 == 0 )); then
        echo "Launched task $TASK_ID / $TOTAL_TASKS"
    fi

done

echo "All tasks launched. Waiting for final jobs to complete..."
wait
echo "Done!"