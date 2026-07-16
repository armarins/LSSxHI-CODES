#!/bin/bash
NAME_RUN="BINGO_mask1"
FILENAME_NEWDIR="mask1"

#######
DIR_OUTSCREEN='/data/AMARINS/CMBWLxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/CMBWLxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT

#######
#TERMINAL INFO
PROJECT='bingo_nch30_980_1260'
DIRPATH_SIMS='/data/AMARINS/CMBWLxHI-DATA/simulations/bingo_nch30_980_1260'
FILEPATH_FIELD1='/data/AMARINS/CMBWLxHI-CODES/theoretical/bingo_HI_cl_nch30_980_1260.txt'
FILEPATH_FIELD2='/data/AMARINS/CMBWLxHI-CODES/theoretical/bingo_CMBWL_cl_nch30_980_1260.txt'
FILEPATH_CROSS='/data/AMARINS/CMBWLxHI-CODES/theoretical/bingo_CMBWLxHI_cl_nch30_980_1260.txt'
DIRPATH_SAVEDATA='/data/AMARINS/CMBWLxHI-DATA/FGremoval'
VERBOSE=1
NSIMS=100

############
#MASK INFO
APPLY_MASK=1
DIRPATH_MASK='/data/AMARINS/CMBWLxHI-DATA/MAPS/MASK'
FILENAME_MASK='mask_256_apodcutless_overlap_allra_mask1.fits'
#filename_mask='mask_256_apodcutless_overlap_allra_mask1.fits'
#filename_mask='mask_256_50fgcut_apodcutless_overlap_allra_mask2.fits'
#filename_mask='mask_256_apodcutless_overlap_allra_mask3.fits'
#filename_mask='mask_256_70fgcut_apodcutless_overlap_allra_mask4.fits'
#filename_mask='mask_256_apodcutless_overlap_allra_mask5.fits'
##############
RUN_NAME='generating_APS_FGremoval'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

#############
#source /home/fba/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate amarins

for i in 2 3 4 5
do
    python $PATHFILE --verbose $VERBOSE --project $PROJECT --dirpath_sims $DIRPATH_SIMS --filename_newdir $FILENAME_NEWDIR\
                     --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                     --dirpath_savedata $DIRPATH_SAVEDATA --apply_mask $APPLY_MASK --dirpath_mask $DIRPATH_MASK --filename_mask $FILENAME_MASK\
                     --ns $i --nrealizations $NSIMS | tee -a $FILEOUT

done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT

