#!/bin/bash

DIR_SCRIPTS='/home/amarins/Programmation/cross/scripts'

#######
#TERMINAL INFO
#PROJECT='bingo_nch30_980_1260'
#DIRPATH_SIMS='/data2/cross/simulations/bingo_nch30_980_1260'
#filepath_field1 = /home/amarins/Programmation/cross/scripts/auxiliary_data/bingo_HI_cl_nch30_980_1260.txt
#filepath_field2 = /home/amarins/Programmation/cross/scripts/auxiliary_data/bingo_CMBWL_cl_nch30_980_1260.txt
#filepath_cross  = /home/amarins/Programmation/cross/scripts/auxiliary_data/bingo_CMBWLxHI_cl_nch30_980_1260.txt
VERBOSE=1
NSIMS=5
APPLY_MASK=1
#DIRNAME_MASK='/data1/dataset/MASK'
FILENAME_MASK='mask_256_apodcutless_overlap_allra_mask1.fits'
#filename_mask='mask_256_apodcutless_overlap_allra_mask1.fits'
#filename_mask='mask_256_50fgcut_apodcutless_overlap_allra_mask2.fits'
#filename_mask='mask_256_apodcutless_overlap_allra_mask3.fits'
#filename_mask='mask_256_70fgcut_apodcutless_overlap_allra_mask4.fits'
#filename_mask='mask_256_apodcutless_overlap_allra_mask5.fits'
#########
RUN_NAME='generating_APS_FGremoval'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

source /home/amarins/anaconda3/bin/activate
conda activate myconda

#cd /home/amarins
#$# must get the ns integers as command-line input
#for i in $@
for i in {2..3}
do
#    echo $PATHFILE --verbose $VERBOSE --ns $i --nrealizations $NSIMS
#    python $PATHFILE --verbose $VERBOSE --ns $i --nrealizations $NSIMS
    python $PATHFILE --verbose $VERBOSE --ns $i --nrealizations $NSIMS --apply_mask $APPLY_MASK --filename_mask $FILENAME_MASK 
done
#mkcd

