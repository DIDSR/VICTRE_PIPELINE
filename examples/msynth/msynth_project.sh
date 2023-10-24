SEED=1
RESULTS_DIRECTORY="./results/"
BREAST_DENSITY='dense'
DOSE=2.04e+09
LESION_DENSITY=1.1
LESION_FILE="/route/to/spiculated/mass_2_size5.0.h5"
FLATFIELD_FILE="/route/to/pregenerated/${BREAST_DENSITY}/flatfield.raw"
PHANTOM_FILE="/route/to/pregenerated/phantoms/${BREAST_DENSITY}/pc_${SEED}_crop.raw.gz"


python msynth_VICTRE.py --results $RESULTS_DIRECTORY \
                        --seed $SEED \
                        --density $BREAST_DENSITY \
                        --dose $DOSE \
                        --simulate_lesions \
                        --lesion_file $LESION_FILE \
                        --lesion_density $LESION_DENSITY \
                        --phantom_file $PHANTOM_FILE \
                        --flatfield_file $FLATFIELD_FILE

# add --lesion_only flag for a lesion only/segmentation projection