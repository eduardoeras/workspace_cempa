#!/bin/bash

#Usage:
#./MONAN_Bias.sh 2025120100 120
#where 120 is the forecast lead time in hours

INIT_DATE=$1
MAX_LEAD=$2
OUTPUT_PATH=$3

echo "Generating bias figures for: "$INIT_DATE" "$MAX_LEAD

INIT_DATE_EPOCH=$(date -u -d "${INIT_DATE:0:4}-${INIT_DATE:4:2}-${INIT_DATE:6:2} ${INIT_DATE:8:2}:00:00" +%s)

BASE_MONAN="${OUTPUT_PATH}/MONAN/${INIT_DATE:0:4}${INIT_DATE:4:2}/${INIT_DATE}"
BASE_BIAS="${OUTPUT_PATH}/Bias/${INIT_DATE:0:4}${INIT_DATE:4:2}/${INIT_DATE}"

WORK_DIR=$(pwd)
OUTPUT_DIR="${OUTPUT_PATH}/Bias/${INIT_DATE:0:4}${INIT_DATE:4:2}/${INIT_DATE}/"

DOMAINS=("GLB" "AMS" "ACC")

SIZE="800x600"

for dom in "${DOMAINS[@]}"; do

    for lead in $(seq 24 24 ${MAX_LEAD}); do

        LEAD3=$(printf "%03d" ${lead})

        END_DATE=$(date -u -d "@$((INIT_DATE_EPOCH + lead*3600))" +%Y%m%d%H)

        echo "Processing lead time ${LEAD3}h"

        # Input files
        #IMG_MONAN="${BASE_MONAN}/MONAN_24precacum_${INIT_DATE}_${END_DATE}_${dom}.png"
        IMG_GPM="${BASE_BIAS}/diff_MONAN_GPM_${dom}_${LEAD3}h.png"
        IMG_GSMAP=${BASE_BIAS}/"diff_MONAN_GSMAP_${dom}_${LEAD3}h.png"
        IMG_MSWEP=${BASE_BIAS}/"diff_MONAN_MSWEP_${dom}_${LEAD3}h.png"

        # Temp files
        #TMP1="tmp_monan_${LEAD3}.png"
        TMP2="tmp_gpm_${LEAD3}.png"
        TMP3="tmp_gsmap_${LEAD3}.png"
        TMP4="tmp_mswep_${LEAD3}.png"

        # Output file
        OUT="Bias_MONAN_${dom}_${LEAD3}.png"

        # Skip if MONAN image doesn't exist
        #if [ ! -f "$IMG_MONAN" ]; then
        #    echo "Skipping montage for ${LEAD3}h ${dom}: MONAN image not found at $IMG_MONAN"
        #    continue
        #fi

        # Standardize 
        #convert ${IMG_MONAN} \
        #-trim \
        #-bordercolor white \
        #-border 50x50 \
        #-background white \
        #-gravity center \
        #+repage \
        #"${TMP1}"

        convert ${IMG_GPM} \
        -trim \
        -bordercolor white \
        -border 50x50 \
        -background white \
        -gravity center \
        +repage \
        "${TMP2}"

        convert ${IMG_GSMAP} \
        -trim \
        -bordercolor white \
        -border 50x50 \
        -background white \
        -gravity center \
        +repage \
        "${TMP3}"

        convert ${IMG_MSWEP} \
        -trim \
        -bordercolor white \
        -border 50x50 \
        -background white \
        -gravity center \
        +repage \
        "${TMP4}"

        # 2x2 Montage
        #montage "${TMP1}" "${TMP2}" "${TMP3}" "${TMP4}" \
        montage "${TMP2}" "${TMP3}" "${TMP4}" \
        -tile 2x2 \
        -geometry +1+1 \
        -background white \
        "${OUTPUT_DIR}/${OUT}"

        # Remove temp files
        #rm -f ${TMP1} ${TMP2} ${TMP3} ${TMP4}
        rm -f ${TMP2} ${TMP3} ${TMP4}
        rm -f ${IMG_GPM} ${IMG_GSMAP} ${IMG_MSWEP}

    done
done