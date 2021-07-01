cd $(dirname "$0")/..

SAMPLE_ID=$1
IN_DIR=$2
OUT_DIR=$3

if [ "$SAMPLE_ID" = "lung_1_1" ]; then

  cp "$IN_DIR/1cycle_GFP_MAX_LungB2 20x closed PH Tumor1.1 zoom1.oib - C=1.tif" "$OUT_DIR/lung_1_1-cycle-01-channel-01.tif"
  cp "$IN_DIR/1cycle_RFP_MAX_LungB2 20x closed PH Tumor1.1 zoom1.oib - C=2.tif" "$OUT_DIR/lung_1_1-cycle-01-channel-02.tif"
  cp "$IN_DIR/1cycle_Cy5_MAX_LungB2 20x closed PH Tumor1.1 zoom1.oib - C=3.tif" "$OUT_DIR/lung_1_1-cycle-01-channel-03.tif"
  
  cp "$IN_DIR/2cycle_DAPI_MAX_LungB2 cycle 2 20x closed PH Tumor 1.1 zoom1.oib - C=0.tif" "$OUT_DIR/lung_1_1-cycle-02-channel-01.tif"
  cp "$IN_DIR/2cycle_GFP_MAX_LungB2 cycle 2 20x closed PH Tumor 1.1 zoom1.oib - C=1.tif" "$OUT_DIR/lung_1_1-cycle-02-channel-02.tif"
  cp "$IN_DIR/2cycle_RFP_MAX_LungB2 cycle 2 20x closed PH Tumor 1.1 zoom1.oib - C=2.tif" "$OUT_DIR/lung_1_1-cycle-02-channel-03.tif"
  cp "$IN_DIR/2cycle_Cy5_MAX_LungB2 cycle 2 20x closed PH Tumor 1.1 zoom1.oib - C=3.tif" "$OUT_DIR/lung_1_1-cycle-02-channel-04.tif"
  
  cp "$IN_DIR/3cycle_DAPI_MAX_LungB2 20x closed PH 3cycle Tumor1.1 zoom1.oib - C=0.tif" "$OUT_DIR/lung_1_1-cycle-03-channel-01.tif"
  cp "$IN_DIR/3cycle_GFP_MAX_LungB2 20x closed PH 3cycle Tumor1.1 zoom1.oib - C=1.tif" "$OUT_DIR/lung_1_1-cycle-03-channel-02.tif"
  cp "$IN_DIR/3cycle_RFP_MAX_LungB2 20x closed PH 3cycle Tumor1.1 zoom1.oib - C=2.tif" "$OUT_DIR/lung_1_1-cycle-03-channel-03.tif"
  cp "$IN_DIR/3cycle_Cy5_MAX_LungB2 20x closed PH 3cycle Tumor1.1 zoom1.oib - C=3.tif" "$OUT_DIR/lung_1_1-cycle-03-channel-04.tif"

elif [ "$SAMPLE_ID" = "lung_2_1" ]; then

  cp "$IN_DIR/1cycle_DAPI_MAX_LungB2 20x closed PH Tumor2.1 zoom1.oib - C=0.tif" "$OUT_DIR/lung_2_1-cycle-01-channel-01.tif"
  cp "$IN_DIR/1cycle_GFP_MAX_LungB2 20x closed PH Tumor2.1 zoom1.oib - C=1.tif" "$OUT_DIR/lung_2_1-cycle-01-channel-02.tif"
  cp "$IN_DIR/1cycle_RFP_MAX_LungB2 20x closed PH Tumor2.1 zoom1.oib - C=2.tif" "$OUT_DIR/lung_2_1-cycle-01-channel-03.tif"
  cp "$IN_DIR/1cycle_cy5_MAX_LungB2 20x closed PH Tumor2.1 zoom1.oib - C=3.tif" "$OUT_DIR/lung_2_1-cycle-01-channel-04.tif"
  
  cp "$IN_DIR/2cycle_DAPI_MAX_LungB2 cycle 2 20x closed PH Tumor 2.1 zoom1.oib - C=0.tif" "$OUT_DIR/lung_2_1-cycle-02-channel-01.tif"
  cp "$IN_DIR/2cycle_GFP_MAX_LungB2 cycle 2 20x closed PH Tumor 2.1 zoom1.oib - C=1.tif" "$OUT_DIR/lung_2_1-cycle-02-channel-02.tif"
  cp "$IN_DIR/2cycle_RFP_MAX_LungB2 cycle 2 20x closed PH Tumor 2.1 zoom1.oib - C=2.tif" "$OUT_DIR/lung_2_1-cycle-02-channel-03.tif"
  cp "$IN_DIR/2cycle_Cy5_MAX_LungB2 cycle 2 20x closed PH Tumor 2.1 zoom1.oib - C=3.tif" "$OUT_DIR/lung_2_1-cycle-02-channel-04.tif"
  
  cp "$IN_DIR/3cycle_DAPI_MAX_LungB2 20x closed PH 3cycle Tumor2.1 zoom1.oib - C=0.tif" "$OUT_DIR/lung_2_1-cycle-03-channel-01.tif"
  cp "$IN_DIR/3cycle_GFP_MAX_LungB2 20x closed PH 3cycle Tumor2.1 zoom1.oib - C=1.tif" "$OUT_DIR/lung_2_1-cycle-03-channel-02.tif"
  cp "$IN_DIR/3cycle_RFP_MAX_LungB2 20x closed PH 3cycle Tumor2.1 zoom1.oib - C=2.tif" "$OUT_DIR/lung_2_1-cycle-03-channel-03.tif"
  cp "$IN_DIR/3cycle_Cy5_MAX_LungB2 20x closed PH 3cycle Tumor2.1 zoom1.oib - C=3.tif" "$OUT_DIR/lung_2_1-cycle-03-channel-04.tif"

elif [ "$SAMPLE_ID" = "lung_2_2" ]; then

  cp "$IN_DIR/1cycle_DAPI_MAX_LungB2 20x closed PH Tumor2.2 zoom1.oib - C=0.tif" "$OUT_DIR/lung_2_2-cycle-01-channel-01.tif"
  cp "$IN_DIR/1cycle_GFP_MAX_LungB2 20x closed PH Tumor2.2 zoom1.oib - C=1.tif" "$OUT_DIR/lung_2_2-cycle-01-channel-02.tif"
  cp "$IN_DIR/1cycle_RFP_MAX_LungB2 20x closed PH Tumor2.2 zoom1.oib - C=2.tif" "$OUT_DIR/lung_2_2-cycle-01-channel-03.tif"
  cp "$IN_DIR/1cycle_Cy5_MAX_LungB2 20x closed PH Tumor2.2 zoom1.oib - C=3.tif" "$OUT_DIR/lung_2_2-cycle-01-channel-04.tif"
  
  cp "$IN_DIR/2cycle_DAPI_MAX_LungB2 cycle 2 20x closed PH Tumor 2.2 zoom1.oib - C=0.tif" "$OUT_DIR/lung_2_2-cycle-02-channel-01.tif"
  cp "$IN_DIR/2cycle_GFP_MAX_LungB2 cycle 2 20x closed PH Tumor 2.2 zoom1.oib - C=1.tif" "$OUT_DIR/lung_2_2-cycle-02-channel-02.tif"
  cp "$IN_DIR/2cycle_RFP_MAX_LungB2 cycle 2 20x closed PH Tumor 2.2 zoom1.oib - C=2.tif" "$OUT_DIR/lung_2_2-cycle-02-channel-03.tif"
  cp "$IN_DIR/2cycle_Cy5_MAX_LungB2 cycle 2 20x closed PH Tumor 2.2 zoom1.oib - C=3.tif" "$OUT_DIR/lung_2_2-cycle-02-channel-04.tif"
  
  cp "$IN_DIR/3cycle_DAPI_MAX_LungB2 20x closed PH 3cycle Tumor2.2 zoom1.oib - C=0.tif" "$OUT_DIR/lung_2_2-cycle-03-channel-01.tif"
  cp "$IN_DIR/3cycle_GFP_MAX_LungB2 20x closed PH 3cycle Tumor2.2 zoom1.oib - C=1.tif" "$OUT_DIR/lung_2_2-cycle-03-channel-02.tif"
  cp "$IN_DIR/3cycle_RFP_MAX_LungB2 20x closed PH 3cycle Tumor2.2 zoom1.oib - C=2.tif" "$OUT_DIR/lung_2_2-cycle-03-channel-03.tif"
  cp "$IN_DIR/3cycle_Cy5_MAX_LungB2 20x closed PH 3cycle Tumor2.2 zoom1.oib - C=3.tif" "$OUT_DIR/lung_2_2-cycle-03-channel-04.tif"

else

  echo "Sample not found."

fi
