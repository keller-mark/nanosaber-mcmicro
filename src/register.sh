cd $(dirname "$0")/..

IJDIR=/Applications/Fiji.app
SAMPLE_ID=$1
  
# Register all cycle 1 images onto cycle 2 DAPI channel

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-01.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-01-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-01-registered.tif"

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-02.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-01-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-02-registered.tif"

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-03.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-01-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-03-registered.tif"

if [ "$SAMPLE_ID" != "lung_1_1" ]; then
  java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
    "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
    "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-04.tif" \
    "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-01-channel-01_direct_transf.txt" \
    "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-01-channel-04-registered.tif"
fi
  

# Register all cycle 3 images onto cycle 2 DAPI channel

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-01.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-03-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-01-registered.tif"

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-02.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-03-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-02-registered.tif"

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-03.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-03-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-03-registered.tif"

java -Xmx512m -cp $IJDIR/jars/ij-1.53c.jar:$IJDIR/plugins/bUnwarpJ_-2.6.13.jar bunwarpj.bUnwarpJ_ -elastic_transform \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-02-channel-01.tif" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-04.tif" \
  "data/$SAMPLE_ID/bunwarpj/$SAMPLE_ID-cycle-03-channel-01_direct_transf.txt" \
  "data/$SAMPLE_ID/pre-raw/$SAMPLE_ID-cycle-03-channel-04-registered.tif"