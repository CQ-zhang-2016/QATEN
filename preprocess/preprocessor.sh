#!/bin/bash

# feature extractor executable
BIN="./get_features"

PDBDIR="$1"
JSONDIR="$2"

# create output folder
mkdir -p $JSONDIR

# loop over all PDBs in the folder
# and extract features
for f in $PDBDIR/*.pdb
do
    ID=`basename $f .pdb`
    echo $ID
    $BIN -i $f -j $JSONDIR/$ID.json -d 10.0
done
