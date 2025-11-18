#!/usr/bin/env bash
set -ue -o pipefail

OUT_DIR="${PWD}/data/raw_data"

# ask confirmation to create the row data directory under PWD/data/raw_data
echo "This script will create the raw data directory under ${OUT_DIR}."

read -r -p "Do you want to proceed? [y/N]: " choice

case "$choice" in
  [yY]|[yY][eE][sS]) ;;
  *) echo "Aborting."; exit 1 ;;
esac

mkdir -p ${OUT_DIR}

# Donwload updated data from Mendeley Data https://data.mendeley.com/datasets/yzxtxn4nmd/4
# in this releasae they included Relapse-free interval data
cd ${OUT_DIR}

wget https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/yzxtxn4nmd-4.zip
unzip *.zip
rm *.zip
# rename the folder
mv "RNA Sequencing-Based Single Sample Predictors of Molecular Subtype and Risk of Recurrence for Clinical Assessment of Early-Stage Breast Cancer" npj_2022