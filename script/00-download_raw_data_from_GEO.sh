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


# [GSE81538] SCAN-B training cohort: 405 patients ------------------------------------------------------------------
mkdir -p "${OUT_DIR}/GSE81538"
cd "${OUT_DIR}/GSE81538"

wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_UCSC_Human_hg19_knownGenes_GTF_appended_10sep2012.gtf.gz"

# download the gene expression file: log2(FPKM + 0.1)
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_gene_expression_405_transformed.csv.gz"

# download the transcript expression file: raw FPKM
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_map_transcriptID_geneSymbol.csv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_transcript_expression_405.csv.gz"

# download pathology consensus key file
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/suppl/GSE81538_pathology_consensus_key.xlsx"

# download the soft file
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/soft/GSE81538_family.soft.gz"
gunzip "GSE81538_family.soft.gz"

# download the family file
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/miniml/GSE81538_family.xml.tgz"
# remove after the file
tar -xzf "GSE81538_family.xml.tgz"

# download the series matrix files
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81538/matrix/GSE81538_series_matrix.txt.gz"
gunzip "GSE81538_series_matrix.txt.gz"

echo "Download of raw data from GSE81538 completed."

# [GSE96058] SCAN-B validation cohort: 3273 samples (136 have technical replicates) -------------------------------
mkdir -p "${OUT_DIR}/GSE96058"
cd "${OUT_DIR}/GSE96058"

# download the gene expression file: log2(FPKM + 0.1)
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"

# download the transcript expression file: raw FPKM
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_transcript_expression_3273_samples_and_136_replicates.csv.gz"

# download the soft file
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/soft/GSE96058_family.soft.gz"
gunzip "GSE96058_family.soft.gz"

# download the family file
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/miniml/GSE96058_family.xml.tgz"
tar -xzf "GSE96058_family.xml.tgz"

# download the series matrix files
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL11154_series_matrix.txt.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL18573_series_matrix.txt.gz"
gunzip "GSE96058*series_matrix.txt.gz"

# echo "Download of raw data from GSE96058 completed."


