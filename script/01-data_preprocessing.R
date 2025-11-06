library(GEOquery)
library(genefu)
library(biomaRt)
library(tidyverse)
library(data.table)

# Load training cohort (GSE81538) -------------------------------------------------------
exp_set_GSE81538 <- getGEO(filename="data/raw_data/GSE81538/GSE81538_series_matrix.txt")

# Full sample metadata as a data.frame
metadata_GSE81538 <- pData(exp_set_GSE81538)

gene_expression_GSE81538 <- fread("data/raw_data/GSE81538/GSE81538_gene_expression_405_transformed.csv.gz")


# Load validation cohort (GSE96058) -------------------------------------------------------
# subset of patients from the small microarray pilot study
exp_set_GSE96058 <- getGEO(filename="data/raw_data/GSE96058/GSE96058-GPL18573_series_matrix.txt")

# Full sample metadata as a data.frame
metadata_GSE96058_1 <- pData(exp_set_GSE96058)

# Load GSE96058 metadata 
exp_set_GSE96058 <- getGEO(filename="data/raw_data/GSE96058/GSE96058-GPL11154_series_matrix.txt")

# Full sample metadata as a data.frame
metadata_GSE96058_2 <- pData(exp_set_GSE96058)

gene_expression_GSE96058 <- fread("/Volumes/TucciSSD/bioinformatics/public_datasets/Breast_cancer_cohorts/SCAN_B-L_Saal_cohort/data/raw_data/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")

metadata_GSE96058 <- rbind(metadata_GSE96058_1, metadata_GSE96058_2)

# Preprocess expression data ---------------------------------------------------
# 3814 patients = 405 training + 3273 validation + 136 technical replicates (from validation cohort)
table(str_detect(colnames(gene_expression_GSE96058), "repl"))

# log2(fpkm + 0.1) 
log_fpkm_mat <- full_join(gene_expression_GSE81538, gene_expression_GSE96058, by = "V1") %>%
  # remove the 136 technical replicates
  dplyr::select(-contains("repl")) %>% 
  column_to_rownames("V1") %>%
  as.matrix()

# validatate that the minimum value is log2(0.1)
near(min(log_fpkm_mat, na.rm = T), log2(0.1), tol = 1e-12)

# replace NA
log_fpkm_mat[is.na(log_fpkm_mat)] <- log2(0.1)

# convert fpkm to tpm according to Bullard et al. 2010 (10.1186/1471-2105-11-94)
logFPKM2TPM <- function(log2fpkm) {
  # revert log2(FPKM + 1)
  fpkm <- 2^log2fpkm - 0.1
  tpm <- fpkm / sum(fpkm) * 1e6
  log2(tpm + 1)
}

log_tpm_mat <- log_fpkm_mat %>%
  apply(2, logFPKM2TPM)
 
fwrite(log_fpkm_mat, "data/processed/log_FPKM.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE, row.names = TRUE)
fwrite(log_tpm_mat, "data/processed/log_TPM.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE, row.names = TRUE)

# define SCMOD2 class -----------------------------------------------------------
# load symbol to entrez mapping 
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# entrez2symbol <- getBM(mart = mart,
#                        attributes = c("external_gene_name", "entrezgene_id"))
# 
# # cache biomart
# fwrite(entrez2symbol, "data/cache/entrez2symbol.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
entrez2symbol <- fread("data/cache/entrez2symbol.tsv")

table(rownames(log_fpkm_mat) %in% c(entrez2symbol$external_gene_name))
entrez2symbol_vector <- set_names(entrez2symbol$entrezgene_id, entrez2symbol$external_gene_name)


annot <- data.frame(
  probe        = rownames(log_fpkm_mat),    # required but can be any unique id
  Gene.Symbol  = rownames(log_fpkm_mat),    # <-- key column for mapping
  EntrezGene.ID  = entrez2symbol_vector[rownames(log_fpkm_mat)],    # <-- key column for mapping
  stringsAsFactors = FALSE,
  row.names = rownames(log_fpkm_mat)
) 

data(scmod2.robust)

scmod2_classification <- molecular.subtyping(
  sbt.model = "scmod2",          # the pretrained model
  data      = t(log_fpkm_mat),     # your expression matrix
  annot     = annot,         
  do.mapping = T
)

scmod2_class <- scmod2_classification[["subtype"]] %>% 
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID") %>%
  rename("SCMOD2" =  ".")

# parse metadata ---------------------------------------------------------------
metadata_GSE96058_clean <- metadata_GSE96058 %>%
  select(title, geo_accession, ends_with(":ch1"), -contains("reading")) %>%
  rename_with(~str_remove(.x, ":ch1$")) %>%
  rename_with(~str_replace_all(.x, " ", "_")) %>%
  # replace "NA" string with NA
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>%
  # remove all leading and trailing spaces
  mutate(across(where(is.character), ~ str_trim(.x))) %>%
  transmute(
    PATIENT_ID = title,
    Cohort = "validation",
    GEO_accession = geo_accession,
    External_id = `scan-b_external_id`,
    Age = as.numeric(age_at_diagnosis),
    Age_10_years_increase = Age/10,
    ER_IHC = case_match(er_status,
                        "0" ~	"Negative",
                        "1" ~	"Positive",
    ) %>% factor(levels = c("Negative", "Low", "Positive")),
    PR_IHC = case_match(pgr_status,
                        "0" ~	"Negative",
                        "1" ~	"Positive",
    ) %>% factor(levels = c("Negative", "Low", "Positive")),
    HER2 = case_when(
      her2_status == "0" ~ "Normal",
      her2_status == "1" ~ "Amplified"
    ),
    Ki67 = case_when(
      ki67_status == "1" ~ "> 20%",
      ki67_status == "0" ~ "0-20%"
    ),
    Nottingham_grade = nhg,
    Pam50 = pam50_subtype,
    Chemotherapy = case_match(chemo_treated,
                              "0" ~ "No",
                              "1" ~ "Yes",
    ),
    Hormone_therapy  = case_match(endocrine_treated,
                                  "0" ~ "No",
                                  "1" ~ "Yes",
    ),
    pT = case_when(
      tumor_size == 0 ~ "Tis",
      # not using pT staging, but rather Tumor_size because no infomration for pT4 are available
      tumor_size > 0 & tumor_size <= 20 ~ "≤ 2 cm",
      tumor_size > 20 & tumor_size <= 50 ~ "> 2 - 5 cm",
      tumor_size > 50 ~ "> 5 cm",
      TRUE ~ "TX") %>%
      factor(levels = c("≤ 2 cm",  "> 2 - 5 cm", "> 5 cm", "Tis", "TX")),
    pN = case_match(lymph_node_group,
                    "NodeNegative" ~ "N0",
                    "SubMicroMet" ~ "N0(i+)",
                    "1to3" ~ "N1",
                    "4toX" ~ "N2-3",
                    .default = "NX"
    ) %>%
      factor(levels = c("N0", "N0(i+)", "N1", "N2-3", "NX")),
    OS = as.numeric(overall_survival_event),
    OS_years = as.numeric(overall_survival_days)/365.25,
    Instrument_model = instrument_model
  )

# second group ---------------------------------------------
metadata_GSE81538_clean <- metadata_GSE81538 %>%
  select(title, geo_accession, ends_with(":ch1"), -contains("reading")) %>%
  rename_with(~str_remove(.x, ":ch1$")) %>%
  rename_with(~str_replace_all(.x, " ", "_")) %>%
  # replace "NA" string with NA
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>%
  # remove all leading and trailing spaces
  mutate(across(where(is.character), ~ str_trim(.x))) %>%
  transmute(
    PATIENT_ID = title,
    Cohort = "Training",
    GEO_accession = geo_accession,
    # Age = as.numeric(age_at_diagnosis),
    # Age_10_years_increase = Age/10,
    ER_IHC = case_match(er_consensus,
                        "0" ~	"Negative",
                        "1" ~	"Low",
                        "2"	~ "Positive",
                        "3"	~ "Positive",
    ) %>% factor(levels = c("Negative", "Low", "Positive")),
    PR_IHC = case_match(pgr_consensus,
                        "0" ~	"Negative",
                        "1" ~	"Low",
                        "2"	~ "Positive",
                        "3"	~ "Positive",
    ) %>% factor(levels = c("Negative", "Low", "Positive")),
    HER2 = case_when(
      her2_consensus == "0" ~ "Normal",
      her2_consensus == "1" ~ "Amplified"
    ),
    Ki67_pct = as.numeric(ki67_consensus),
    Ki67 = case_when(
      as.numeric(ki67_consensus) > 20 ~ "> 20%",
      as.numeric(ki67_consensus) <= 20 ~ "0-20%"
    ),
    Nottingham_grade = paste0("G", nhg_consensus),
    Pam50 = pam50_subtype,
    instrument_model = "HiSeq 2000"
  )

metadata <- bind_rows(metadata_GSE81538_clean, metadata_GSE96058_clean) %>%
  left_join(scmod2_class) %>%
  relocate(SCMOD2, .after = Pam50)

# export ----------------
write_tsv(metadata, "data/processed/SCANB-metadata.tsv")
saveRDS(metadata, "data/processed/SCANB-metadata.rds")
