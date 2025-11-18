#library(GEOquery)
library(genefu)
library(biomaRt)
library(tidyverse)
library(data.table)
library(readxl)

days_to_year <- function(days) {
  return(days / 365.25)
}

# load NPJ data ----------------------------------------------------------------
metadata <- read_excel("data/raw_data/npj_2022/Supplemental Data Table/Supplementary Data Table 1 - 2023-01-13.xlsx",
                           sheet = "SCANB.9206", guess_max = 5000) %>%
  mutate(
    Age = `Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)`,
    Age_10_years_increase = Age/10,
    Tumor_size_mm = Size.mm,
    pN = case_match(LN.spec,
                         "N0" ~ "N0",
                         "SubMicroMet" ~ "N0(i+)",
                         "1to3" ~ "N1",
                         "4toX" ~ "N2-3",
                         .default = "NX"
        ) %>% factor(levels = c("N0", "N0(i+)", "N1", "N2-3", "NX")),
    Nottingham_grade = paste0("G", NHG) %>% ifelse(. == "GNA", NA, .)
    )  %>%
  # convert days in years
  mutate(across(ends_with("days"), days_to_year)) %>%
  rename_with(~str_replace(.x, "days", "years"), ends_with("days")) %>%
  rename_with(~str_to_upper(.x), ends_with("event")) %>%
  rename_with(~str_to_upper(.x), ends_with("years")) %>%
  # keep only the 6660 patients that passed all the quality check and skip technical duplicates
  filter(Follow.up.cohort == T) %>%
  select(Patient, Case, Sample, GEX.assay,
         Age = `Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)`,
         Age_10_years_increase,
         Histotype = InvCa.type,
         ER, ER_pct = ER.pct, PR, HER2, Ki67,
         Nottingham_grade, pT, Tumor_size_mm, pN,
         treatment_group = TreatGroup, Clinical_group = ClinGroup,
         # inferred labels
         SSP_Pam50 = SSP.PAM50, SSP_Subtype = SSP.Subtype,
         SSP_Nottingham_grade = SSP.NHG, SSP_Ki67 = SSP.Ki67,
         ends_with("years"), ends_with("event")) %>%
  rename_all(~str_remove(.x, "_event$"))
  
# check no duplicated patients
any(metadata_NPJ$Patient %>% table() > 1)

exp_data <- fread("/Volumes/TucciSSD/bioinformatics/public_datasets/Breast_cancer_cohorts/SCAN_B-analysis_boilerplate/data/raw_data/npj_2022/StringTie FPKM Gene Data LibProtocol adjusted/SCANB.9206.genematrix_noNeg.txt") %>%
  column_to_rownames("V1") %>%
  as.matrix()

# keep only the 6660 filtered patients
exp_data <- exp_data[,metadata$GEX.assay]

dim(exp_data)

# log transform FPKM
# fpkm scale requires 0.1 offset 
log_fpkm <- log2(exp_data + 0.1)
fwrite(log_fpkm, "data/processed/log_FPKM.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE, row.names = TRUE)



# perform SCMOD2 imputation ----------
gene_anno <- fread("data/raw_data/npj_2022/Gene Annotations/Gene.ID.ann.txt")

ensembl2symbol <- set_names(gene_anno$HGNC, gene_anno$Gene.ID, )
ensembl2entrez <- set_names(gene_anno$EntrezGene, gene_anno$Gene.ID,)

annot <- data.frame(
  probe        = rownames(log_fpkm),    # required but can be any unique id
  Gene.Symbol  = ensembl2symbol[rownames(log_fpkm)],    # <-- key column for mapping
  EntrezGene.ID  = ensembl2entrez[rownames(log_fpkm)],    # <-- key column for mapping
  stringsAsFactors = FALSE,
  row.names = rownames(log_fpkm)
) 

data(scmod2.robust)

scmod2_classification <- molecular.subtyping(
  sbt.model = "scmod2",          # the pretrained model
  data      = t(log_fpkm),     # your expression matrix
  annot     = annot,         
  do.mapping = T
)

scmod2_class <- scmod2_classification[["subtype"]] %>% 
  as.data.frame() %>%
  rownames_to_column("GEX.assay") %>%
  rename("SCMOD2" =  ".") %>%
  mutate(SCMOD2 = factor(SCMOD2, levels = c("ER+/HER2- Low Prolif", "ER+/HER2- High Prolif",  "HER2+", "ER-/HER2-")))

data(pam50.robust)

pam50_classification <- molecular.subtyping(
  sbt.model = "pam50",          # the pretrained model
  data      = t(log_fpkm),     # your expression matrix
  annot     = annot,         
  do.mapping = T
)

pam50_class <- pam50_classification[["subtype"]] %>% 
  as.data.frame() %>%
  rownames_to_column("GEX.assay") %>%
  rename("Pam50" =  ".") %>%
  mutate(Pam50 = factor(Pam50, levels = c("LumA", "LumB",  "Her2", "Basal", "Normal")))

# merge SCMOD2 class into metadata
metadata.1 <- metadata %>%
  left_join(scmod2_class, by = "GEX.assay") %>%
  left_join(pam50_class, by = "GEX.assay")

# export ----------------
write_tsv(metadata.1, "data/processed/SCANB-metadata.tsv")
saveRDS(metadata.1, "data/processed/SCANB-metadata.rds")



