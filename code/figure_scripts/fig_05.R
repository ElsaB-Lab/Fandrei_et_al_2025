.cran_packages = c("scDNA", "SingleCellExperiment", "Seurat", "readxl")


## Loading library
for (pack in .cran_packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

source("helper/my_toolbox.R")
source("helper/ggstyles.R")
theme_set(gtheme(14))

Sys.setlocale(locale = "fr_FR.UTF-8")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

h5path = '/Users/davidfandrei/Desktop/PPM1D/tapestri/data/merged.h5'

dsb_norm_prot <- read.csv("../../analysis/norm_prot_matrix.csv") %>%
  column_to_rownames("X") %>% 
  dplyr::select(-IgG2a, -IgG2b) %>%
  t() %>% 
  as.data.frame()

# Clin data
dd = read_excel("data/dd.xlsx")

# Molecular data
maf <- read_xlsx("data/maf_final.tsv")

# detected variants
final_detected = readRDS("../../analysis/final_detected.RDS")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function sce to Seurat
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sce_to_seurat = function(sce, adt_norm) {
  
  sce<-scDNA::enumerate_clones(sce)
  
  exclude = which(rownames(altExp(sce, "Protein")) %in% "CD22")
  
  altExp(sce, "Protein") = altExp(sce, "Protein")[-exclude,]
  
  assay(sce, "DSB_norm") <- adt_norm
  
  cdata_prot = droplet_metadata %>% dplyr::filter(Cell %in% colnames(protein_mat)) %>% column_to_rownames("Cell") %>% S4Vectors::DataFrame()
  protein_sce@colData <- cdata_prot
  
  SingleCellExperiment::altExp(sce, "Protein") <- protein_sce
  
  # altExp(sce, "Protein") = altExp(sce, "Protein")[which(!rownames(altExp(sce, "Protein")) %in% c("IgG2a", "IgG2b")),]
  
  s<-Seurat::as.Seurat(x = altExp(sce),counts="Protein",data="DSB_norm")
  s<-SeuratObject::RenameAssays(object = s, originalexp = 'Protein')
  DefaultAssay(s) <- "Protein"
  
  return(s)
  
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Get variants
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

include_variants = final_detected %>%
  group_by(id) %>%
  filter(row_number() == 1, single_cell_id != "M181")

sce_merged <-tapestri_h5_to_sce(
  file=h5path,
  variant_set = include_variants,
  GT_cutoff = 0,
  DP_cutoff = 0,
  GQ_cutoff = 0,
  protein = F
)

# Compare ADT levels between groups
ppm1d_var = rowData(sce_merged) %>%
  as.data.frame() %>%
  filter(SYMBOL == "PPM1D") %>% 
  rownames_to_column("ID") %>% 
  pull(ID)

ppm1d_ngt = t(assay(sce_merged, "NGT")[ppm1d_var,]) %>%
  data.frame() %>%
  mutate(Het = "", Hom = "")

for (var in colnames(ppm1d_ngt)[1:6]) {
  for (cell in rownames(ppm1d_ngt)) {
    if (ppm1d_ngt[cell, var] == 1) {
      ppm1d_ngt[cell,"Het"] <- paste0(ppm1d_ngt[cell,]$Het, var, ";") 
    } else if (ppm1d_ngt[cell, var] == 2) {
      ppm1d_ngt[cell,"Hom"] <- paste0(ppm1d_ngt[cell,]$Hom, var, ";") 
    }
  }
}