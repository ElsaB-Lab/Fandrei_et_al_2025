.cran_packages = c("tidyverse", "scales", "psych", "cowplot", "ggsci",
                   "MetBrewer", "ggstar", "ggpubr", "ggh4x", "ggnewscale",
                   "patchwork", "readxl")

## Loading library
for (pack in .cran_packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

# source("code/helper/ggstyles.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# UMI-seq
ww = read_xlsx("data/ww.xlsx")

# Clinical data
d = read.delim("../../../data/ccohort/new_base.tsv", dec = ",")[c(1:41, 142:155)]
dd <- d[which(d$NIP %in% nip.seq),]








