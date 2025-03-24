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

source("code/helper/my_toolbox.R")
source("code/helper/ggstyles.R")
source("~/.R/ggstyles.R")
theme_set(gtheme(14))

Sys.setlocale(locale = "fr_FR.UTF-8")
set.seed(1234)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sce = readRDS("data/sce_merged.RDS")

sample_clonal_hierarchy = data.frame(
  sample = c("M180", "M182", "M183", "M184", "M185", "M186", "M187"),
  ppm1d_clonal_hierarchy = c("subclonal", "co-/dominant", "bystander", "co-/dominant", "co-/dominant", "bystander", "co-/dominant")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 6a: UMAP patient GR021
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

final_detected = read_excel("data/single_cell_data.xlsx", sheet=1)
s_m180 = import_sample_norm_prot("M180")

m180_ngt_umap = DimPlot(s_m180,reduction = "UMAP", pt.size=.8, group.by = c("PPM1D.510", "TET2.I274I",  "DNMT3A.R635W", "RUNX1.D198N")) +
  plot_layout(guides="collect", nrow=2) & theme(legend.position="bottom", legend.text = element_text(size=16), strip.text.x = element_text(size=18), axis.text = element_text(size=14), title = element_text(size=21)) & scale_ngt

m180_pheno_umap = FeaturePlot_scCustom(
  s_m180,
  reduction = "UMAP",
  features = c("CD34", "CD38", "HLA-DR", "CD117", "CD7", "CD13"),
  pt.size = .1, na_cutoff = 1, max.cutoff = "q99",
  colors_use = rev(MetBrewer::met.brewer("Hokusai1",n=100)),
  raster = F,
  num_columns = 3
) &
  theme(
    legend.position="right",
    legend.text = element_text(size=16),
    strip.text.x = element_text(size=18),
    axis.text = element_text(size=14),
    title = element_text(size=21)
  )

fig_6a_pl = (m180_ngt_umap | m180_pheno_umap) + plot_layout(widths=c(1,1.9)) + plot_annotation(title="GR021") & dimtheme

ggsave2(
  "figures/main/figure_6/fig_6a.png",
  fig_6a_pl,
  height=105, width=210, scale=2, unit="mm", dpi=400, bg="white"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 6b: Differential protein abundance according to PPM1D status
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ppm1d_var = rowData(sce) %>%
  as.data.frame() %>%
  filter(SYMBOL == "PPM1D") %>% 
  rownames_to_column("ID") %>% 
  pull(ID)

ppm1d_ngt = t(assay(sce, "NGT")[ppm1d_var,]) %>%
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

ppm1d_mut = ppm1d_ngt %>% 
  mutate(PPM1D_BIN = ifelse(Het != "" | Hom != "", 1, 0)) %>% 
  pull(PPM1D_BIN)

ppm1d_mapping = colData(sce) %>%
  as.data.frame() %>%
  mutate(PPM1D_BIN = ppm1d_mut) %>%
  rownames_to_column("Barcode")

processed_df =  t(assay(altExp(sce))) %>%
  data.frame() %>%
  rownames_to_column("Barcode") %>%
  left_join(ppm1d_mapping, by="Barcode") %>%
  left_join(sample_clonal_hierarchy, by="sample") %>%
  relocate(sample, PPM1D_BIN, ppm1d_clonal_hierarchy, .before="CD10") %>%
  pivot_longer(c(5:23)) %>%
  filter(!name %in% c("IgG2a", "IgG2b", "CD22")) %>%
  mutate(PPM1D_BIN = factor(ifelse(PPM1D_BIN==1, "MUT", "WT"), levels = c("WT", "MUT") ))

counts = processed_df %>%
  group_by(ppm1d_clonal_hierarchy, name, PPM1D_BIN) %>%
  summarize(mean_expression = mean(value)) %>%
  pivot_wider(names_from = 'PPM1D_BIN', values_from = 'mean_expression') %>% 
  ungroup() %>%
  mutate(log2FoldChange = log2(MUT/WT)) %>%
  mutate(zscore = scale(log2FoldChange))

p_values = processed_df %>%
  group_by(ppm1d_clonal_hierarchy, name ) %>%
  rstatix::wilcox_test(value ~ PPM1D_BIN)

clonal_hierarchy_plot = left_join(counts, p_values) %>%
  mutate(padj = p.adjust(p), ) %>% 
  mutate(ppm1d_clonal_hierarchy = factor(paste0("PPM1D\n", ppm1d_clonal_hierarchy), levels = c("PPM1D\nbystander", "PPM1D\nsubclonal", "PPM1D\nco-/dominant")))

dotplot = ggplot(clonal_hierarchy_plot, aes(x=name, y=ppm1d_clonal_hierarchy, color=log2FoldChange)) +
  geom_point(aes(shape = ifelse(padj >= 0.05, NA, "s")), size = 3.5, stroke = 2, color = "grey10" ) +
  geom_point(shape=16, size=4) +
  scico::scale_color_scico(
    palette = "vik", midpoint = 0,
    begin = .1, end = .9, breaks = scales::breaks_pretty(3),
    limits = c(-max(abs(clonal_hierarchy_plot$log2FoldChange)), max(abs(clonal_hierarchy_plot$log2FoldChange)) )
  ) +
  scale_shape_manual(values = c(21), name = "FDR < 0.05", labels = c("yes", "")) +
  mytheme(8) +
  theme(
    panel.grid.major = element_line(color="grey80", linetype=2, linewidth = .3),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    plot.title = element_text(face = "bold", size=rel(1.25))
  ) +
  guides(
    color = guide_colorbar(
      title = "Log2FC (Mut vs. WT)", order=1,
      tticks.linewidth= .75/.pt, frame.linewidth = .5/.pt, title.vjust=.56, title.position = "top",
      barwidth = unit(.35, 'lines'), barheight= unit(4, 'lines')
    )
  ) +
  ggtitle("Single cell differential expression of surface proteins")

ggsave2(
  "figures/main/figure_6/fig_6b.png",
  dotplot,
  width=120, height=35, dpi=400, bg="white", unit="mm", scale=1.6
)

