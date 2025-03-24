.cran_packages = c("tidyverse", "readxl", "ggsci", "RColorBrewer", "ggpubr", 
                   "tibble", "prodlim", "kableExtra", "ggforce", "gridExtra",
                   "patchwork", "table1", "survival", "survminer", "rstatix")

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
theme_set(gtheme(14))

Sys.setlocale(locale = "fr_FR.UTF-8")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Clin data
dd = read_excel("data/dd.xlsx") %>%
  mutate(DIAG2 = factor(DIAG2, levels = c("CH", "CCUS", "MDS", "AML", "MPN")))

# Keep PPM1D at baseline
fsamples <- dd %>% dplyr::filter(seq_0none_1first_2others_Xnoppm1d %in% c(0,1)) %>% pull(GRId_sample)
ddb = dd[dd$GRId_sample %in% fsamples, ]

# Mol Data
ddmut <- read.table("data/ddmut.tsv", header=T, sep="\t", stringsAsFactors=F)
ffmut <- read.table("data/ffmut.tsv", header=T, sep="\t", stringsAsFactors=F)
maf <- read_xlsx("data/maf_final.xlsx")

# Subset molecular data at baseline
ddmutb <- ddmut[ddmut$GRId_sample %in% fsamples , ]
ffmutb <- ffmut[ffmut$GRId_sample %in% fsamples , ]
mafb <- maf[maf$GRId_sample %in% fsamples , ]

# Add some clinical information to maf
mafb <- left_join(mafb, ddb[,c("GRId_sample", "DIAG2","Sex", "t21")], by="GRId_sample")

# Order genes by frequency
mafb$Gene <- factor(mafb$Gene, levels=rev(names(sort(table(mafb$Gene)))))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Inference of PPM1D clonal hierarchy
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Correct VAF for known cytogenetic aberration and sex
mafb <- mafb %>% 
  mutate(VAF_corrected = case_when(
    Sex == "M" & Chr == "chrX" ~ VAF / 2,
    Gene %in% c("TP53", "JAK2", "EZH2") & VAF > 60 ~ VAF / 2,
    Gene == "RUNX1" & t21 == T ~ VAF * .75,
    Gene == "U2AF2" & VAF > 70 ~ VAF * .75,
    Sample_Gene_Pchange == "GR081_1_RUNX1_L148Q" ~ VAF / 2,
    .default = VAF
  ))

# Function to get table with mutations per sample that allow decision of clonal 
# hierarchy from clinical panel NGS sequencing at baseline for each sample
estimate_clonal_hierarchy <- function(sample, maf) {
  
  # Filter sample
  sample_df <- maf[maf$GRId == sample,]
  
  # Check if there are comutations
  comutation_df <- sample_df %>%
    filter(Gene != "PPM1D")
  
  # Get PPM1D mutation and co-mutation with highest VAF
  ppm1d_max_vaf <- sample_df %>%
    dplyr::filter(Gene == "PPM1D") %>%
    dplyr::filter(VAF_corrected == max(VAF_corrected))
  
  # Save
  hierarchy_df = setNames(
    data.frame(matrix(ncol = 13, nrow = 1), row.names=NULL),
    c("GRId", "diagnosis", "PPM1D_Pchange", "PPM1D_VAF", "largest_comut_gene",
      "largest_comut_pchange", "largest_comut_VAF", "largest_comut_VAF_corrected",
      "p_value", "OR",	"CI95",	"q_value", "ppm1d_clonal_hierarchy")
  )
  
  hierarchy_df = hierarchy_df %>%
    mutate(GRId = sample, diagnosis = sample_df$DIAG2[1], PPM1D_Pchange = ppm1d_max_vaf$Pchange,
           PPM1D_VAF = ppm1d_max_vaf$VAF)
  
  # If there are not any co-mutations, classify as 'only_ppm1d'
  if (nrow(comutation_df) == 0) {
    
    hierarchy_df$ppm1d_clonal_hierarchy <- "Only PPM1D"
    
    # If there are co-mutations, perform Fisher's exact test if differences
    # in VAF are significant
  } else {
    
    comutation_max_vaf <- comutation_df %>%
      dplyr::filter(VAF_corrected == max(VAF_corrected)) %>%
      dplyr::filter(row_number() == 1)
    
    hierarchy_df = hierarchy_df %>%
      mutate(largest_comut_gene = comutation_max_vaf$Gene, largest_comut_pchange = comutation_max_vaf$Pchange,
             largest_comut_VAF = comutation_max_vaf$VAF, largest_comut_VAF_corrected = comutation_max_vaf$VAF_corrected)
    
    # If sum of VAF < 35% or PPM1D mutation VAF < 3%, PPM1D mutation is considered
    # a bystander mutation
    if (ppm1d_max_vaf$VAF_corrected + comutation_max_vaf$VAF_corrected < 35 | ppm1d_max_vaf$VAF_corrected < 5) {
      
      hierarchy_df$ppm1d_clonal_hierarchy <- "Bystander"
      
    } else {
      
      # Perform Fisher's exact test
      mat <- data.frame(
        "alt_alleles" = c((ppm1d_max_vaf$VAF_corrected/100) * ppm1d_max_vaf$Depth, (comutation_max_vaf$VAF_corrected/100) * comutation_max_vaf$Depth),
        "wt_alleles" = c(ppm1d_max_vaf$Depth - ((ppm1d_max_vaf$VAF_corrected/100) * ppm1d_max_vaf$Depth), comutation_max_vaf$Depth - (comutation_max_vaf$Depth*(comutation_max_vaf$VAF_corrected/100))),
        row.names = c(ppm1d_max_vaf$Gene, comutation_max_vaf$Gene)
      )
      sig <- fisher.test(mat)
      
      hierarchy_df = hierarchy_df %>%
        mutate(p_value =  signif(sig$p.value, 2), OR = round(sig$estimate, 2),
               CI95 = paste0("[", round(sig$conf.int[1], 2), "-", round(sig$conf.int[2], 2), "]"))
    }
  }
  
  return(hierarchy_df)
}

# Iterate over samples
samples = ddb$GRId
hierarchy_table = do.call(
  rbind,
  lapply(samples, estimate_clonal_hierarchy, maf = mafb)
)

# p value adjustment
hierarchy_table$q_value = p.adjust(hierarchy_table$p_value,  method = "bonferroni")

# Make decisions for subclonal, dominant, co-dominant cases
hierarchy_table = hierarchy_table %>%
  mutate(ppm1d_clonal_hierarchy = case_when(
    q_value >= .05 ~ "Co-dominant",
    q_value < .05 & PPM1D_VAF > largest_comut_VAF ~ "Dominant",
    q_value < .05 & largest_comut_VAF > PPM1D_VAF ~ "Subclonal",
    .default = ppm1d_clonal_hierarchy
  ))

table(hierarchy_table$ppm1d_clonal_hierarchy, useNA = "ifany")

# Result --> supplementary table
write.csv(
  hierarchy_table,
  "figures/supplement/table_1.csv",
  row.names = F,
  quote=F
)

# Add to clinical and molecular data
mafb = left_join(mafb, hierarchy_table[,c("GRId", "ppm1d_clonal_hierarchy")])
ddb = left_join(ddb, hierarchy_table[,c("GRId", "ppm1d_clonal_hierarchy")], by="GRId")

# proportions amongst diagnostic groups
freq_diag <- data.frame(table(ddb$DIAG2)) %>% mutate(labels = paste0(Var1, "\n", "n= ", Freq))
freq_list <- freq_diag$labels
freq_list

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Characterization of groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cl_boxplot = ddb %>%
  mutate(ppm1d_clonal_hierarchy = factor(
    recode(ppm1d_clonal_hierarchy, Bystander="Bystander/uncertain \n"),
    levels = c("Co-dominant", "Dominant", "Subclonal", "Bystander/uncertain \n", "Only PPM1D")
  )) %>%
  group_by(ppm1d_clonal_hierarchy) %>%
  summarize(n=n()) %>% 
  mutate(perc=n/sum(n)*100) %>%
  ggplot(aes(x="", y=perc, fill = ppm1d_clonal_hierarchy)) +
  geom_bar(position="stack", stat="identity", alpha=.75, color="black") +
  guides(fill="none") +
  geom_text(aes(label=paste0(ppm1d_clonal_hierarchy, " n=", n, " (", round(perc, 0), "%)")), position = position_stack(vjust = 0.5), color="black", size=2.6) +
  scale_fill_manual(values = unlist(unname(anno.clonality.all[-3]))) + 
  coord_trans(ylim=c(0,100)) +
  theme_classic() + gtheme(12) +
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x = element_blank()) +
  ylab("Groups in %")
cl_boxplot

ggsave2(
  "figures/main/figure_4/fig4b.png",
  cl_boxplot,
  width = 45, height = 60, dpi = 500, bg = "white", units = "mm", scale = 1.3
)

# Regroup dominant and co-dominant
ddb = ddb %>%
  mutate(
    ppm1d_clonal_hierarchy = factor(
      ifelse(ppm1d_clonal_hierarchy=="Dominant" | ppm1d_clonal_hierarchy == "Co-dominant", "Co-/dominant", ppm1d_clonal_hierarchy),
      levels = c("Co-/dominant", "Subclonal", "Bystander", "Only PPM1D")
    ))

mafb = mafb %>%
  mutate(
    ppm1d_clonal_hierarchy = factor(
      ifelse(ppm1d_clonal_hierarchy=="Dominant" | ppm1d_clonal_hierarchy == "Co-dominant", "Co-/dominant", ppm1d_clonal_hierarchy),
      levels = c("Co-/dominant", "Subclonal", "Bystander", "Only PPM1D")
    ))

# PPM1D VAF amongst groups
cl_ppm1d_vaf = mafb %>%
  mutate(ppm1d_clonal_hierarchy = factor(
    recode(ppm1d_clonal_hierarchy, Bystander="Bystander/uncertain"),
    levels = c("Co-/dominant", "Subclonal", "Bystander/uncertain", "Only PPM1D")
  )) %>%
  group_by(GRId) %>%
  filter(Gene == "PPM1D") %>%
  filter(VAF == max(VAF)) %>%
  ggplot(aes(x=ppm1d_clonal_hierarchy, y=VAF_corrected)) + 
  geom_violin(alpha=.8, draw_quantiles = .5, size=.5) + 
  geom_jitter(aes(color=ppm1d_clonal_hierarchy), height = 0, width = 0.2, size=1, alpha=.7) + 
  theme_classic() + gtheme(12) +
  theme(axis.title.x = element_blank()) +
  ylab("PPM1D VAF (%)") +
  angle45 +
  scale_y_continuous(breaks = seq(0, 50, by=10)) +
  scale_color_manual(values = unlist(unname(anno.clonality.all[3:6]))) + 
  stat_compare_means(label = "p.format", label.y = 80, label.x = 1, size = 4) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T) +
  # facet_grid(.~paste0("PPM1D")) +
  guides(fill = "none", color="none")

tapply(subset(mafb, Gene == "PPM1D")$VAF, subset(mafb, Gene == "PPM1D")$ppm1d_clonal_hierarchy, summary)

ggsave2(
  "figures/supplement/ppm1d_vaf_clonal_groups.png",
  cl_ppm1d_vaf,
  width = 55, height = 70, dpi = 500, bg = "white", units = "mm", scale = 1.4
)

# Distribution of groups in diagnostic categories
ddb %>%
  group_by(ppm1d_clonal_hierarchy, DIAG2) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n))

cl_diag_bars = 
  ddb %>%
  mutate(ppm1d_clonal_hierarchy = factor(
    recode(ppm1d_clonal_hierarchy, Bystander="Bystander/uncertain"),
    levels = c("Co-/dominant", "Subclonal", "Bystander/uncertain", "Only PPM1D")
  )) %>%
  ggplot(aes(x=DIAG2, fill=ppm1d_clonal_hierarchy)) +
  geom_bar(color="black", position="fill", alpha=.8) +
  theme_classic() + gtheme(12) + 
  xlab("Diagnosis") + ylab("Proportion of patients") + topleg +
  scale_fill_manual(name = "Group", values = unlist(unname(anno.clonality.all[3:6]))) +
  coord_flip() +
  scale_x_discrete(labels = freq_list) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(vjust=8))

# Co-mutations barplot
# Get 15 most commonly mutated genes and their status for each patient
patient_gene_table = mafb %>%
  mutate(GRId = factor(GRId)) %>%
  group_by(GRId, Gene) %>%
  summarize(n=n()) %>%
  mutate(n = ifelse(n>1, 1, n)) %>%
  complete(Gene = c(levels(mafb$Gene)[2:length(levels(mafb$Gene))]), fill = list(n = 0)) %>%
  left_join(ddb[,c("GRId", "ppm1d_clonal_hierarchy")]) %>%
  filter(!ppm1d_clonal_hierarchy %in% "Only PPM1D") %>%
  mutate(ppm1d_clonal_hierarchy = factor(ppm1d_clonal_hierarchy,levels = c("Co-/dominant", "Subclonal", "Bystander"))) %>%
  ungroup()

binary_cl = ddb %>%
  filter(DIAG2 != "MPN") %>%
  mutate(
    DIAG3 = ifelse(DIAG2 %in% c("CH", "CCUS"), "CH/CCUS", "MDS/AML" ),
    codom_binary = ifelse(ppm1d_clonal_hierarchy == "Co-/dominant", 1, 0),
    codom_subclonal = ifelse(ppm1d_clonal_hierarchy %in% c("Co-/dominant", "Subclonal"), 1, 0)
  )

fisher.test(binary_cl$DIAG3, binary_cl$codom_binary)
fisher.test(binary_cl$DIAG3, binary_cl$codom_subclonal)

ftable(patient_gene_table$Gene, patient_gene_table$n, patient_gene_table$ppm1d_clonal_hierarchy)

# proportions amongst hierarchy groups
hier_diag <- data.frame(table(ddb$ppm1d_clonal_hierarchy))

get_p_value = function(x) {
  gene_table = patient_gene_table[patient_gene_table$Gene == x,]
  cont_table = table(gene_table$n, gene_table$ppm1d_clonal_hierarchy)
  test = rstatix::fisher_test(cont_table)
  return(data.frame(Gene = x, p_value = test$p))
}

signif_table = do.call(
  rbind,
  lapply(levels(mafb$Gene)[2:10], get_p_value)
)
signif_table$q_value = p.adjust(signif_table$p_value)
signif_table$p_signif = pval_to_signif_code(signif_table$q_value)

prop_df = ddmutb %>%
  pivot_longer(4:ncol(ddmutb)) %>%
  left_join(ddb[,c("GRId", "ppm1d_clonal_hierarchy")]) %>%
  mutate(name = factor(name, levels = levels(mafb$Gene))) %>%
  group_by(name, ppm1d_clonal_hierarchy) %>%
  summarize(n = sum(value)) %>%
  mutate(prop = case_when(
    ppm1d_clonal_hierarchy == "Co-/dominant" ~ n / hier_diag[hier_diag$Var1 == "Co-/dominant",]$Freq,
    ppm1d_clonal_hierarchy == "Subclonal" ~ n / hier_diag[hier_diag$Var1 == "Subclonal",]$Freq,
    ppm1d_clonal_hierarchy == "Bystander" ~ n / hier_diag[hier_diag$Var1 == "Bystander",]$Freq,
    ppm1d_clonal_hierarchy == "Only PPM1D" ~ n / hier_diag[hier_diag$Var1 == "Only PPM1D",]$Freq
  ))

comut_bars = prop_df %>%
  mutate(ppm1d_clonal_hierarchy = factor(
    recode(ppm1d_clonal_hierarchy, Bystander="Bystander/uncertain"),
    levels = c("Co-/dominant", "Subclonal", "Bystander/uncertain", "Only PPM1D")
  )) %>%
  filter(name %in% levels(mafb$Gene)[2:10]) %>%
  ggplot(aes(y=100*prop, x=name, fill=ppm1d_clonal_hierarchy)) +
  geom_bar( position=position_dodge(preserve="single"), stat="identity", color="black", alpha=.8, show.legend = F) +
  geom_text(data = signif_table, aes(x=Gene, y = 52, label = p_signif), size=6.5, inherit.aes = F) +
  theme_classic() + gtheme(12) +
  theme(axis.title.x = element_blank(), legend.position = c(0.83, 0.7), axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  ylab("% of patients per group") +
  scale_fill_manual(name="PPM1D clonal hierarchy", values=unlist(unname(anno.clonality.all[3:6])))

clonal_char = cl_diag_bars + comut_bars + plot_layout(guides="collect", widths=c(1, 1.25)) & theme(legend.position = "bottom")

ggsave2(
  "figures/main/figure_4/fig4c_d.png",
  clonal_char,
  width = 180, height = 90, dpi = 500, bg = "white", units = "mm", scale = 1
)

# Description co-mutations
# TP53
patients_tp53 = mafb %>%
  filter(Gene == "TP53") %>%
  group_by(GRId) %>%
  filter(row_number() == 1) %>%
  pull(GRId)

tp53_binary = ddb %>%
  mutate(tp53_status = ifelse(GRId %in% patients_tp53, 1, 0), ppm1d_clonal_hierarchy2 = ifelse(ppm1d_clonal_hierarchy=="Bystander", "bystander", "other"))

tp53_binary$ppm1d_clonal_hierarchy2 = factor(tp53_binary$ppm1d_clonal_hierarchy2, levels=c("other", "bystander"))

tp53_binary %>%
  group_by(ppm1d_clonal_hierarchy2, tp53_status) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n))

fisher.test(tp53_binary$ppm1d_clonal_hierarchy2, tp53_binary$tp53_status)

# DNMT3A
patients_dnmt3a = mafb %>%
  filter(Gene == "DNMT3A") %>%
  group_by(GRId) %>%
  filter(row_number() == 1) %>%
  pull(GRId)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Survival according to TP53 status
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

surv_df = ddb %>%
  filter(!is.na(Date_LFU)) %>%
  left_join(ddmutb[,c("GRId", "TP53")]) %>%
  mutate(
    os = ifelse(DIAG2 == "MPN",
                interval(Date_NGS, Date_LFU) %/% months(1),
                interval(Date_diagnosis, Date_LFU) %/% months(1)
    ),
    Death = ifelse(Death == "oui", 1 , 0),
    TP53 = ifelse(TP53 == 1, "mut", "WT")
  )

tp53_mds = survdiff(Surv(os, Death) ~ TP53, data = subset(surv_df, DIAG2 == "MDS"))
tp53_aml = survdiff(Surv(os, Death) ~ TP53, data = subset(surv_df, DIAG2 == "AML"))

surv <- as.formula(paste0("Surv(", "os", ",", "Death", ") ~", "TP53 + DIAG2"))
surv_fit <- survfit(surv, data=subset(surv_df, DIAG2 %in% c("MDS", "AML")))
surv_fit$call$formula <- surv

g_surv <- survminer::ggsurvplot(
  surv_fit,
  data=subset(surv_df, DIAG2 %in% c("MDS", "AML")), 
  color= "DIAG2", 
  legend.title="Diagnosis",
  linetype = "TP53",
  size=.75,
  conf.int= F,
  xlim=c(0,33),
  break.time.by=6, 
  pval=F,
  surv.median.line="none",
  risk.table=TRUE,
  # tables.col="TP53",
  # tables.y.text=FALSE,
  risk.table.title="No. at risk",
  risk.table.fontsize=6,
  ggtheme = theme_classic(),
)

cols <- c(MDS = anno_colour$Diagnosis[[3]], AML = anno_colour$Diagnosis[[4]])

g_pl = g_surv$plot +
  annotate("text", x=28, y = 0.98, label = paste0("log-rank test: p=", signif(tp53_mds$pvalue, 1)), size=4.5, color = anno_colour$Diagnosis[[3]]) +
  annotate("text", x=28, y = 0.9, label = paste0("log-rank test: p=", signif(tp53_aml$pvalue, 1)), size=4.5, color = anno_colour$Diagnosis[[4]]) +
  scale_colour_manual(values = cols) +
  scale_linetype_manual(values = c( "longdash", "solid")) +
  theme(legend.position = "top") +
  guides(linetype = guide_legend(nrow=2), color = guide_legend(nrow=2)) +
  theme_classic(12) +
  gtheme() +
  theme(legend.key.width = unit(1.1,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_rect(linewidth=1, fill=NA),
        strip.text.x = element_text(size=10),
        legend.position = "top",
        axis.title.y = element_text(vjust=-12.5)
  ) + 
  ylab("OS Probability") +
  xlab("Months after diagnosis")

gt <- g_surv$table +
  theme(
    axis.line=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none",
    plot.title = element_text(size=18),
    axis.text.y = ggtext::element_markdown(color = rep(c(anno_colour$Diagnosis[[4]], anno_colour$Diagnosis[[3]]), 2), size=14)
  ) +
  xlab("") +
  ylab("") +
  scale_y_discrete(labels =rev(c("MDS: TP53 mut", "AML: TP53 mut", "MDS: TP53 WT", "AML: TP53 WT")))

g_pl / gt + plot_layout(ncol = 1, heights = c(3, 1))

ggsave2(
  "figures/main/figure_4/fig_4e.png",
  width=110, height=120, dpi=400, bg="white", units = "mm", scale = 1.5
)

# COX PH
surv_df$ppm1d_clonal_hierarchy = relevel(surv_df$ppm1d_clonal_hierarchy, ref="Only PPM1D")
summary(coxph(Surv(os, Death) ~ ppm1d_clonal_hierarchy, surv_df))

summary(coxph(Surv(os, Death) ~ Age_NGS + ppm1d_clonal_hierarchy + TP53, subset(surv_df, DIAG2 %in% c("CH", "CCUS", "MPN")) ) )
summary(coxph(Surv(os, Death) ~ Age_NGS + ppm1d_clonal_hierarchy + TP53, subset(surv_df, DIAG2 %in% c("MDS", "AML")) ) )
