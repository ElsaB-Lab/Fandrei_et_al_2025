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

source("code/helper/ggstyles.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# UMI-seq
ww = read_xlsx("data/ww.xlsx")
samples <- unique(ww$GRId)

# Clinical data
d = read_xlsx("data/dd.xlsx")
fsamples = d %>% filter(seq_0none_1first_2others_Xnoppm1d %in% c(0,1)) %>% pull(GRId_sample)
dd <- d[which(d$GRId %in% samples),]
dd = dd %>% mutate(cohort = case_when(
  grepl("GR", GRId) ~ "GR",
  grepl("ALFA", GRId) ~ "ALFA"
)) %>%
  filter(GRId_sample %in% fsamples) %>%
  mutate(cohort = ifelse(GRId %in% samples, "OvBIOMark", cohort), DIAG2 = factor(DIAG2, levels=c('CH','CCUS','MDS','AML','MPN'))) 

# Mol Data
ddmut <- read.table("data/ddmut.tsv", header=T, sep="\t", stringsAsFactors=F)
ffmut <- read.table("data/ffmut.tsv", header=T, sep="\t", stringsAsFactors=F)
maf <- read_xlsx("data/maf_final.xlsx")

# Import binary treatment information
treatment_bin <- read_excel(
  "data/timepoints_treatment.xlsx"
)
colnames(treatment_bin)[1] = "Sample"
treatment_bin = treatment_bin %>% mutate(Patient = gsub("\\_[0-9]", "", Sample))

ww <- left_join(ww, treatment_bin[, c("Sample", "timepoint_binary")], by="Sample")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Process treatment dataframe
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

treatments_all_path = "data/treatment_groups.xlsx"

get_tt_df <- function(patient) {
  df <- read_excel(treatments_all_path, sheet = patient) %>% tibble()
  df <- df[,!apply(df,2, function(x) all(is.na(x)))]
  df %>% pivot_longer(!Date, names_to="line", values_to="treatment") -> df
  df$line <- car::recode(df$line, " 'Alkylant/Combi'='Pt/Alkylating' ")
  df$GRId <- paste(patient)
  df$treatment[is.na(df$treatment)] <- "None"
  df$treatment <- factor(df$treatment)
  
  return(df)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 3a: Panel PPM1D cohort
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

PlotPatient = function(
    patient,
    col.treatment=c('#749e89','#abccbe','#d9768a','#c399a2','#9f6e71','#41507b',"#e9e9e9",'#7d87b2'),
    col_number = 1,
    verbose = F
) {
  
  # UMI data
  wp <- ww %>% filter(GRId == patient)
  
  vmax <- max(100*wp$VAF)*1.05
  vmin <- min(100*wp$VAF)

  lp = wp %>%
    arrange(cum_days) %>%
    last()

  # Color interpolation for each gene
  mut <- unique(wp$Gene_Pchange)
  glist <- factor(unique(wp$Gene))
  coll = vector("character", 0)
  panel_map <- lapply(glist, function(x) mut[ grep(x, mut) ])
  names(panel_map) <- glist

  for (g in glist) {
    coll <- c(coll, .colorInterpolation(g, mutations=mut, map=panel_map))
  }
  
  # Day 0
  tt = read_excel("data/ww.xlsx", sheet=2)
  ttp = tt %>% filter(GRId == patient)
  
  # ## Panel with cumulative days
  dd.treatment = get_tt_df(patient)
  dd.treatment = dd.treatment %>% mutate(cum_days=as.numeric(difftime(Date, ttp$date_sample), unit="days"))

  # Common x-axis
  common_x_scale = coord_cartesian(xlim = c(0, max(wp$cum_days, dd.treatment$cum_days)))
  
  g41 <- wp %>%
    ggplot(aes(x=cum_days,y=100*VAF,color=Gene_Pchange,group=Gene_Pchange)) +
    geom_vline(data=lp, aes(xintercept=cum_days), color="grey", linetype="dashed", alpha=.8) +
    geom_point(size=3) +
    geom_line(size=1) +
    xlab("Days") + ylab("VAF (%)") +
    ylim(c(0,vmax)) +
    facet_grid(.~paste(patient)) +
    scale_color_manual(values=coll) +
    guides(color=guide_legend(ncol=col_number)) +
    theme_classic() +
    gtheme() +
    theme(
      axis.title.x = element_blank(),
      legend.title=element_blank(),
      legend.position="top",
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_blank(),
      plot.tag = element_text(size=rel(2), face="bold")
    ) +
    common_x_scale
  
  # Treatment plot
  lev = c('Pt/Alkylating','PARPi','VEGFa','Other Chemo')
  
  dd.treatment = dd.treatment %>%
    filter(line %in% lev) %>%
    mutate(
      line = factor(line, levels=lev)
    ) %>%
    dplyr::select(cum_days, line, treatment) %>%
    complete(cum_days, line, fill=list(treatment="None"))
  
  gtr4 = dd.treatment %>%
    ggplot(aes(x=cum_days, y=line, col=treatment, group=line)) +
    geom_line(linewidth=4.25) +
    scale_y_discrete(drop=F) +
    scale_colour_manual(name="", values=col.treat, drop=F) +
    ylab("Treatment") + xlab("Days after inclusion") +
    theme_bw() +
    gtheme() +
    noleg +
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    ) +
    common_x_scale
  
  return(list(g41 = g41, gtr4 = gtr4, dd = dd.treatment))
}

# Plot timepoints PPM1D
res_gr022 = PlotPatient("GR022")
res_gr035 = PlotPatient("GR035")
res_gr067 = PlotPatient("GR067")
res_gr068 = PlotPatient("GR068")
res_gr070 = PlotPatient("GR070")
res_gr051 = PlotPatient("GR051")
res_gr054 = PlotPatient("GR054")
res_gr062 = PlotPatient("GR062")
res_gr063 = PlotPatient("GR063")
res_gr065 = PlotPatient("GR065")

## Compose panel figure 3
format_treat_panel <- theme(axis.text.y = element_blank(), legend.position = "none")
panel = 
  # first row
  (res_gr022$g41 + theme(axis.title.y=element_text(vjust=-20), plot.tag = element_text(face="bold") )) + (res_gr035$g41+noytitle) + (res_gr051$g41+noytitle)  + (res_gr054$g41+noytitle) + (res_gr062$g41+noytitle) +
  
  # treat row
  (res_gr022$gtr4 + noleg) + (res_gr035$gtr4 + format_treat_panel) +  (res_gr051$gtr4 + format_treat_panel) + (res_gr054$gtr4 + format_treat_panel) + (res_gr062$gtr4 + format_treat_panel) +
  
  # second row
  (res_gr063$g41+theme(axis.title.y=element_text(vjust=-20))) + (res_gr065$g41+noytitle) + (res_gr067$g41+noytitle) + (res_gr068$g41+noytitle) + (res_gr070$g41+noytitle) +
  
  # treat row
  (res_gr063$gtr4 + noleg) + (res_gr065$gtr4 + format_treat_panel) + (res_gr067$gtr4 + format_treat_panel) + (res_gr068$gtr4 + format_treat_panel) + (res_gr070$gtr4 + format_treat_panel) +
  
  plot_layout(ncol=5, heights = c(2.8, 1, 2.8, 1)) &
  gtheme(13) +
  theme(legend.text = element_text(size=11), legend.key.size = unit(.35, 'cm'), strip.background = element_rect(fill="white"))

ggsave2(
  "figures/main/figure_3/fig_3a.png",
  panel,
  width = 160, height = 130, dpi = 500, bg = "white", units = "mm", scale = 2.2
)

## Legend
legend_labels = c('Carboplatin','Oxaliplatin',"peg. liposomal doxorubicin",'Bevacizumab', 'Niraparib','Olaparib','Paclitaxel')

# treatment levels
patients = c("GR022", "GR035", "GR051", "GR054", "GR062", "GR063", "GR065", "GR067", "GR068", "GR070")
treat_levels_df = do.call(rbind, lapply(patients, get_tt_df)) %>%
  filter(line %in% c('Pt/Alkylating','PARPi','VEGFa','Other Chemo'), treatment != "None", GRId %in% patients ) %>% 
  rowwise() %>%
  mutate(
    x = rnorm(1),
    y = rnorm(1)
  )
unique(treat_levels_df$treatment)
treat_levels=c("Carbo", "Oxaliplatine", "Caelyx", "Avastin", "Niraparib", "Olaparib", "Taxol")
treat_levels_df$treatment = factor(treat_levels_df$treatment, levels = treat_levels)

treat_scatter = ggplot(treat_levels_df, aes(x = x, y = y, color = treatment)) +
  theme_classic() + gtheme(12) +
  geom_point(size=2) +
  scale_color_manual(name="Treatment", values = col.treat, labels = legend_labels) +
  guides(color = guide_legend(nrow=1, title.position="left", override.aes = list(shape=20, size=5)))
treat_scatter

leg = get_legend(treat_scatter)

ggsave2(
  "figures/main/figure_3/fig_3a_legend.png",
  as_ggplot(leg),
  height=10, width=205, dpi=400, bg="white", units="mm", scale=1.4
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 3b: Mutation growth rate per gene
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fc_mut_groups = www %>%
  dplyr::select(nkey , Patient, Sample, Gene, Gene_Pchange, time_diff = step_days, VAF, timepoint_binary) %>%
  group_by(Patient, Gene_Pchange) %>%
  mutate(lag_VAF = lag(VAF, default = 0)) %>%
  mutate(VAF_diff =  100 * (VAF - lag_VAF), VAF_fc =  VAF / lag_VAF) %>%
  filter(time_diff > 5) %>%
  mutate(
    treatment = factor(ifelse(timepoint_binary == 1, "yes", "no"), levels = c("yes", "no")),
    mut_growth_rate = VAF_diff / time_diff
  )

# Stats on Genes
fc_stats = fc_mut_groups %>%
  group_by(Gene, treatment) %>%
  filter(Gene %in% c("PPM1D", "DNMT3A", "TET2", "TP53")) %>%
  summarize(
    n = dplyr::n(),
    median_mut_growth = median(mut_growth_rate, na.rm = T),
    min_mut_growth = min(mut_growth_rate, na.rm = T),
    max_mut_growth = max(mut_growth_rate, na.rm = T),
    median_VAF = median(mut_growth_rate, na.rm = T),
    min_VAF = min(VAF_diff, na.rm = T),
    max_VAF = max(VAF_diff, na.rm = T),
    median_t = median(time_diff, na.rm = T),
    min_t = min(time_diff, na.rm = T),
    max_t = max(time_diff, na.rm = T)
  )
fc_stats

comp_vert_tp53 = fc_mut_groups %>% filter(Gene %in% c("PPM1D", "TP53"), timepoint_binary == 1)
p_val_ppm1d_vs_tp53_on = wilcox.test(comp_vert_tp53$mut_growth_rate, ifelse(comp_vert_tp53$Gene=="PPM1D", 1, 0))
comp_vert = fc_mut_groups %>% filter(Gene %in% c("PPM1D", "DNMT3A"), timepoint_binary == 1)
p_val_ppm1d_vs_dnmt3a_on = wilcox.test(comp_vert$mut_growth_rate, ifelse(comp_vert$Gene=="PPM1D", 1, 0))

mgr_pl = fc_mut_groups %>%
  filter(Gene %in% c("PPM1D", "DNMT3A", "TP53")) %>%
  mutate(Gene = factor(Gene, levels = c("PPM1D", "DNMT3A", "TP53"))) %>% 
  ggplot(aes(x=Gene, y=mut_growth_rate, fill=treatment)) +
  geom_hline(yintercept = 0, alpha=.7, linetype="dashed") +
  geom_boxplot( color="black", alpha=.5, outlier.shape = NA, linewidth=.5) +
  geom_bracket(
    xmin = 0.75, xmax = 1.75, y.position =-0.09,
    label = "",
    tip.length = -1,
    inherit.aes = F
  ) +
  geom_point(position=position_jitterdodge(jitter.width = .25), size=1, shape=21, alpha=.78) +
  scale_color_manual(name="Exposure to alkylating agent", values=met.brewer("Ingres")[c(5,4)]) +
  scale_fill_manual(name="Exposure to alkylating agent", values=met.brewer("Ingres")[c(5,4)]) +
  annotate("text", x=1.25, y=-.1, label = paste0("p = ", signif(p_val_ppm1d_vs_dnmt3a_on$p.value, 1) )) +
  scale_y_continuous(limits = c(-0.14,0.15)) +
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format))) ,
    vjust=-.1, method = "wilcox.test"
  ) +
  ylab("% VAF/day") +
  theme_classic() +
  gtheme(13) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    plot.title = element_text(face="bold"),
    plot.tag = element_text(size=rel(2), face="bold")
  ) +
  guides(fill = guide_legend(nrow=1, title.position = "left")) +
  ggtitle("Mutation growth rate per interval")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 3c/d: Bar plots treatment intervals
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pal = scico::scico(n=3, begin=.125, end=.9, palette="lipari")
names(pal) = c("Decrease > 5%","Increase > 5%", "Stable")

###################### UMI-sequencing
fc_mut_intervals <- fc_mut_groups %>% 
  filter(Gene == "PPM1D") %>%
  mutate(
    diff_group = factor(
      case_when( VAF_fc < 0.95  ~ "Decrease > 5%", VAF_fc > 1.05 ~ "Increase > 5%", .default = "Stable"),
      levels = c("Increase > 5%", "Stable", "Decrease > 5%")),
    treatment = factor(ifelse(timepoint_binary == 1, "On", "Off"), levels = rev(c("On", "Off")))
  ) %>%
  group_by(treatment, diff_group) %>%
  summarize(n = dplyr::n()) %>%
  mutate(perc = n/sum(n))

mut_groups_bar1 <- fc_mut_intervals %>%
  mutate(treatment = factor(ifelse(treatment=="On", "On\nn=46", "Off\nn=30"), levels = c("On\nn=46", "Off\nn=30"))) %>%
  ggplot(aes(x = treatment, y = perc, fill = diff_group)) +
  geom_bar(stat="identity", color = "black") +
  theme_classic() + gtheme(12) +
  scale_fill_manual(name = "Change of PPM1D VAF", values = pal) +
  ylab("Proportion of intervals") +
  ggtitle("UMI cohort") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face="bold"),
    axis.title.x = element_blank(),
    plot.tag = element_text(size=rel(2), face="bold")
  )

###################### Longitudinal from bulk

bulk_tp = read_excel(
  "data/bulk_longitudinal.xlsx"
)

bulk_tp_new = maf %>%
  mutate(nkey = paste0(GRId, "_", Gene_Pchange)) %>%
  left_join(bulk_tp[,c("GRId_sample", "alkylating/radiotherapy/parpi/topoisomerase")]) %>%
  filter(Gene == "PPM1D", GRId %in% names(which(table(d$GRId)>1))) %>%
  group_by(nkey) %>%
  arrange(GRId_sample) %>%
  mutate(VAF_diff = (VAF - lag(VAF)) , VAF_fc = (VAF / lag(VAF))) %>%
  ungroup() %>%
  dplyr::select(GRId_sample, Pchange, `alkylating/radiotherapy/parpi/topoisomerase`, VAF, VAF_diff, VAF_fc) %>%
  filter(!is.na(VAF_fc), !is.na(`alkylating/radiotherapy/parpi/topoisomerase`)) %>%
  mutate(
    diff_group = factor(case_when(
      VAF_fc < 0.95  ~ "Decrease > 5%", VAF_fc > 1.05 ~ "Increase > 5%", .default = "Stable"),
      levels = c("Increase > 5%", "Stable", "Decrease > 5%")
    ),
    treatment = factor(ifelse(`alkylating/radiotherapy/parpi/topoisomerase` == 1, "On\nn=11", "Off\nn=51"), levels = rev(c("On\nn=11", "Off\nn=51")))
  )

bulk_tp_binary = bulk_tp_new %>%  
  mutate(diff_group = factor(ifelse(
    diff_group  == "Increase > 5%", "Increase > 5%", "Dec"
  ), levels = c("Increase > 5%", "Dec")))

bulk_tp_binary %>%
  group_by(treatment, diff_group) %>%
  summarize(n = dplyr::n()) %>%
  mutate(perc = n/sum(n))

fisher.test(bulk_tp_binary$treatment, bulk_tp_binary$diff_group)

# length of intervals
intervals_clinical = d %>%
  filter(GRId_sample %in% bulk_tp_binary$GRId_sample) %>%
  group_by(GRId) %>%
  mutate(days = as.numeric(difftime(Date_NGS, lag(Date_NGS), units = "days"))) %>%
  dplyr::select(GRId_sample, days)

summary(intervals_clinical$days)

mut_groups_2 <- bulk_tp_new %>%  
  group_by(treatment, diff_group) %>%
  summarize(n = dplyr::n()) %>%
  mutate(perc = n/sum(n))  %>%
  mutate(treatment = factor(treatment, levels = c("On\nn=11", "Off\nn=51")))

mut_groups_bar2 = ggplot(mut_groups_2, aes(x = treatment, y = perc, fill = diff_group)) +
  geom_bar(stat="identity", color = "black") +
  theme_classic() + gtheme(12) +
  scale_fill_manual(name = "Change of PPM1D VAF", values = pal) +
  ylab("Proportion of intervals") +
  ggtitle("Clinical sequencing") +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.title.x = element_blank(),
    plot.tag = element_text(size=rel(2), face="bold")
  )

fig_3_b_c_d = plot_grid(
  mgr_pl + theme(legend.position="none"), NULL, mut_groups_bar1 + theme(legend.position="none"), NULL, mut_groups_bar2, 
  ggpubr::get_legend(mgr_pl), NULL, NULL, ggpubr::get_legend(mut_groups_bar1), NULL, 
  ncol=5, nrow=2,
  rel_widths = c( 1.6,.3,1,.3, 1), rel_heights= c(5,1)
)

ggsave2(
  "figures/main/figure_3/fig_3_b_c_d.png",
  fig_3_b_c_d,
  height=65, width=160, unit="mm", dpi=600, bg="white", scale=2
)

