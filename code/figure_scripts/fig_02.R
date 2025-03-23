.cran_packages = c("tidyverse", "scales", "psych", "cowplot", "ggsci", "RColorBrewer",
                   "MetBrewer", "readxl", "ggpubr", "ggh4x", "ggnewscale",
                   "survival", "survminer", "patchwork", "ComplexHeatmap")

## Loading library
for (pack in .cran_packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

Sys.setlocale(locale = "fr_FR.UTF-8")

source("code/helper/ggstyles.R")
source("code/helper/my_toolbox.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD AND PREPROCESS DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Clinical data
dd = read_excel("data/dd.xlsx")
dd$Numseq <- gsub("^[A-Z]{2,4}[0-9]{3}_", "", dd$GRId_sample)
dd$Age_NGS <- as.numeric(dd$Age_NGS)
dd$DIAG2 <- factor(dd$DIAG2, levels=c('CH','CCUS','MDS','AML','MPN'))

# Enumerate time points per patient per DIAG2
dd = dd %>% arrange(gsub(" ", "", GRId_sample))
dd = dd %>% group_by(GRId, DIAG2) %>% mutate(diag_id = paste0(DIAG2, "_", 1:n())) %>% ungroup()

# Keep PPM1D at baseline
fsamples = dd %>% filter(seq_0none_1first_2others_Xnoppm1d %in% c(0,1)) %>% pull(GRId_sample)
length(fsamples)

#Subset clinical data at baseline
ddb <- dd[dd$GRId_sample %in% fsamples, ]

# Mol Data
ddmut <- read.table("data/ddmut.tsv", header=T, sep="\t", stringsAsFactors=F)
ffmut <- read.table("data/ffmut.tsv", header=T, sep="\t", stringsAsFactors=F)
maf <- read_xlsx("data/maf_final.xlsx")

# Remove ALFA009_1	c.1261-8T>G	10
maf = maf[which(!maf$Sample_Gene_Pchange == "ALFA009_1_PPM1D_splice_exon5" ),]

# Subset molecular data at baseline
ddmutb = ddmut[ddmut$GRId_sample %in% fsamples , ]
ffmutb = ffmut[ffmut$GRId_sample %in% fsamples , ]
mafb = maf[maf$GRId_sample %in% fsamples , ]

nrow(mafb) == 365

mafb = left_join(mafb, ddb[,c("GRId", "DIAG2", "Sex", "t21")], by="GRId")
mafb$DIAG2 <- factor(mafb$DIAG2, levels=levels(ddb$DIAG2))

# PPM1D mutation only
ppmaf = mafb[which(mafb$Gene=="PPM1D"),]

# Enumerate
ddb = left_join(ddb, ppmaf[,c("GRId_sample","Gene")] %>%
                  group_by(GRId_sample) %>%
                  summarize(num_ppm1d=n()), by="GRId_sample")

ddb$DIAG2 <- factor(ddb$DIAG2, levels=c("CH","CCUS","MDS","AML","MPN"))

ddb = ddb %>% mutate(
  ATCD_F=recode(ATCD_K, "Ovaire"="Gynecologic", "Endometre"="Gynecologic", "Trompe"="Gynecologic", "Sein" = "Breast cancer",
                "Tumeur neuroendocrine"="Other", "Col de l'uterus" = "Gynecologic", 
                "MDS/LA" = "Other", "Hemopathie lymphoide" = "Lymphoid", "SMP" = "Other",
                "Aucun" = "de novo hematological",
                "Paragangliome" = "Other", "Lynch" = "Other", "Autre" = "Other",
                "Prostate" = "Prostate cancer", "Osteosarcome" = "Osteosarcoma", "Poumon" = "Lung cancer"),
  tmn = ifelse(DIAG2 %in% c("MDS", "AML"), "tMN", "CH/CCUS/MPN"),
  latency = ifelse(
    DIAG2 == "MPN",
    as.numeric(difftime(Date_NGS, Date_diagnosis, unit="days"))/30.417,
    as.numeric(difftime(Date_diagnosis, Date_cancer, unit="days"))/30.417
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Description
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

summary(ppmaf$VAF)

ddb = ddb %>% mutate(STUDY = ifelse(grepl("ALFA", GRId), "ALFA", "GR"))

ddb %>%
  group_by(STUDY, DIAG2) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n))

mafb = mafb %>% mutate(STUDY = ifelse(grepl("ALFA", GRId), "ALFA", "GR"))

mafb %>%
  filter(Gene == "PPM1D") %>%
  group_by(STUDY) %>%
  summarize(median = median(VAF), min = min(VAF), max=max(VAF))

mafb %>%
  mutate(GRId = factor(GRId)) %>%
  filter(Gene != "PPM1D") %>%
  summarize(n = n(), .by="GRId") %>%
  complete(GRId, fill=list(n=0)) %>%
  mutate(STUDY = ifelse(grepl("ALFA", GRId), "ALFA", "GR")) %>%
  group_by(STUDY) %>%
  summarize(median=median(n), min=min(n), max=max(n))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2a: Overall survival in diagnostic groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ddb_surv = ddb %>%
  filter(!is.na(Date_LFU)) %>%
  mutate(
    os = ifelse(DIAG2 == "MPN",
                interval(Date_NGS, Date_LFU) %/% months(1),
                interval(Date_diagnosis, Date_LFU) %/% months(1)
    ),
    Death = ifelse(Death == "oui", 1 , 0)
  )

os_diag = km_plot(
  ddb_surv,
  "os",
  "Death",
  "DIAG2",
  anno_colour$Diagnosis,
  "Diagnosis",
  c(27.5, .99),
  limits=c(0,33),
  break.by = 6
)
os_diag

ggsave2(
  "figures/main/figure_2/fig_2a.png",
  width=110, height=120, dpi=400, bg="white", units = "mm", scale = 1.5
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2b: Lolliplot of all PPM1D mutations
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Add AA position
ppmaf$AAnumber <- sapply(ppmaf$Pchange, function(x) {
  y <- strsplit(x, split="_")[[1]][1]
  yy <- as.numeric(strsplit(y, split="")[[1]])
  if (any(!is.na(y))) {
    res <- as.numeric(paste(yy[!is.na(yy)],collapse=""))
  } else {
    res <- NA
  }
})

ppmaf %>% filter(Pchange == "R552*") %>% group_by(DIAG2) %>% summarize(n=n())

ppmaf_ch_ccus <- ppmaf %>% filter(DIAG2 %in% c("CH", "CCUS"))
ppmaf_ch_ccus$Ylolli <- sapply(1:nrow(ppmaf_ch_ccus), function(i) {
  tmp <- ppmaf_ch_ccus %>% filter(AAnumber==(ppmaf_ch_ccus[i,] %>% pull(AAnumber)))
  tmp <- tmp[order(tmp$VAF),]
  # jres <- match( (ppmaf_ch_ccus_mpn[i,] %>% pull(Pchange)) , tmp$Pchange )
  jres <- length(tmp$Pchange)
  return(jres)
})

ppmaf_ch_ccus <- ppmaf_ch_ccus %>% mutate(size = ifelse(Ylolli>2, Ylolli, NA), label = ifelse(Ylolli>2, Pchange, NA)) %>% group_by(AAnumber) %>% filter(row_number() == 1)

ppmaf_tmn <- ppmaf %>% filter(DIAG2 %in% c("MPN", "MDS", "AML"))
ppmaf_tmn $Ylolli <- sapply(1:nrow(ppmaf_tmn), function(i) {
  tmp <- ppmaf_tmn %>% filter(AAnumber==(ppmaf_tmn[i,] %>% pull(AAnumber)))
  tmp <- tmp[order(tmp$VAF),]
  # jres <- match( (ppmaf_tmn[i,] %>% pull(Pchange)) , tmp$Pchange )
  jres <- length(tmp$Pchange)
  return(jres)
})
ppmaf_tmn <- ppmaf_tmn %>% mutate(size = ifelse(Ylolli>1, Ylolli, NA), label = ifelse(Ylolli>1, Pchange, NA)) %>% group_by(AAnumber) %>% filter(row_number() == 1)

# PPM1D GENE STRCUCTURE
ppppm1d = data.frame(
  x1=c(390,421),
  x2=c(421, 605),
  y1=c(-.4,-.4),
  y2=c(.4, .4),
  domain=c("", "C-Terminus (Exon 6)")
)

#~~~~~  PLOT V1
lolliplot = ggplot() + theme_classic() +# PPM1D gene and domains
  geom_rect(data=ppppm1d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = domain, color = domain), color = NA,
            alpha=0.9) +
  geom_text(data=ppppm1d, aes(label=domain, x = (x1 + x2)/2, y = 0), color="white", size=4) +
  scale_fill_manual(values = c("grey", derain[7]), guide=FALSE) +
  scale_color_manual(values = c("grey", derain[7]), guide=FALSE) +
  ggnewscale::new_scale_fill() +
  
  # CH/CCUS/MPN
  geom_segment(data=ppmaf_ch_ccus , aes(x=AAnumber, xend = AAnumber, y = .4, yend = Ylolli, color=Effect)) +
  geom_point(data=ppmaf_ch_ccus , aes(x=AAnumber, y=Ylolli, fill=Effect, size = Ylolli), pch=21) +
  ggrepel::geom_text_repel(data=ppmaf_ch_ccus, aes(x=AAnumber, y=Ylolli, label = label), direction = "y" , nudge_y =1.2, segment.alpha = 0 ) +
  geom_text(data=ppmaf_ch_ccus, aes(x=AAnumber, y=Ylolli, label = size), color ="white") +
  
  # MPN/MDS/AML
  geom_segment(data=ppmaf_tmn , aes(x=AAnumber, xend = AAnumber, y = -.4, yend = -Ylolli, color=Effect)) +
  geom_point(data=ppmaf_tmn , aes(x=AAnumber, y=-Ylolli, fill=Effect, size = Ylolli, label = size), pch=21) +
  # ggrepel::geom_text_repel(data=ppmaf_tmn, aes(x=AAnumber, y=-Ylolli, label = label), direction = "y" , nudge_y = -1, segment.alpha = 0 ) +
  geom_text(data=ppmaf_tmn, aes(x=AAnumber, y=-Ylolli, label = size), color ="white") +
  
  # Adjust aes
  scale_size_continuous(range = c(5, 9)) +
  scale_fill_manual(name = "PPM1D mutation", values = c("#6699CC", "#997700")) +
  scale_color_manual(values = c("#6699CC", "#997700")) +
  guides(color="none", size="none", fill = guide_legend(override.aes=list(size=4.5))) +
  
  scale_x_continuous(limits=c(390, 624), breaks = seq(400, 605, 25), expand = c(0, 0)) +
  ylim(-9, 9) +
  xlab("Amino acid position") +
  theme_classic() +
  theme(axis.text=element_text(size=12,colour="black"),
        axis.title=element_text(size=13,colour="black"),
        axis.title.y = element_blank(),
        axis.line.y=element_line(),
        legend.text=element_text(size=12,colour="black"),
        legend.title=element_text(size=15,colour="black"),
        legend.position = c(0.1, 0.9)
        
  ) +
  theme(axis.ticks.y=element_blank(), axis.line.y=element_blank(), axis.text.y=element_blank()) +
  annotate("label", x = 600, y = 7, label = "CH/CCUS", size=5) +
  annotate("label", x = 600, y = -4, label = "MPN/MDS/AML", size=5)

ggsave2(
  "figures/main/figure_2/fig_2b.png",
  lolliplot,
  width = 240, height = 150, bg = "white", units = "mm", dpi=400, scale=1.1
)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2c: VAF of PPM1D mutations in diagnostic groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

num_ppm1d_df = ppmaf %>%
  group_by(GRId) %>%
  summarize(num_ppm1d = n())

ddb = ddb %>%
  left_join(num_ppm1d_df) %>%
  mutate(range_ppm1d = factor(ifelse(num_ppm1d>2, 1, 0), levels = c(0, 1)))

ddb %>%
  mutate(DIAG3 = ifelse(DIAG2 %in% c("CH", "CCUS"), "CH/CCUS", "MPN/MDS/AML")) %>%
  group_by(DIAG3, range_ppm1d) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n))

vaf_diag_pl = ggplot(ppmaf, aes(x=DIAG2, y=VAF)) +
  geom_violin(fill=derain[7], alpha=.3, linewidth=.5) + 
  geom_jitter(height = 0, width = 0.2, size=1, color=derain[7], alpha=0.8) +
  stat_summary(
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    color = "darkorange", size=.55, linewidth=0.5) +
  theme_classic() + gtheme(11) + 
  ylab("VAF (%) PPM1D") + xlab("Diagnosis") + 
  stat_compare_means(label = "p.format", size = 4, labels = function(x) format(x, scientific = TRUE)) +
  scale_y_continuous(limits = c(0,50)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

ggsave2(
  "figures/main/figure_2/fig_2c.png",
  vaf_diag_pl,
  width = 55, height = 60, dpi = 500, bg = "white", units = "mm", scale = 1.5
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2d: Number of PPM1D mutations in diagnostic groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

n_mut_diag_pl = ddb %>%
  mutate(DIAG2 = factor(DIAG2, levels = rev(levels(DIAG2)))) %>%
  ggplot(aes(x=num_ppm1d, fill = DIAG2)) + geom_bar(color="black") + 
  theme_classic() + gtheme(20) + 
  ylab("Number of patients") + xlab("Number of PPM1D mutations") +
  scale_fill_manual(values = rev(pal_jama("default")(5)), name = "Diagnosis") +
  theme(legend.position = c(0.5, 0.8))

ggsave2(
  "figures/main/figure_2/fig_2d.png",
  n_mut_diag_pl,
  width = 100, height = 112.5, dpi = 500, bg = "white", units = "mm", scale = 1.4
)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2e: Number of co-mutations in diagnostic groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Comutation number
mm = mafb[mafb$Gene!="PPM1D",] %>% group_by(GRId) %>% summarize(num_comut=n())
ddb = left_join(ddb, mm[,c("GRId","num_comut")], by="GRId")
ddb[is.na(ddb$num_comut),"num_comut"] = 0

# Range of co-mutations
ddb$range_comut2 = ddb$num_comut
ddb$range_comut2[ddb$num_comut>2] <- ">2"
ddb$range_comut2 = factor(ddb$range_comut2, levels=c("0","1","2",">2"))

# Add frequencies of group
freq_list <- data.frame(table(ddb$DIAG2)) %>% mutate(labels = paste0(Var1, "\n", "n= ", Freq)) %>%
  pull(labels)

n_comut_diag_pl = ggplot(ddb, aes(fill=range_comut2, x=DIAG2)) +
  geom_bar(color="black", position= position_fill(reverse = T), size=.5) +
  theme_classic() + gtheme(12) + 
  theme(axis.title.y = element_blank()) +
  xlab("Diagnosis") + 
  ylab("Proportion of patients") + 
  topleg +
  scale_x_discrete(labels = freq_list) +
  scale_fill_manual(name="Number of co-mut.", values=demuth[5:8]) + coord_flip() +
  guides(fill =guide_legend(override.aes = list(size=3)))

ggsave2(
  "figures/main/figure_2/fig_2e.png",
  n_comut_diag_pl,
  width = 65, height = 55, dpi = 500, bg = "white", units = "mm", scale = 1.6
)

# Stats
ddb_comp = ddb %>%
  filter(DIAG2 != "MPN") %>%
  mutate(
    tmn = ifelse(DIAG2 %in% c("CH", "CCUS"), "CH/CCUS", "MDS/AML"),
    comut_binary = ifelse(num_comut>0, 1, 0)
  )
fisher.test(ddb_comp$tmn, ddb_comp$comut_binary)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2e: Oncoplot at baseline
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Splitting per diagnosis
cl_df <- ddb %>% group_by(GRId) %>% arrange(DIAG2) %>% filter(row_number()==1) %>%
  ungroup() %>% select(GRId, DIAG2) %>% tibble() %>% column_to_rownames(var = "GRId")
colnames(cl_df)[1] <- "Diagnosis"
cl_df$Diagnosis <- factor(cl_df$Diagnosis)

anno_colour <- list(
  Diagnosis = c(CH=jama_pal[1], CCUS=jama_pal[2], MDS=jama_pal[3], AML=jama_pal[4], MPN=jama_pal[5])
)

rownames(ddmutb) <- NULL
genes = ddmutb[c(1,4:ncol(ddmutb))] %>% column_to_rownames("GRId") %>%
  t() %>%
  data.frame() %>%
  mutate(sumVar = rowSums(across(1:112))) %>%
  mutate(per=sumVar/112) %>%
  arrange(-per) %>%
  dplyr::select(per) %>%
  rownames_to_column("Gene") %>%
  filter(per>0.03 | Gene == "CALR") %>%
  select("Gene")

mafb$range_vaf <- cut(mafb$VAF, breaks = c(0,5,10,50,100), labels = c("1 to 5%","5 to 10%","10 to 50%","more than 50%"), include.lowest = F)
mafb$range_vaf <- factor(mafb$range_vaf, levels=rev(c("1 to 5%","5 to 10%","10 to 50%","more than 50%")))

## Range for palette
pal_grad = rev(brewer.pal(4, "RdYlBu"))

maf_heat_2 = mafb %>% filter(Gene %in% genes$Gene) %>% group_by(GRId)
maf_heat_2$Gene = as.character(maf_heat_2$Gene)
maf_heat_2$Gene = factor(maf_heat_2$Gene, levels=rev(names(sort(table(maf_heat_2$Gene)))))

# Matrix for Heatmap
heat_2 <- matrix(nrow=length(levels(maf_heat_2$Gene)),
                 ncol = length(unique(mafb$GRId)), dimnames=list(levels(maf_heat_2$Gene), rownames(cl_df)))
heat_2[is.na(heat_2)] <- ""

## Max VAF
mafc_2 <- maf_heat_2 %>% group_by(GRId, Gene) %>% slice(which.max(VAF))

## VAF range
for (row in c(1:nrow(mafc_2))) {
  heat_2[mafc_2[row,]$Gene, mafc_2[row,]$GRId] <- paste0(heat_2[mafc_2[row,]$Gene, mafc_2[row,]$GRId], mafc_2[row,]$range_vaf)
}

meanVAF_df_2 = maf_heat_2 %>% group_by(Gene) %>% summarize(meanVAF=mean(VAF))

col = c("1 to 5%" = pal_grad[1], "5 to 10%" = pal_grad[2], "10 to 50%" = pal_grad[3], "more than 50%" = pal_grad[4])
# col = c("0 to 5%" = "blue", "5 to 10%" = "red", "10 to 50%" = "green", "more than 50%" = "orange")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#ececec", col = NA))
  },
  # big blue
  "1 to 5%" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["1 to 5%"], col = NA))
  },
  # big red
  "5 to 10%" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["5 to 10%"], col = NA))
  },
  "10 to 50%" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["10 to 50%"], col = NA))
  },
  "more than 50%" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["more than 50%"], col = NA))
    
  }
)

# Cytogenetic annotation
cytog <- ddb %>% dplyr::select(GRId, Complexe, del17_17p, del7_7q, del5q) %>%
  column_to_rownames("GRId") %>% data.frame() # %>% replace(is.na(.), 0)
colnames(cytog) <- c("CK", "del17/17p", "del7/7q", "del5q")
cytog = cytog %>% transmute_all(funs(ifelse(. == 1, deparse(substitute(.)), NA)))

# Get order of Genes
anno_order_2 <- oncoPrint(as.matrix(heat_2), alter_fun=alter_fun)@row_order

anno_colour$Cytog <- c(derain[3], derain[4], derain[5], derain[7])
names(anno_colour$Cytog) <- c("CK", "del17/17p", "del7/7q", "del5q")

top_anno <- HeatmapAnnotation(
  Diagnosis = as.matrix(cl_df),
  col=anno_colour,
  border=T,
  na_col="white",
  # show_legend=c(Diagnosis = F),
  show_annotation_name = c(Diagnosis=FALSE),
  annotation_name_gp=gpar(fontface=1, fontsize=16),
  annotation_legend_param = list(
    Diagnosis = list( 
      title_gp = gpar(cex=1.25, fontface = 1), 
      at = levels(cl_df$Diagnosis),
      labels_gp = gpar(fontsize = 16))
  ),
  simple_anno_size = unit(6, "mm")
)

bottom_anno <- HeatmapAnnotation(
  Cytog = as.matrix(cytog[rownames(cl_df),]),
  border=T,
  col = anno_colour,
  na_col="white",
  show_legend=c(Cytog=F),
  annotation_name_gp=gpar(fontface=1, fontsize=16),
  simple_anno_size = unit(6, "mm")
)

heatmap_legend_param = list(title = "MAX VAF", at = c('1 to 5%', '5 to 10%', '10 to 50%', 'more than 50%'), 
                            labels = c('1 to 5%', '5 to 10%', '10 to 50%', 'more than 50%'),
                            labels_gp = gpar(fontsize=16),
                            title_gp=gpar(fontsize=16)
)

h2 = oncoPrint(as.matrix(heat_2),
               column_split = cl_df,
               top_annotation = top_anno,
               bottom_annotation = bottom_anno,
               right_annotation = rowAnnotation(
                 `% Mean VAF` = anno_points(meanVAF_df_2[anno_order_2,]$meanVAF,         
                                            axis_param = list(side = "top", labels_rot = 45, at = c(0, 25)),
                                            width=unit(2.5, "cm"),
                 ),
                 row_barplot = anno_oncoprint_barplot(
                   border = F, width = unit(3.5, "cm"), show_fraction = T,
                   axis_param = list(side = "top", labels_rot = 45),
                 ),
                 annotation_name_rot = 0,
                 annotation_name_side="top"
               ),
               show_column_names=F,
               column_title_gp = gpar(fontsize = 16, fontface="bold"),
               row_names_gp = gpar(fontface = "italic"),
               alter_fun = alter_fun, 
               alter_fun_is_vectorized = T,
               col = col, 
               remove_empty_columns = F,
               remove_empty_rows = F,
               heatmap_legend_param = heatmap_legend_param,
               show_heatmap_legend = T,
               pct_gp = gpar(fontsize = 16)
)

lgd = Legend(
  title = "MAX VAF",
  at = c('1 to 5%', '5 to 10%', '10 to 50%', 'more than 50%'),
  legend_gp = gpar(fill = col),
  labels_gp = gpar(fontsize=16),
  title_gp=gpar(fontsize=16, fontface="bold"),
)

h2_pl = draw(h2, merge_legends=T)

png("figures/main/figure_2/fig_2f.png", width = 300, height = 150, bg = "white", units = "mm", res=600)
h2_pl
dev.off()
