require(ggplot2)
require(MetBrewer)
require(ggsci)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Basic styles
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### ggplot styles
theme1 = theme_bw() + theme(text=element_text(face="bold", size=20))
theme0 = theme_bw() + theme(text=element_text(face="bold", size=10))
noleg = theme(legend.position = "none")
topleg = theme(legend.position = "top")
angle45 = theme(axis.text.x = element_text(angle = 45, hjust = 1))
noxtitle = theme(axis.title.x = element_blank())
noytitle = theme(axis.title.y = element_blank())
noxlabel = theme(axis.text.x = element_blank())
noylabel = theme(axis.text.y = element_blank())
nolegtitle = theme(legend.title = element_blank())

legtwoshape = guides(shape=guide_legend(ncol=2))
legtwofill = guides(fill=guide_legend(ncol=2))


bold15  = theme(axis.text = element_text(size=15,face="bold"))
bold20  = theme(axis.text = element_text(size=20,face="bold"))

gtheme <- function(size=18) {
    theme(text=element_text(size=size,colour="black"),
          axis.text=element_text(size=size,colour="black"),
          axis.title=element_text(size=size,colour="black"),
          strip.text.x = element_text(size=size),
          legend.text=element_text(size=size))
}

#colour palettes
derain <- met.brewer("Derain")
demuth <- met.brewer("Demuth")
degas <- met.brewer("Degas")
monet <- met.brewer("Monet")
hokusai <- met.brewer("Hokusai1")
cassatt <- met.brewer("Cassatt1")

jama_pal <- pal_jama("default")(7)


## Specific annotations PPM1D project
anno_colour <- list(
    Diagnosis = c(CH=jama_pal[1], CCUS=jama_pal[2], MDS=jama_pal[3], AML=jama_pal[4], MPN=jama_pal[5])
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Longitudinal data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# PPM1D 12
# TP53 3
# TET2 5
# DNMT3A 3
ippm1d <- c( met.brewer("Manet")[6:11] , rev(met.brewer("Hokusai2")) ) # Blue
itp53 <- c( met.brewer("Benedictus")[1:6] , met.brewer("Cassatt1")[1:4] ) # Pink
itet2 <- met.brewer("VanGogh3") # Green
idnmt3a <- met.brewer("Tam") # Yellow
iother <- c( met.brewer("Isfahan1")[1:4] , met.brewer("Demuth")[6:10] ) # Brown and Grey

# Col dict
col.treat = c(
  "Carbo" = "#9f6e71", "Oxaliplatine" = "#6f4e50", "Taxol" = "#7d87b2", "Caelyx" = "#629ae7" ,"Avastin" = "#2a385e",
  "Gemzar" = "#cf9023", "Atezolizumab" = "#cf7023", "Niraparib" = "#d9768a", "Olaparib" = "#c399a2",
  "Fulvestrant" = "#abccbe", "Letrozole" = "#749e89", "Olaparib Placebo" = "#3c3c3c", "None" = "#e9e9e9"
)

## Colours treatplots
col.treat2 = c(
    # Pt-based - brown
    "CB" = "#9f6e71", "OP" = "#6f4e50", "PC" = "#a98b8b",
    # Alkylating - blue/green
    "DXR" = "#629ae7", "VYXEOS" = "#566eaf", "3+7 (IDA)"="#a6d383", "CYC"="#0493a9", "IFO"="#62a768",
    # Immunotherapy purple
    "ATZ"="#010652",
    # VEGFa
    "AVA" = "#542d6c", 
    # Targeted red/pink
    "NIR" = "#d9768a", "OLA" = "#bc415c",
    # Nucleobase analogues red
    "GEM"="#9d0514", "Hydroxyurea" = "#e70a20",
    # Hormonotherapy
    "FULV" = "#abccbe", "LET" = "#749e89",
    # Leukemia-directed yellow
    "AZA"="#a85c06", "AZA/VEN"="#e0be12", "MEL"="#b0a56a",
    # Other grey
    "Trabectedin"="#6b6061", "MTX"="#6e7775", "ETOP"="#5b5c60",  "ETOP/Ara-C"="#7e7f84", 
    "None"="#e9e9e9"
)

ExpandColors <- function(colors, n, steps = 11){
  if(n <= steps){
    suppressWarnings({
      sapply(colors, function(x){colorRampPalette(c(x, "#000000"))(steps)}) %>% 
        as.data.frame() %>% 
        filter(row_number() <= n) %>% 
        gather(key = original.color, value = expanded.color)
    })
  }else{
    warning("Select n < steps!")
  }
}

# Gene Color Interpolation
gene_colors = c(
  "#FC822AFF",
  ATM="pink2",
  "linen",
  CBL="grey30",
  "black",
  DNMT3A="goldenrod1",
  PPM1D="#6bc5e3",
  TET2="#749e89",
  TP53="#CD534CFF"
)

.colorInterpolation <- function(x, mutations, map) {
  col <- ExpandColors(gene_colors[[x]], n=length(mutations), steps=7)$expanded.color
  names(col) <- map[[x]]
  return(col)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clonal groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Colors PPM1D clonal groups
pal_npg <- pal_npg("nrc")(4)

anno.clonality.all <- c(
  "Dominant" = pal_npg[1],
  "Co-dominant" = "#e47a7a",
  "Co-/dominant"=pal_npg[1],
  "Subclonal" =pal_npg[2],
  "Bystander" = pal_npg[4],
  "Only PPM1D" = pal_npg[3] 
)

