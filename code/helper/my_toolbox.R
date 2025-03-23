require(reshape2)

#-------------#
# utils
#-------------#
# pvalue to code
pval_to_signif_code <- function(p_vals) {
  return(ifelse(p_vals < 0.001, "***",
                ifelse(p_vals < 0.01, "**",
                       ifelse(p_vals < 0.05, "*", ""))))
}

#-------------#
# KM curve
#-------------#
print_and_flush <- function(string) {
    # Print a string in stdout and flush the result (useful to print in Jupyter Notebook in for loops)
    # → Arguments
    #     - string: string to print

    cat(string)
    flush.console()
}

print_p_value_from_survival_curve <- function(data, column_name, value_1_name, value_2_name) {
    # work in progress...
    stmp <- survdiff(as.formula(sprintf('Surv(os_diag_years, os_status) ~ %s', column_name)), data = data[data[,column_name] %in% c(value_1_name, value_2_name),])
    pval <- 1 - pchisq(stmp$chisq, length(stmp$n) - 1)
    
    print_and_flush(sprintf('%-21s and %21s: p = %.3f\n', value_1_name, value_2_name, pval))
}

print_p_value_from_boxplot <- function(data, column_name, value_1_name, value_2_name, clinical_name) {
    # work in progress...
    
    vec1 = data[ which(data[,column_name] %in% value_1_name), clinical_name]
    vec2 = data[ which(data[,column_name] %in% value_2_name), clinical_name]

    pval <- t.test(vec1,vec2)$p.val 
    
    print_and_flush(sprintf('%-21s and %21s: p = %.3f\n', value_1_name, value_2_name, pval))
}

km_plot <- function(
    x,
    time,
    status,
    groups,
    palette,
    title,
    pval_coord,
    limits=c(0,21),
    conf = F,
    break.by=10,
    plot_title = NULL
  ) {
  
  surv <- as.formula(paste0("Surv(", time, ",", status, ") ~", groups))
  surv_fit <- survfit(surv, data=x)
  surv_fit$call$formula <- surv
  
  print(surv_fit)
  
  g_surv <- ggsurvplot(surv_fit, data=x, 
                       color= groups, 
                       legend.title=title,
                       size=.75,
                       conf.int= conf,
                       xlim=limits,
                       break.time.by=break.by, 
                       pval=T,  pval.coord=pval_coord, pval.size=5,
                       surv.median.line="hv",
                       palette=palette,
                       risk.table=TRUE,
                       tables.col=groups, tables.y.text=FALSE,
                       risk.table.title="No. at risk",risk.table.fontsize=6,
                       ggtheme = theme_classic(),
  )
  
  g <- g_surv$plot +
    theme_classic() +
    gtheme() +
    theme(legend.key.width = unit(1.1,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_rect(linewidth=1, fill=NA),
          strip.text.x = element_text(size=10),
          plot.title = element_text(face="bold")
    ) + 
    ylab("OS Probability") + 
    guides(color = guide_legend(override.aes = list(shape = NA), nrow=2)) +
    xlab("Months after diagnosis") +
    scale_color_manual(values=palette, limits=levels(x[[groups]])) +
    ggtitle(plot_title)
  
  gt <- g_surv$table +
    theme(
      axis.line=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      legend.position="none",
      plot.title = element_text(size=18)
    ) +
    xlab("") +
    ylab("") +
    scale_color_manual(values=palette, limits=levels(x[[groups]]))
  
  g_g <- ggarrange(g + theme(legend.position="top"), gt +theme(plot.margin = unit(c(-0.5, 0, .5, .5), "cm")),
                   ncol=1, nrow=2, heights=c(5,1.5))
  
  return(g_g)
}
