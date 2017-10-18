library(ggpubr)

add.difference.bar <- function(g1, x1, x2, y.min, y.max, text, ylim.max = 0.5) {
  ## 提取 y 轴信息，并适当延长 y 轴
  axis.info <- ggplot_build(g1)$panel$ranges[[1]]
  y.step <- axis.info$y.major_source[2] - axis.info$y.major_source[1]
  if(is.null(axis.info)){
    axis.info <- ggplot_build(g1)$layout$panel_ranges
    y.step <- axis.info[[1]]$y.major_source[2] - axis.info[[1]]$y.major_source[1]
  }
  if ((0.5*y.step+y.max) > ylim.max)
    {g1 <- g1 + scale_y_continuous(limits =  c(y.min - y.step, y.max + y.step))}
  g1 + geom_path(data = data.frame(x = c(x1, x1, x2, x2), 
                                   y = c(0.3*y.step+y.max, 0.4*y.step+y.max, 
                                         0.4*y.step+y.max, 0.3*y.step+y.max)),
                 aes(x = x, y = y)) +
    annotate("text", label = text, x = (x1 + x2) / 2, y = 0.5*y.step+y.max)
}

indexdata = read.table("alpha_diversity_groupby.xls", header=T, sep="\t")

if (length(levels(indexdata$Description)) > 2) {
  for (i in c("Chao1", "Shannon", "Observed", "se.chao1", "ACE", "se.ACE", "Simpson", "InvSimpson", "Fisher")) {

    index_stats  <- aov(indexdata[,i] ~ indexdata$Description)
    Tukey_HSD    <- TukeyHSD(index_stats, ordered=FALSE, conf.level=0.95)
    Tukey_HSD_df <- as.data.frame(Tukey_HSD_index$`indexdata$Description`)
    Tukey_HSD_df <- cbind.data.frame(groupvs = rownames(Tukey_HSD_index_df), Tukey_HSD_index$`indexdata$Description`)
    write.table(Tukey_HSD_df[order(Tukey_HSD_index_df$p, decreasing=FALSE), ], file=paste(i, "_Tukey_stats.xls", sep=""),
                  append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=TRUE)

    wilcox_test    <- compare_means(i ~ Description, data=indexdata, method = "wilcox.test")
    wilcox_test_df <- as.data.frame(wilcox_test[,2:5])
    if (0.05 %in% wilcox_test_df$p) {print("!!!!!")}
    write.table(wilcox_test_df[order(wilcox_test_df$p, decreasing=FALSE), ], file=paste(i, "_wilcox_stats.xls", sep=""),
                  append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=F,col.names=TRUE)

    my_comparisons <- list()
    p1 <- ggboxplot(indexdata, x="Description", y=i, color = "Description", palette = rainbow(length(levels(indexdata$Description))), add = "jitter")
    p2 <- p1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns=TRUE, method = "wilcox.test") #label = "p.signif"
    p3 <- p2 +  stat_compare_means(method = "kruskal.test")
    ggsave(paste(i, "_index.pdf", sep=""))
#  dunn_Test <- dunnTest(indexdata[,i], indexdata$Description, method="bonferroni")
#  pvalue <- dunn_Test$res$P.unadj
#
#      if( pvalue < 0.05 ){
#    ## 添加差异线
#    x1 <- 1
#    x2 <- 2
#    y.min <- min(indexdata[,i])
#    y.max <- max(indexdata[,i])
#    g1 <- add.difference.bar(g1, x1, x2, y.min, y.max, 
#                             paste("p = ", as.character(round(sample.diversity.pvalue, 6), sep = "")))
#  }
  }
}

if (length(levels(indexdata$Description)) == 2) {
  for (i in c("Chao1", "Shannon", "Observed", "se.chao1", "ACE", "se.ACE", "Simpson", "InvSimpson", "Fisher")) {
        wilcox_test    <- compare_means(i ~ Description, data=indexdata, method = "wilcox.test")
    wilcox_test_df <- as.data.frame(wilcox_test[,2:5])
    if (0.05 %in% wilcox_test_df$p) {print("!!!!!")}
    write.table(wilcox_test_df[order(wilcox_test_df$p, decreasing=FALSE), ], file=paste(i, "_wilcox_stats.xls", sep=""),
                  append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=F,col.names=TRUE)
    my_comparisons <- list()
    p1 <- ggboxplot(indexdata, x="Description", y=i, color = "Description", palette = rainbow(length(levels(indexdata$Description))), add = "jitter")
    p2 <- p1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns=TRUE, method = "wilcox.test") #label = "p.signif"
    p3 <- p2 +  stat_compare_means(method = "wilcox.test")
    ggsave(paste(i, "_index.pdf", sep=""))

  }
}

