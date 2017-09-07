version <- "V1.0.0"
Usage_and_Arguments <- function(author, version){
      cat(paste("
      Author: ", author, " 
      Version: ", version, "
      Description: A function to draw alpha diversity index plot that offers more control over appearance.
                   data         : (required)Phyloseq format file.
                   group        : Group the data based on a sample variable.
                   tax.aggregate: The taxonomic level that the data should be aggregated to. 
                   tax.add      : The taxonomic level that the data should be aggregated to.
                   sort.by      : Sort the barplot by a specific value of the \"group\" parameter.
      Usage:
         Rscript micro_BARplot.R <otu_table_tax_even.matrix> <group.info> <rep_set.tre>
      Example:
         Rscript micro_BARplot.R <otu_table_tax_even.txt> <map.txt> <rep_set.tre>
"))
}

args <- commandArgs(TRUE)
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help"){
   Usage_and_Arguments(author, version)
   q()
}

library("phyloseq"); packageVersion("phyloseq") #[1] ‘1.20.0’
library("dplyr"); packageVersion("dplyr") #[1] ‘0.7.2’
library("reshape2")
library("data.table")
library("ggplot2")
library("dunn.test")
library("RColorBrewer")
otufile <- args[1] #"otu_table_tax_even.txt"
mapfile <- args[2] #"map.txt"
trefile <- args[3] #"rep_set.tre"
qiimedata <- import_qiime(otufile, mapfile, trefile)

amp_diversity <- function(data, group="Sample", tax.aggregate="Phylum", tax.add=NULL, sort.by=NULL){
  data <- list(abund=as.data.frame(otu_table(data)@.Data),
               tax=data.frame(tax_table(data)@.Data, OTU=rownames(tax_table(data))),
               meta=suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  return(data)
}

data <- amp_diversity(qiimedata, group="Description", sort.by=NULL)
abund <- data[["abund"]]
tax <- data[["tax"]]
meta <- data[["meta"]]


alpha_diversity <- estimate_richness(qiimedata)
alpha_diversity <- data.frame(Sample=rownames(alpha_diversity), alpha_diversity)  #diversity(t(category.abund), index="shannon")
indexdata <-  merge(alpha_diversity, meta, by.x="Sample", by.y="X.SampleID")
write.table(indexdata, file="alpha_diversity_groupby.xls", sep ="\t", row.names=F,col.names=T, quote=F)

exit()

for (i in c("Observed", "Chao1", "ACE", "Shannon", "Simpson") {
  p <- ggplot(indexdata, aes(x=Description, y=i, color=Description)) +
              geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
              geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
              labs(x="Groups", y=paste(i, " index", sep="\t"))
  ggsave(paste(i, "_index.pdf", sep=""))
}
# kruskal.test vs aov
test_type="aov"
if (length(levels(indexdata$Description)) >= 3) {
  if(test_type="kruskal") {
    index_stats <- kruskal.test(indexdata[, i] ~ indexdata$Description)
    if (index_stats$statistic <- qchisq(0.950, index_stats$parameter)) {
      print("Refuse null hypothesis H0!!!")
    } else {
      print("Receive null hypothesis H0!!!")
    }
 #  gene_pvalue <- ddply(.data=indexdata, .variable=.(Description), .fun=function(subdata) 
 #                      { kruskal.test(value ~ Description, data=subdata )$p.value })
 # gene_pvalue <- gene_pvalue[order(gene_pvalue[,2]),]
   dunn_Test <- dunnTest(indexdata[,i], indexdata$Description, method="bonferroni")
   punadjvalue <- dunn_Test$res$P.unadj
  }

  if(test_type="aov") {
    index_stats <-  aov(indexdata[,i] ~ indexdata$Description)
    Tukey_HSD_index <- TukeyHSD(index_stats, ordered=FALSE, conf.level=0.95)
    Tukey_HSD_index_df <- as.data.frame(Tukey_HSD_index$`indexdata$Description`)
    write.table(Tukey_HSD_index_df[order(Tukey_HSD_index_df$p, decreasing=FALSE), ], file=paste(i, "alpha_stats.xls", sep=""),
                append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=TRUE,col.names=TRUE)
  }
}
#}
# wilcox.test VS t.test
#if (length(levels(indexdata$Description)) == 2 ) {
#
#}
#abund_dist <- vegdist(t(abund))
#abund_anosim <- anosim(abund_dist, meta[match(rownames(t(abund)), rownames(meta)),]$Description)
#abund_anosim <- summary(abund_anosim)

## Reference_info
#http://joey711.github.io/phyloseq/index.html
