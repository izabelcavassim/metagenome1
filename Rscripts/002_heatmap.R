version <- "V1.0.0"
Usage_and_Arguments <- function(author, version){
      cat(paste("
      Author: ", author, " 
      Version: ", version, "
      Description: A function to draw a barplot that offers more control over appearance.
                   data         : (required)Phyloseq format file.
                   group        : Group the data based on a sample variable.
                   tax.aggregate: The taxonomic level that the data should be aggregated to. 
                   tax.add      : The taxonomic level that the data should be aggregated to.
                   tax.show     : The number of taxa to show or a vector of taxa names (default: 30).
                   calc         : The method of calculate across the groups (default: \"mean\").
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
library("RColorBrewer")
library(pheatmap)

otufile <- args[1] #"otu_table_tax_even.txt"
mapfile <- args[2] #"map.txt"
trefile <- args[3] #"rep_set.tre"
qiimedata <- import_qiime(otufile, mapfile, trefile)
dsn <- transform_sample_counts(qiimedata, function(x) x / sum(x))

amp_barplot <- function(data, group="Sample", tax.aggregate="Phylum", tax.add=NULL, tax.show=10, calc="mean", sort.by=NULL){
  data <- list(abund=as.data.frame(otu_table(data)@.Data),
               tax=data.frame(tax_table(data)@.Data, OTU=rownames(tax_table(data))),
               meta=suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  ## Clean up the taxonomy
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  meta <- data[["meta"]]

    suppressWarnings(
    if (!is.null(tax.add)){
      if (tax.add != tax.aggregate) {
        tax <- data.frame(tax, Display=apply(tax[,c(tax.add,tax.aggregate)], 1, paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display=tax[,tax.aggregate])
    }
  )

  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display=tax[,"Display"], abund) %>%
    melt(id.var="Display", value.name= "Abundance", variable.name="Sample")
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>%
    as.data.frame()

  ## Add group information
  suppressWarnings(
    if (group != "Sample"){
      if (length(group) > 1){
        grp <- data.frame(Sample=rownames(meta), Group=apply(meta[,group], 1, paste, collapse=" "))
      } else{
        grp <- data.frame(Sample=rownames(meta), Group=meta[,group])
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else{ abund5 <- data.frame(abund3, Group=abund3$Sample)}
  )

  ## Take the average to group level
  if (calc == "mean"){
    abund6 <- data.table(abund5)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>%
      as.data.frame()
  }

  ## Find the X most abundant levels
  if (calc == "mean"){
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance=sum(Abundance)) %>%
      arrange(desc(Abundance))
  }

  ## Subset to X most abundant levels
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){
      tax.show <- nrow(TotalCounts)
    }
    abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show])
  }

  ## Subset to a list of level names
  if (!is.numeric(tax.show)){
    if (tax.show != "all"){
      abund7 <- filter(abund6, Display %in% tax.show)
    }
    ### Or just show all  
    if (tax.show == "all"){
      tax.show <- nrow(TotalCounts)
      abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show])
    }
  }
  abund7 <- as.data.frame(abund7)
  return(abund7)
}


################
## heatmaplot_bygroup

for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  heatmapdata <- amp_barplot(qiimedata, group="Description", tax.aggregate=i, tax.add=NULL, tax.show=100, calc="mean", sort.by=NULL)
  heatmapdata$Display <- factor(heatmapdata$Display, levels=c(levels(heatmapdata$Display), "Others"))
  heatmapdata$Display[is.na(heatmapdata$Display)] <- "Others"
  heatmapdata <- heatmapdata[order(heatmapdata$Sample),]
  heatmapdatawide <-dcast(heatmapdata, Display ~ Group, value.var = "sum", fun.aggregate=sum)
  rownames(heatmapdatawide) <- heatmapdatawide[,1]
  hm_matrix <- heatmapdatawide[, -1]
  hm <- apply(hm_matrix,2, function(x){log10(x/sum(x) + 1)})
  write.csv(hm, paste("normalized",i, "_heatmap.xls", sep = ""))
  system("sed -i 's/,/\t/g' *_heatmap.xls")
  if (nrow(hm_matrix) <= 50) {
    pheatmap(hm, scale="row", clustering_distance_rows="correlation", cluster_cols=F,
             fontsize = 4, clustering_method="complete", filename=paste(i, "_heatmap.pdf", sep = ""))
    pheatmap(hm, scale="row", clustering_distance_rows="correlation", cluster_cols=T,
             fontsize = 4,clustering_method="complete", filename=paste(i, "_heatmap_clustercols.pdf", sep = ""))
  }

  if(nrow(hm_matrix) > 50) {
    pheatmap(hm, scale="row", clustering_distance_rows="correlation", cluster_cols=F,
             clustering_method="complete", filename=paste(i, "_heatmap.pdf", sep = ""),
             fontsize = 4, cellheight=5)
    pheatmap(hm, scale="row", clustering_distance_rows="correlation", cluster_cols=T,
             clustering_method="complete", filename=paste(i, "_heatmap_clustercols.pdf", sep = ""),
             fontsize = 4, cellheight=5)
  }
}




## Parameter
## heatmap():
##  color = colorRampPalette(c("blue4","blue","white", "red","red4"))(nrow(hm))
##   cellwidth=60, cellheight=5

## Reference_info
#http://joey711.github.io/phyloseq/index.html

