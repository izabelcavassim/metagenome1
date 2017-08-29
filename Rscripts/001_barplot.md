## Packages
```
library("phyloseq"); packageVersion("phyloseq") #[1] ‘1.20.0’
library(dplyr); packageVersion("dplyr") #[1] ‘0.7.2’
library("reshape2")
library("data.table")
library(ggplot2)
```

## Data
```
otufile = "otu_table_tax_even.txt"
mapfile = "map.txt"
trefile = "rep_set.tre"
qiimedata <- import_qiime(otufile, mapfile, trefile)
```
```
dsn <- transform_sample_counts(qiimedata, function(x) x / sum(x))
```

  
## BARPLOT FUNCTION
```
amp_barplot <- function(data, group = "Sample", tax.aggregate = "Phylum", tax.add = NULL,
  tax.show = 10, calc = "mean", sort.by = NULL){
  data <- list(abund = as.data.frame(otu_table(data)@.Data),
               tax = data.frame(tax_table(data)@.Data, OTU = rownames(tax_table(data))),
               meta = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))

  ## Extract the data into separate objects for readability

  abund <- data[["abund"]]  
  tax <- data[["tax"]]
  meta <- data[["meta"]]
  
    suppressWarnings(
    if (!is.null(tax.add)){
      if (tax.add != tax.aggregate) {
        tax <- data.frame(tax, Display = apply(tax[,c(tax.add,tax.aggregate)], 1, paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[,tax.aggregate])
    }
  )
  
  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display = tax[,"Display"], abund) %>%
    melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()

  ## Add group information
  suppressWarnings(
    if (group != "Sample"){
      if (length(group) > 1){
        grp <- data.frame(Sample = rownames(meta), Group = apply(meta[,group], 1, paste, collapse = " "))     
      } else{
        grp <- data.frame(Sample = rownames(meta), Group = meta[,group]) 
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else{ abund5 <- data.frame(abund3, Group = abund3$Sample)}
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
      summarise(Abundance = sum(Abundance)) %>%
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
```

```
c_color <- c(brewer.pal(9,"Set1"), brewer.pal(12,"Paired"), brewer.pal(12,"Set3"))

if (nsamples(qiimedata) < 60) {
  axis_x_size=8
  width      =8
  height     =6
}
micro_theme <- function(..., base_size=12, bg='white') {
  require(grid)
  theme_classic(...) +
    theme_bw() + theme(panel.grid=element_blank()) +
    theme(panel.border=element_blank()) +
    theme(axis.line=element_line(size=0.5, colour="black"))
}

print("barplot_bygroup....")
for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  plotdata <- amp_barplot(dsn, group="Description", tax.aggregate=i, tax.add=NULL, tax.show=20, calc="mean", sort.by=NULL)
  plotdata$Display <- factor(plotdata$Display, levels=c(levels(plotdata$Display), "Others"))
  plotdata$Display[is.na(plotdata$Display)] <- "Others"
  plotdata <- plotdata[order(plotdata$Sample),]

  p <- ggplot(plotdata, aes(x= Group, y=Abundance, fill=Display)) +
              geom_bar(position="fill", stat="identity") +
              scale_y_continuous(expand=c(0, 0)) +
              labs(x="", y="Relative abundance") +
              scale_fill_manual(values=c_color[1:30]) +
              geom_hline(yintercept=c(0.25, 0.5, 0.75),color="grey50",linetype=2, size=0.3) +
              theme(axis.text.x=element_text(color="black")) +
              micro_theme() +
              guides(fill=guide_legend(ncol=1, title = NULL)) +
              theme(legend.key.size=unit(.4,'cm'), legend.key=element_rect(colour='white', fill='pink', size=1, linetype=1)) +
              theme(axis.text.x=element_text(angle=60, size=axis_x_size, vjust=0.5, color="black"))
  ggsave(paste(i, "_groupby.pdf", sep=""))
}

## barplot_bysample
print("barplot_bysample....")
for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  plotdata <- amp_barplot(dsn, group="Sample", tax.aggregate=i, tax.add=NULL, tax.show=20, calc="mean", sort.by=NULL)
  plotdata$Display <- factor(plotdata$Display, levels=c(levels(plotdata$Display), "Others"))
  plotdata$Display[is.na(plotdata$Display)] <- "Others"
  plotdata <- plotdata[order(plotdata$Sample),]

  p <- ggplot(plotdata, aes(x= Group, y=Abundance, fill=Display)) +
              geom_bar(position="fill", stat="identity") +
              scale_y_continuous(expand=c(0, 0)) +
              labs(x="", y="Relative abundance") +
              scale_fill_manual(values=c_color[1:30]) +
              geom_hline(yintercept=c(0.25, 0.5, 0.75),color="grey50",linetype=2, size=0.3) +
              micro_theme() +
              guides(fill=guide_legend(ncol=1, title = NULL)) +
              theme(legend.key.size=unit(.4,'cm'), legend.key=element_rect(colour='white', fill='pink', size=1, linetype=1)) +
              theme(axis.text.x=element_text(angle=60, size=axis_x_size, vjust=0.5, color="black"))
  ggsave(paste(i, "_sampleby.pdf", sep=""), width=width, height=height)
}
```
