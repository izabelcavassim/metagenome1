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
```
microdata <- filter_taxa(qiimedata, function(x) max(x) >= 1, TRUE) %>%
rarefy_even_depth()

microdata <- prune_taxa(taxa_sums(microdata) > 0, microdata)
```

```
data <- qiimedata
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

```