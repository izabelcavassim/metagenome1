setwd("C:/Users/Administrator/Desktop/xiuamp")
library("phyloseq")
library("magrittr")
library("vegan")
library("ape")
library("ggdendro")
library("ggplot2")
library("dplyr")
library("reshape2")
library("data.table")
load("DNAext_1.0.RData")


#rep <- subset_samples(V13, Exp.Biorep == "YES") %>%
#  rarefy_even_depth(sample.size = 17000, rngseed = 712) %>%
#  filter_taxa(function(x) max(x) >= 10, TRUE)
  
bbn <- subset_samples(V13,Exp.beadbeating == "YES") %>%
  transform_sample_counts(function(x) x / sum(x) * 100)
  

  amp_barplot <- function(data, group = "Sample", tax.aggregate = "Phylum", tax.add = NULL, tax.show = 10, calc = "mean", sort.by = NULL){
  data <- list(abund = as.data.frame(otu_table(data)@.Data),
               tax = data.frame(tax_table(data)@.Data, OTU = rownames(tax_table(data))),
               meta = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))

  ## Clean up the taxonomy
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
        grp <- data.frame(Sample = rownames(sample), Group = apply(sample[,group], 1, paste, collapse = " "))     
      } else{
        grp <- data.frame(Sample = rownames(sample), Group = sample[,group]) 
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
amp_barplot(bbn, group = "Sample", tax.aggregate = "Phylum", tax.add = NULL, tax.show = 10, calc = "mean", sort.by = NULL)
ggplot(abund7, aes(x= Group, y = Abundance, fill = Display)) + geom_bar(position = "fill", stat="identity")

  
