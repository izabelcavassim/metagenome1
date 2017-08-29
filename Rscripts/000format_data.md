## Data
```
> library("phyloseq"); packageVersion("phyloseq")
[1] ‘1.20.0’
> library(dplyr); packageVersion("dplyr")
[1] ‘0.7.2’
> library(ggplot2)
> library(vegan)
> library(scales)
> library(grid)
> library(reshape2)
```

## 实战

```
data(GlobalPatterns)
```

### 
```
sample_names(GlobalPatterns)
ntaxa(GlobalPatterns)
nsamples(GlobalPatterns)
rank_names(GlobalPatterns)
sample_variables(GlobalPatterns)
otu_table(GlobalPatterns)
tax_table(GlobalPatterns)
phy_tree(GlobalPatterns)
taxa_names(GlobalPatterns)
```

```
## 统计样本reads数目
sample_sums(qiimedata)
## 提取某物种层次
tax_glom(GlobalPatterns, taxrank="Family")
## 
subset_samples
##
prune_taxa(taxa_sums(microdata) > 0, )

```
The included data have different number of sequences per sample. If the purpose is to make heatmaps or boxplots, it is reccomended to convert the abundances to percentages instead.
```
dfn <- transform_sample_counts(GlobalPatterns, function(x) {x/sum(x)})
```
If the purpose is to make ordination it is recommended to rarefy to the same number of reads.f doing ordination, it is recommended to first remove low abundant species.Here we keep OTUs that have been seen more than 9 times (of 10000) in at least 1 sample.
```
dsf <- filter_taxa(GlobalPatterns, function(x) max(x) >= 10, TRUE)
dfr <- rarefy_even_depth(GlobalPatterns, sample.size = 17000, rngseed = 712)
```

### Chaining all together
```
library(magrittr)
subset_samples(GlobalPatterns, SampleType != "Sediment (estuary)") %>%
  transform_sample_counts(function(x) x / sum(x) * 100) %>%
  rarefy_even_depth(sample.size = 17000, rngseed = 712) %>%
  filter_taxa(function(x) max(x) >= 10, TRUE)

GP <- GlobalPatterns %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)}) %>%  # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)  
```

```
ggplot(GP, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  facet_grid(Description~.) +
  geom_bar(stat = "identity")
```

## alpha
```
microdf <- prune_taxa(taxa_sums(qiimedata) > 0, qiimedata)
plot_richness(microdf, measures=c("Chao1", "Shannon"))
plot_richness(microdf, x="Description", measures=c("Chao1", "Shannon"))
```

## Unconstrained Ordinations

```r
## Remove OTUs that do not show appear more than 2 times in more than half the samples
wb <- genefilter_sample(qiimedata, filterfun_sample(function(x) x > 2), A=0.5*nsamples(qiimedata))
microdf <- prune_taxa(wb, qiimedata)

## Transform to even sampling depth.
microdf <- transform_sample_counts(microdf, function(x) 1E6 * x/sum(x))

dist <- "bray"
ord_methods <- c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
for (i in ord_methods) {
ordi <- ordinate((microdf, method=i, distance=dist)
plot_ordination(physeq, ordi, "samples", color="SampleType")
}

##
ordu <- ordinate(microdf, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(microdf, ordu, color="", shape="")
```




