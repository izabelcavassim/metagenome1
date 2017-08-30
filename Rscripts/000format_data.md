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
wb <- genefilter_sample(GlobalPatterns, filterfun_sample(function(x) x > 2), A=0.5*nsamples(GlobalPatterns))
microdf <- prune_taxa(wb, GlobalPatterns)

## Transform to even sampling depth.
microdf <- transform_sample_counts(microdf, function(x) 1E6 * x/sum(x))

##
set.seed(1)
dist <- "bray"
ord_methods <- c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
for (i in ord_methods) {
ordi <- ordinate((microdf, method=i, distance=dist)
plot_ordination(physeq, ordi, "samples", color="SampleType")
}

##
ordu <- ordinate(physeq = microdf, method = "PCoA", distance = "unifrac", weighted=TRUE)
plot_ordination(microdf, ordu, color="", shape="")
```

## Permanova
```r
set.seed(1)
# Calculate bray curtis distance matrix
microdf_bray <- phyloseq::distance(microdf, method = "bray")
sampledf <- data.frame(sample_data(microdf))
adonis(microdf_bray ~ SampleType, data = sampledf)

Call:
adonis(formula = microdf_bray ~ SampleType, data = sampledf) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SampleType  8    8.2941 1.03676  6.1773 0.74405  0.001 ***
Residuals  17    2.8532 0.16783         0.25595           
Total      25   11.1473                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
beta <- betadisper(microdf_bray, sampledf$SampleType)
permutest(beta)
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
Groups     8 0.51584 0.064480 4.3376    999  0.008 **
Residuals 17 0.25271 0.014865                        
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Reference_info
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC93182/  
https://f1000research.com/articles/5-1492/v2

