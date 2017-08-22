## Data
```
> library("phyloseq"); packageVersion("phyloseq")
[1] ‘1.20.0’
> library(dplyr); packageVersion("dplyr")
[1] ‘0.7.2’
```

## 创建phyloseq文件
```
otufile = "otu_table_taxonomy.txt"
mapfile = "map.txt"
trefile = "rep_set.tre"

qiimedata <- import_qiime(otufile, mapfile, trefile)

sample_names(qiimedata)
ntaxa(qiimedata)
nsamples(qiimedata)
rank_names(qiimedata)
sample_variables(qiimedata)
otu_table(qiimedata)
tax_table(qiimedata)
phy_tree(qiimedata)
taxa_names(qiimedata)
```
## 实战

```
data(GlobalPatterns)
```

### 筛选某组样本|删除丰度小的otu|抽平
```
rep <- subset_samples(V13, Exp.Biorep == "YES") %>%
  filter_taxa(function(x) max(x) >= 10, TRUE) %>%
  rarefy_even_depth(sample.size = 17000, rngseed = 712)
```
### 

```
## alpha

```
microdf <- prune_taxa(taxa_sums(qiimedata) > 0, qiimedata)
plot_richness(microdf, measures=c("Chao1", "Shannon"))
plot_richness(microdf, x="Description", measures=c("Chao1", "Shannon"))
```
## beta
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

```


# Load data
MiDAS data  http://midasfieldguide.org/
```
library("phyloseq"); packageVersion("phyloseq")
GP <- load("https://github.com/xiucz/metagenome1/blob/master/data/GlobalPatterns.RData")
```
# Subset to the relevant samples
```
rep <- subset_samples(GP, SampleType != "Sediment (estuary)")
```

The included data have different number of sequences per sample. If the purpose is to make heatmaps or boxplots, it is reccomended to convert the abundances to percentages instead.
```
dfn <- transform_sample_counts(rep, function(x) x / sum(x) * 100)
```

If the purpose is to make ordination it is recommended to rarefy to the same number of reads.f doing ordination, it is recommended to first remove low abundant species.Here we keep OTUs that have been seen more than 9 times (of 10000) in at least 1 sample.
```
dfr <- rarefy_even_depth(rep, sample.size = 17000, rngseed = 712)
dsf <- filter_taxa(rep, function(x) max(x) >= 10, TRUE)
```

## Chaining all together
```
library(magrittr)
subset_samples(GP, SampleType != "Sediment (estuary)") %>%
transform_sample_counts(function(x) x / sum(x) * 100) %>%
rarefy_even_depth(sample.size = 17000, rngseed = 712) %>%
filter_taxa(function(x) max(x) >= 10, TRUE)
```
# convert to a single phyloseq object
```
otutable <- read.delim(file = "data/otutable.txt", sep = "\t", header = T, check.names = F, row.names = 1)
metadata <- read.delim(file = "data/metadata.txt", header = T, sep = "\t")
refseq <- readDNAStringSet("data/otus.fa", format = "fasta")
d <- amp_load(otutable = otutable, metadata = metadata, refseq = refseq)
```
