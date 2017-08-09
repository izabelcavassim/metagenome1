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
