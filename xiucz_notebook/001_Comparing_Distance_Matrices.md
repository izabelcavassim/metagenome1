## Data
```
wget -c ftp://ftp.microbio.me/qiime/tutorial_files/88_soils.zip -O Comparing_Distance
```

## distance_matrix_from_mapping.py
**– Calculate the pairwise dissimilarity on one column of a mappping file**
```
cd ~/metagenome1/Qiime1/Comparing_Distance
distance_matrix_from_mapping.py -i map.txt -c PH
```
## compare_distance_matrices.py
### Mantel Test
```
compare_distance_matrices.py --method=mantel -i unweighted_unifrac_dm.txt,PH_dm.txt -o mantel_out -n 999
```

```
# Number of entries refers to the number of rows (or cols) retained in each
# distance matrix after filtering the distance matrices to include only those
# samples that were in both distance matrices. p-value contains the correct
# number of significant digits.
DM1	DM2	Number of entries	Mantel r statistic	p-value	Number of permutations	Tail type
unweighted_unifrac_dm.txt	PH_dm.txt	77	0.75592	0.001	999	two-sided
```
The Mantel r statistic of 0.75592 indicates that there is relatively strong positive correlation between the UniFrac and pH matrices. The p-value of 0.001 indicates that our results are statistically significant at an alpha of 0.05. We determined the p-value by specifying 999 permutations with the -n option. By default, the p-value is calculated using a two-tailed test, though this can be changed using the -t option.

### Partial Mantel Test
```
compare_distance_matrices.py --method=partial_mantel \
-i unweighted_unifrac_dm.txt,weighted_unifrac_dm.txt \
-c PH_dm.txt -o partial_mantel_out -n 999
```
```
# Number of entries refers to the number of rows (or cols) retained in each
# distance matrix after filtering the distance matrices to include only those
# samples that were in both distance matrices. p-value contains the correct
# number of significant digits.
DM1	DM2	CDM	Number of entries	Mantel r statistic	p-value	Number of permutations	Tail type
unweighted_unifrac_dm.txt	weighted_unifrac_dm.txt	PH_dm.txt	77	0.68183	0.001	999	greater
```
The Mantel r statistic of 0.68183 indicates that there is relatively strong positive correlation between the unweighted and weighted UniFrac distance matrices while controlling for differences in pH. The p-value of 0.001 indicates that our results are statistically significant at an alpha of 0.05. As with the Mantel test (above), we can also specify more than two distance matrices as inputs, and separate partial Mantel tests will be performed for all pairs of input distance matrices, using the same control matrix for each test.

### Mantel Correlogram

## R
```r
data(varespec)
data(varechem)
veg.dist <- vegdist(varespec) # Bray-Curtis
env.dist <- vegdist(scale(varechem), "euclid")
mantel(veg.dist, env.dist)
mantel(veg.dist, env.dist, method="spear")
```

## Reference_Info
http://qiime.org/tutorials/distance_matrix_comparison.html  
http://qiime.org/scripts/distance_matrix_from_mapping.html  
http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/mantel.html
