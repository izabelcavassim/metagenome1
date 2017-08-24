##

## 
```r
library(vegan)
data(dune)
data(dune.env)
```
```r
adonis(dune ~ Management*A1, data=dune.env, permutations=99)
Call:
adonis(formula = dune ~ Management * A1, data = dune.env, permutations = 99) 

Permutation: free
Number of permutations: 99

Terms added sequentially (first to last)

              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Management     3    1.4686 0.48953  3.2629 0.34161   0.01 **
A1             1    0.4409 0.44089  2.9387 0.10256   0.02 * 
Management:A1  3    0.5892 0.19639  1.3090 0.13705   0.25   
Residuals     12    1.8004 0.15003         0.41878          
Total         19    4.2990                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Management分组之间的 F-value= 3.2629，p-value = 0.01。

```r
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m) {
 library(vegan)
 co = as.matrix(combn(unique(factors),2))
 pairs = c()
 F.Model =c()
 R2 = c()
 p.value = c()

 for(elem in 1:ncol(co)){
 ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
 factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
 pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
 F.Model =c(F.Model,ad$aov.tab[1,4]);
 R2 = c(R2,ad$aov.tab[1,5]);
 p.value = c(p.value,ad$aov.tab[1,6])
 }
 p.adjusted =p.adjust(p.value,method=p.adjust.m)
 pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
 return(pairw.res)
}
pairwise.adonis(dune, dune.env$Management, sim.method="bray", p.adjust.m= "bonferroni")
     pairs  F.Model        R2 p.value p.adjusted
1 SF vs BF 2.514890 0.2643110   0.071      0.426
2 SF vs HF 1.857489 0.1710790   0.115      0.690
3 SF vs NM 3.425694 0.2551595   0.006      0.036
4 BF vs HF 1.567531 0.2071390   0.202      1.000
5 BF vs NM 2.715242 0.2794827   0.024      0.144
6 HF vs NM 3.423068 0.2755413   0.038      0.228
```

## 引用
http://meiweiping.cn/%E5%9C%A8R%E4%B8%AD%E6%AD%A3%E7%A1%AE%E8%BF%90%E8%A1%8CPERMANOVA-and-pairwise-comparison%E5%8F%8A%E6%B3%A8%E6%84%8F%E4%BA%8B%E9%A1%B9/

