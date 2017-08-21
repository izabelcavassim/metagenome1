# LEfSe(linear discriminant analysis(LDA)),线性判别分析
```
Step1：利用Kruskal-Wallis秩和检验检测所有的特征物种，检测不同组间的物种丰度差异，并获得显著差异物种 

Step2：再利用Wilcoxon秩和检验检查在显著差异物种类中的所有亚种比较是否都趋同于同一分类级别 

Step3：对生成的向量集建立一个线性判别分析模型，最终得到特定分类水平下的差异物种列表及effect size 
```
**1.按影响性及关联性，以LEfSe算法计算各组中有显著差异的微生物群落或物种的LDA分值。**

LDA score大于预设值的显著差异物种，默认分值为2.0。柱状图的长短代表的是LDA score，即不同组间显著差异物种的影响程度。 

**2.根据LEfSe结果，以进化分枝树图展示各层次水平中存在组间差异的微生物群落或物种结构。** 

红色与绿色表示不同的分组情况，由内到外表示的是门、纲、目、科、属的物种分类水平。进化树中的红色节点为在红色组别中起重要作用的物生物分类，绿色节点为绿色组别中起到重要作用的物生物分类，黄色均代表无显著差异的物种。 

**3.根据LEfSe结果中检测到的存在组间差异的微生物群落或物种绘制丰度直方图。**

# 常见的biom格式转化

```
biom convert -i table.txt -o table.from_txt_json.biom \
--table-type="OTU table" --to-json
biom convert -i table.txt -o table.from_txt_hdf5.biom \
--table-type="OTU table" --to-hdf5
```
```
biom convert -i table.biom -o table.from_biom_w_taxonomy.txt --to-tsv \
--header-key taxonomy

biom convert -i table.biom -o table.from_biom_w_consensuslineage.txt --to-tsv \
--header-key taxonomy --output-metadata-id "Consensus Lineage"
```
## koeken模块
### 安装koeken模块
### 运行
```
koeken.py -i ./otu_table_even.biom -o ./output_folder/ -m ./lefse_map.txt -cl Description -sp Day
cp ./output_folder/lefse_output/format_lefse/2016_format.txt ./lefse.in
cp ./output_folder/lefse_output/run_lefse/2016.txt ./lefse.res
~/lefse/plot_res.py ./lefse.res ./lefse.pdf --dpi 600 --format pdf
~/lefse/plot_cladogram.py ./lefse.res ./lefse.cladogram.pdf --format pdf --dpi 600
```

## Galaxy Online
[Galaxy](http://huttenhower.sph.harvard.edu/galaxy/)

[Taxonomic biomarker discovery with LEfSe](https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial)

[lefse_wiki](https://bitbucket.org/biobakery/biobakery/wiki/lefse)

### 
```
python ~/metagenome1/Rscripts/LEfSe/format_input.py ./hmp_aerobiosis_small.txt ./hmp_aerobiosis_small.in \
-c 1 -s 2 -u 3 -o 1000000
```
 -c：大分组信息所在行;  
 -s：小分组信息所在行，如果没有小的分组可以不填;  
 -u：样品信息所在行;  
 -o：标准值，输入的丰度值按照该值重新计算，让输入的丰度值变大,，如果设置成负数，则不做任何处理;  

```
python ~/metagenome1/Rscripts/LEfSe/run_lefse.py ./hmp_aerobiosis_small.in ./hmp_aerobiosis_small.res
```
 -a：Kruskal-Wallis秩和检验筛选biomarker的p-value值;  
 -w：两组组间Wilcoxon秩和检验筛选biomarker的p-value值;  
 -l：LDA score
 --wilc：是否需要运行Wilcoxon step 0是运行，1是不运行，默认是运行Output;  

输出.res格式文件内容如下:  
Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae  5.0923016841  Low_O2  4.74694106197  2.91304680962e-07

总共5列，第一列biomarker名称;  
第二列是平均丰度最大的log10的值，如果平均丰度小于10的按照10来计算;  
第三列是差异基因或物种富集的组名称;  
第四列是LDA值;  
第五列是Kruskal-Wallis秩和检验的p值，如果不是biomarker则用“-”表示。

```
python ~/metagenome1/Rscripts/LEfSe/plot_res.py results/hmp_aerobiosis_small.res output_figures/hmp_aerobiosis_small.png
```
 --feature_font_size：设置feature字体的大小;  
 --format：图片输出的格式;  
 --dpi：图片的像素;  
 --title：标题名称，默认为空;  
 --title_font_size：标题字体大小;  
 --class_legend_font_size：图例字体大小;  
 --width：图片宽度;  
 --height：图片高度;
 --left_space：左边距;  
 --right_space：右边距
```
python ~/metagenome1/Rscripts/LEfSe/plot_cladogram.py ./hmp_aerobiosis_small.res ./hmp_aerobiosis_small.cladogram.png --format png
```
 --max_point_size：大点的大小，默认是6;  
 --min_point_size：小点的大小，默认是1;  
 --point_edge_width：圈的边线粗细，默认0.25;  
 --siblings_connector_width：同一级的宽度;  
 --parents_connector_width：上一级连接的宽度;  
 --title：标题;  
 --label_font_size：label字体大小;  
 --background_color：背景颜色
```
python ~/metagenome1/Rscripts/LEfSe/plot_features.py ./hmp_aerobiosis_small.in ./hmp_aerobiosis_small.res biomarkers_raw_images/
```
### 其他
```python
python2 ~/xiucz/metagenome1/biom2lefse.py \
-i ~/xiucz/metagenome1/caporaso-gut.biom \
-m ~/xiucz/metagenome1/caporaso-gut.txt \
-f SEX -o caporaso-lefse.txt
```

# Reference_Info
https://github.com/fw1121/pannenkoek
