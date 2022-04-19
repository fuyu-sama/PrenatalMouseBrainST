# 计算方法
## 比对
### trimmomatic
使用trimmomatic将R1修剪，只保留前28bp。参数`trimmomatic-0.36.jar SE -threads 5 ${idx}_S1_L001_R1_001.fastq.gz ${idx}_S1_L001_R1_001.trimmed.fastq.gz CROP:28`。
### spaceranger
将修剪后的R1与完整的R2作为spaceranger输入，使用默认参数进行比对。

## 标准化、去批次  
### logcpm
对于m个spot s, n个基因g，logcpm计算如下：
$$ logcpm_{s, g}=\ln(\frac{count_{s,g}*10000}{\sum_{g=1}^{n}count_{s,g}}+1) $$
所有操作使用numpy（版本1.20.2）、pandas（版本1.2.4）进行。
### combat
将各个样品logcpm后的矩阵合并为一个完整的矩阵后，输入`scanpy.pp.combat()`函数（版本1.8.1），使用默认参数并保存输出，取赵金倩筛选的spatial pattern基因作为下游分析的输入。

## 聚类
### SC3聚类
使用R包SC3（版本1.18.0）函数`sc3()`，参数`biology = FALSE, gene_filter = FALSE`，将每个样本分别聚5-28类，挑选最好一组结果并手动划分区间。
### 层次聚类
首先计算所选SC3结果中所有类每个基因表达平均值，然后使用R函数`dist()`计算各个类之间的欧式距离并使用`hclust()`进行层次聚类。

## 差异表达
使用R包Seurat（版本4.0.2）函数`FindMarkers()`，默认参数（Wilcox秩和检验）进行。
### 区域特异性差异表达
将各个时间点相同区域的spot合并为同一类，指定比较对象为剩余所有spot运行函数。
### 时间点差异表达
对每个区域进行时间点差异表达。首先读取区域特异性差异表达结果，选取校正后p值小于0.01的基因进行时间点差异表达。将区域内待比较时间点的重复合并，指定比较对象为剩余所有spot运行函数。

## 降维
使用sklearn（版本0.24.1）进行。首先使用`PCA()`将combat结果降维至40维，然后使用`TSNE()`将PCA结果降维2维并作图。

## RCTD
将每个样本的原始count（原文要求）输入`RCTD`,读取对应单细胞数据（下丘脑：[13]，皮层：[14]），使用`create.RCTD()`及`run.RCTD()`函数默认参数进行RCTD.

## 转录因子富集
将各个区域上调的基因（avg_log2FC > 0， p_val_adj < 0.01）作为列表，输入homer软件（版本4.11）提供的`findMotifs.pl`脚本，参数`-start -1000 -end 1000`，进行转录因子富集分析。输出结果作转录因子气泡图、GO富集气泡图。  

# Cauchy Combination
$$ T_{ccr} = \Sigma_{i=1}^{n}w_{i}tan(\pi(0.5-p)) $$

# 使用的主要软件包版本
|   软件名称   |  版本  | 来源  |
| :---------: | :----: | :---: |
| spaceranger | 1.2.2  |  [1]  |
| trimmomatic |  0.36  |  [2]  |
|   python    | 3.8.5  |  [3]  |
| matplotlib  | 3.4.1  |  [4]  |
|    numpy    | 1.20.2 |  [5]  |
|   pandas    | 1.2.4  |  [6]  |
|   scanpy    | 1.8.1  |  [7]  |
|   sklearn   | 0.24.1 |  [8]  |
|      R      | 4.0.3  |  [9]  |
|    RCTD     | 1.2.0  | [10]  |
|     SC3     | 1.18.0 | [11]  |
|   Seurat    | 4.0.2  | [12]  |

[1] https://support.10xgenomics.com/spatial-gene-expression/software/overview/welcome  
[2] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.  
[3] https://www.python.org  
[4] J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.  
[5] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.   
[6] Jeff Reback, Wes McKinney, jbrockmendel, Joris Van den Bossche, Tom Augspurger, Phillip Cloud, gfyoung, Sinhrks, Adam Klein, Matthew Roeschke, Simon Hawkins, Jeff Tratner, Chang She, William Ayd, Terji Petersen, Marc Garcia, Jeremy Schendel, Andy Hayden, MomIsBestFriend, … Mortada Mehyar. (2020). pandas-dev/pandas: Pandas 1.0.3 (v1.0.3). Zenodo. https://doi.org/10.5281/zenodo.3715232  
[7] Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology, 19, 5, Article 15. https://doi.org/10.1186/s13059-017-1382-0  
[8] Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.  
[9] R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.  
[10] Cable, D. M., Murray, E., Zou, L. L. S., Goeva, A., Macosko, E. Z., Chen, F., & Irizarry, R. A.(2021). Robust decomposition of cell type mixtures in spatial transcriptomics [Article; Early Access]. Nature Biotechnology, 15. https://doi.org/10.1038/s41587-021-00830-w  
[11] Kiselev, V. Y., Kirschner, K., Schaub, M. T., Andrews, T., Yiu, A., Chandra, T., Natarajan, K. N., Reik, W., Barahona, M., Green, A. R., & Hemberg, M. (2017). SC3: consensus clustering of single-cell RNA-seq data. Nature Methods, 14(5), 483-+. https://doi.org/10.1038/nmeth.4236
[12]  
[13] Romanov, R. A., Tretiakov, E. O., Kastriti, M. E., Zupancic, M., Haring, M., Korchynska, S., Popadin, K., Benevento, M., Rebernik, P., Lallemend, F., Nishimori, K., Clotman, F., Andrews, W. D., Parnavelas, J. G., Farlik, M., Bock, C., Adameyko, I., Hokfelt, T., Keimpema, E., & Harkany, T. (2020). Molecular design of hypothalamus development. Nature, 82(7811), 246-+. https://doi.org/10.1038/s41586-020-2266-0  
[14] Di Bella, D. J., Habibi, E., Stickels, R. R., Scalia, G., Brown, J., Yadollahpour, P., Yang, S. M., Abbate, C., Biancalani, T., Macosko, E. Z., Chen, F., Regev, A., & Arlotta, P. (2021). Molecular logic of cellular diversification in the mouse cerebral cortex [Article]. Nature, 595(7868), 554-+. https://doi.org/10.1038/s41586-021-03670-5
