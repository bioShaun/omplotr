<!-- README.md is generated from README.Rmd. Please edit that file -->
``` r
knitr::opts_chunk$set(
  comment = "#>",
  fig.path = "show/README-",
  warning = FALSE,
  tidy = TRUE, 
  tidy.opts = list(comment = FALSE)
)
```

omplotr: 'ggplot2' Based RNAseq Plot Function Collection
========================================================

Theme
-----

`theme_onmath` is a ggplot theme used in almost all rnaseq plots.

``` r
library(omplotr)
```

    #> Loading required package: ggplot2

    #> Loading required package: topGO

    #> Loading required package: BiocGenerics

    #> Loading required package: parallel

    #> 
    #> Attaching package: 'BiocGenerics'

    #> The following objects are masked from 'package:parallel':
    #> 
    #>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    #>     clusterExport, clusterMap, parApply, parCapply, parLapply,
    #>     parLapplyLB, parRapply, parSapply, parSapplyLB

    #> The following objects are masked from 'package:stats':
    #> 
    #>     IQR, mad, xtabs

    #> The following objects are masked from 'package:base':
    #> 
    #>     anyDuplicated, append, as.data.frame, cbind, colnames,
    #>     do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    #>     grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
    #>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    #>     Position, rank, rbind, Reduce, rownames, sapply, setdiff,
    #>     sort, table, tapply, union, unique, unsplit, which, which.max,
    #>     which.min

    #> Loading required package: graph

    #> Loading required package: Biobase

    #> Welcome to Bioconductor
    #> 
    #>     Vignettes contain introductory material; view with
    #>     'browseVignettes()'. To cite Bioconductor, see
    #>     'citation("Biobase")', and for packages 'citation("pkgname")'.

    #> Loading required package: GO.db

    #> Loading required package: AnnotationDbi

    #> Loading required package: stats4

    #> Loading required package: IRanges

    #> Loading required package: S4Vectors

    #> 
    #> Attaching package: 'S4Vectors'

    #> The following objects are masked from 'package:base':
    #> 
    #>     colMeans, colSums, expand.grid, rowMeans, rowSums

    #> 

    #> Loading required package: SparseM

    #> 
    #> Attaching package: 'SparseM'

    #> The following object is masked from 'package:base':
    #> 
    #>     backsolve

    #> 
    #> groupGOTerms:    GOBPTerm, GOMFTerm, GOCCTerm environments built.

    #> 
    #> Attaching package: 'topGO'

    #> The following object is masked from 'package:IRanges':
    #> 
    #>     members

``` r
p <- ggplot(mtcars) + geom_point(aes(x = wt, y = mpg, colour = factor(gear)))
p + theme_onmath() + ggtitle("theme_onmath")
```

![](show/README-unnamed-chunk-2-1.png)

Plot
----

functions to generate plot in ngs analysis

### QC

#### Reads GC distribution

``` r
head(gc_test_data, 4)
```

    #>   X.Base     sample variable     value
    #> 1      1 DZ-A-LDY-1        A 0.1644290
    #> 2      2 DZ-A-LDY-1        A 0.2451466
    #> 3      3 DZ-A-LDY-1        A 0.2764013
    #> 4      4 DZ-A-LDY-1        A 0.3381226

``` r
gc_line_plot(gc_test_data)
```

![](show/README-unnamed-chunk-3-1.png)

#### Reads Quality barplot

``` r
head(rq_test_data, 4)
```

    #>   Quality Count   Proportion      color     sample
    #> 1      11    57 1.016564e-06 dodgerblue YP-B-WYX-6
    #> 2      12  7352 1.311189e-04 dodgerblue YP-B-WYX-6
    #> 3      13 40981 7.308736e-04 dodgerblue YP-B-WYX-6
    #> 4      14 57256 1.021129e-03 dodgerblue YP-B-WYX-6

``` r
reads_quality_plot(rq_test_data)
```

![](show/README-unnamed-chunk-4-1.png)

### Quant

#### expression box, violin and density plot

``` r
head(exp_test_data, 4)
```

    #>                 9dpi-pj-jp1 9dpi-pj-jp2 9dpi-pj-jp3 9dpi-pj-jp4
    #> ENSRNA049464904      64.515      48.860      34.595      25.636
    #> ENSRNA049468231   17763.048   28554.280    4878.607   12802.249
    #> ENSRNA049468277     544.106    1152.839     169.713     497.665
    #> ENSRNA049471043    4926.117    7815.150    1198.127    4545.421

``` r
head(test_sample_data, 4)
```

    #>   condition      sample
    #> 1 9dpi-pj-A 9dpi-pj-jp1
    #> 2 9dpi-pj-A 9dpi-pj-jp2
    #> 3 9dpi-pj-B 9dpi-pj-jp3
    #> 4 9dpi-pj-B 9dpi-pj-jp4

``` r
om_boxplot(exp_test_data, test_sample_data, "box")
```

![](show/README-unnamed-chunk-5-1.png)

``` r
om_boxplot(exp_test_data, test_sample_data, "violin")
```

![](show/README-unnamed-chunk-5-2.png)

``` r
om_boxplot(exp_test_data, test_sample_data, "density")
```

![](show/README-unnamed-chunk-5-3.png)

``` r
om_boxplot(exp_test_data, test_sample_data, "all")
```

![](show/README-unnamed-chunk-6-1.png)

#### expression PCA analysis point plot

``` r
om_pca_plot(exp_test_data, test_sample_data)
```

![](show/README-unnamed-chunk-7-1.png)

### expression correlation heatmap

``` r
om_correlation_plot(exp_test_data, test_sample_data)
```

![](show/README-unnamed-chunk-8-1.png)

#### diff expression volcano plot

``` r
head(diff_test_data, 4)
```

    #>                  Gene_ID X9dpi.pj.jp1 X9dpi.pj.jp2 X9dpi.pj.jp3
    #> 21437       Os01g0977250        4.780        4.551        7.806
    #> 33806       Os03g0738600        1.209        0.946        1.338
    #> 30663       Os03g0823900        5.598        5.303        6.322
    #> 32700 EPlOSAG00000051674        0.000        0.000        1.482
    #>       X9dpi.pj.jp4       logFC    PValue       FDR          compare
    #> 21437        6.489 -0.21257899 0.3964249 0.6631611 Case1_vs_Control
    #> 33806        1.392  0.03988567 0.9082303 1.0000000 Case1_vs_Control
    #> 30663        7.596  0.04831892 0.8428859 0.9858939 Case1_vs_Control
    #> 32700        0.000 -2.62518576 1.0000000 1.0000000 Case1_vs_Control

``` r
om_volcano_plot(diff_test_data, "Case_vs_Control")
```

![](show/README-unnamed-chunk-9-1.png)

``` r
om_volcano_plot(diff_test_data, "ALL")
```

![](show/README-unnamed-chunk-9-2.png)

#### expression heatmap

``` r
om_heatmap(exp_test_data, test_sample_data)
```

![](show/README-unnamed-chunk-10-1.png)

#### expression cluster line plot

``` r
head(cluster_test_data, 4)
```

    #>     cluster            Gene_id    variable     value
    #> 1 cluster_1 EPlOSAG00000008604 9dpi-pj-jp1 0.4935853
    #> 2 cluster_1 EPlOSAG00000018483 9dpi-pj-jp1 0.5236810
    #> 3 cluster_1 EPlOSAG00000027962 9dpi-pj-jp1 0.5008437
    #> 4 cluster_1 EPlOSAG00000045132 9dpi-pj-jp1 0.4879498

``` r
om_cluster_plot(cluster_test_data)
```

![](show/README-unnamed-chunk-11-1.png)

### Enrichment analysis

#### GO

``` r
head(go_test_data_list[["test_diff_genes"]], 4)
```

    #> [1] "Os01g0150000" "Os01g0172600" "Os01g0220700" "Os01g0303600"

``` r
head(go_test_data_list[["test_gene_len"]], 4)
```

    #>        gene_id gene_len
    #> 1 Os01g0290700     2289
    #> 2 Os01g0249200     5219
    #> 3 Os01g0152200     3747
    #> 4 Os01g0295700     5836

``` r
head(go_test_data_list[["test_go_anno"]], 4)
```

    #>        gene_id      go_id
    #> 1 Os01g0100100 GO:0005622
    #> 2 Os01g0100100 GO:0006886
    #> 3 Os01g0100100 GO:0017137
    #> 4 Os01g0100100 GO:0005096

``` r
goseq_output <- om_goseq(go_test_data_list[["test_diff_genes"]], go_test_data_list[["test_gene_len"]], 
    go_test_data_list[["test_go_anno"]])
```

    #> Using manually entered categories.

    #> Calculating the p-values...

    #> 'select()' returned 1:1 mapping between keys and columns

![](show/README-unnamed-chunk-12-1.png)

``` r
head(goseq_output, 4)
```

    #>       category over_represented_pvalue    qvalue numDEInCat numInCat
    #> 571 GO:0010287             0.008881002 0.3812842          2        2
    #> 848 GO:0042803             0.008933210 0.3812842          2        2
    #> 76  GO:0001172             0.024128606 0.3812842          2        3
    #> 123 GO:0003968             0.024128606 0.3812842          2        3
    #>                                     term ontology
    #> 571                        plastoglobule       CC
    #> 848    protein homodimerization activity       MF
    #> 76          transcription, RNA-templated       BP
    #> 123 RNA-directed RNA polymerase activity       MF
    #>                         DE_id
    #> 571 Os01g0118000,Os01g0173000
    #> 848 Os01g0235200,Os01g0251000
    #> 76  Os01g0198000,Os01g0198100
    #> 123 Os01g0198000,Os01g0198100

``` r
gene_go_map <- system.file("extdata", "topgo_test_data.txt", package = "omplotr")
run_topgo(gene_go_map, go_test_data_list[["test_diff_genes"]], goseq_output)
```

    #> 
    #> Building most specific GOs .....

    #>  ( 384 GO terms found. )

    #> 
    #> Build GO DAG topology ..........

    #>  ( 648 GO terms and 805 relations. )

    #> 
    #> Annotating nodes ...............

    #>  ( 779 genes annotated to the GO terms. )

    #> Loading required package: Rgraphviz

    #> Loading required package: grid

    #> 
    #> Attaching package: 'grid'

    #> The following object is masked from 'package:topGO':
    #> 
    #>     depth

    #> 
    #> Attaching package: 'Rgraphviz'

    #> The following objects are masked from 'package:IRanges':
    #> 
    #>     from, to

    #> The following objects are masked from 'package:S4Vectors':
    #> 
    #>     from, to

    #> 
    #> Building most specific GOs .....

    #>  ( 164 GO terms found. )

    #> 
    #> Build GO DAG topology ..........

    #>  ( 331 GO terms and 679 relations. )

    #> 
    #> Annotating nodes ...............

    #>  ( 778 genes annotated to the GO terms. )

![](show/README-unnamed-chunk-12-2.png)

    #> 
    #> Building most specific GOs .....

    #>  ( 517 GO terms found. )

    #> 
    #> Build GO DAG topology ..........

    #>  ( 1566 GO terms and 3214 relations. )

    #> 
    #> Annotating nodes ...............

    #>  ( 726 genes annotated to the GO terms. )

![](show/README-unnamed-chunk-12-3.png)![](show/README-unnamed-chunk-12-4.png)
