<!-- README.md is generated from README.Rmd. Please edit that file -->


omplotr: 'ggplot2' Based RNAseq Plot Function Collection
========================================================

Theme
-----

`theme_onmath` is a ggplot theme used in almost all rnaseq plots.

``` r
library(omplotr)
```

    #> Loading required package: ggplot2

``` r
p <- ggplot(mtcars) + geom_point(aes(x = wt, y = mpg,colour = factor(gear)))
p + theme_onmath() + ggtitle("theme_onmath")
```

![](show/README-unnamed-chunk-2-1.png)

Plot
----

functions to generate plot in ngs analysis

### QC

#### Reads GC distribution

``` r
# Fastqc GC result
head(gc_test_data, 4)
```

    #>   X.Base     sample variable     value
    #> 1      1 DZ-A-LDY-1        A 0.1644290
    #> 2      2 DZ-A-LDY-1        A 0.2451466
    #> 3      3 DZ-A-LDY-1        A 0.2764013
    #> 4      4 DZ-A-LDY-1        A 0.3381226

``` r
# lineplot of GC distribution across Fastq file
gc_line_plot(gc_test_data)
```

![](show/README-unnamed-chunk-3-1.png)

#### Reads Quality barplot

``` r
# Reads Quality result
# Bars of Quality <= 30 were marked with color 'dodgerblue', 
# Bars of Quality > 30 were marked with color 'navy'.
head(rq_test_data, 4)
```

    #>   Quality Count   Proportion      color     sample
    #> 1      11    57 1.016564e-06 dodgerblue YP-B-WYX-6
    #> 2      12  7352 1.311189e-04 dodgerblue YP-B-WYX-6
    #> 3      13 40981 7.308736e-04 dodgerblue YP-B-WYX-6
    #> 4      14 57256 1.021129e-03 dodgerblue YP-B-WYX-6

``` r
# Reads Quality barplot
reads_quality_plot(rq_test_data)
```

![](show/README-unnamed-chunk-4-1.png)

### Quant

#### expression heatmap

``` r
# expression matrix
head(exp_test_data, 4)
```

    #>                 9dpi-pj-jp1 9dpi-pj-jp2 9dpi-pj-jp3 9dpi-pj-jp4
    #> ENSRNA049464904      64.515      48.860      34.595      25.636
    #> ENSRNA049468231   17763.048   28554.280    4878.607   12802.249
    #> ENSRNA049468277     544.106    1152.839     169.713     497.665
    #> ENSRNA049471043    4926.117    7815.150    1198.127    4545.421

``` r
# sample information
head(test_sample_data, 4)
```

    #>   condition      sample
    #> 1 9dpi-pj-A 9dpi-pj-jp1
    #> 2 9dpi-pj-A 9dpi-pj-jp2
    #> 3 9dpi-pj-B 9dpi-pj-jp3
    #> 4 9dpi-pj-B 9dpi-pj-jp4

``` r
# plot expression heatmap
om_heatmap(exp_test_data, test_sample_data)
```

![](show/README-unnamed-chunk-5-1.png)
