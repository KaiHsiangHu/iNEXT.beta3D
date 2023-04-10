<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.beta3D (R package)

<h5 align="right">
Latest version: 2023-04-10
</h5>
<font color="394CAE">
<h3 color="394CAE" style="font-weight: bold">
Introduction to iNEXT.beta3D (R package): Excerpt from iNEXT.beta3D
User’s Guide
</h3>
</font> <br>
<h5>
<b>Hsieh, Anne Chao, Kai-Hsiang Hu</b> <br><br> <i>Institute of
Statistics, National Tsing Hua University, Hsin-Chu, Taiwan 30043</i>
</h5>

<br> iNEXT.beta3D (iNterpolation and EXTrapolation for three dimensions
of Beta diversity) is an original R package available in
[Github](https://github.com/AnneChao). An online version of
[iNEXT.beta3D Online](https://chao.shinyapps.io/iNEXT_beta3D/) is also
available for users without an R background. Detailed information about
all functions in iNEXT.beta3D is provided in the iNEXT.beta3D Manual in
[iNEXT.beta3D_vignettes](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/A%20Quick%20Introduction%20to%20iNEXT.beta3D.html),
which is available from [Anne Chao’s
website](http://chao.stat.nthu.edu.tw/wordpress/software_download/).

`iNEXT.beta3D` focuses on Chiu’s et al.’s multiple assemblages
decomposition (Chiu, Jost, and Chao (2014)) based on three measures of
Hill numbers of order q: species richness (`q = 0`), Shannon diversity
(`q = 1`, the exponential of Shannon entropy) and Simpson diversity
(`q = 2`, the inverse of Simpson concentration). Besides, `iNEXT.beta3D`
also promote phylogenetic diversity level and functional diversity
level. For each diversity measure, `iNEXT.beta3D` uses the observed
sample of abundance or incidence data (called the “reference sample”) to
compute diversity estimates and the associated 95% confidence intervals
for the following two types of rarefaction and extrapolation (R/E):

1.  Sample-size-based (or size-based) R/E sampling curves:
    `iNEXT.beta3D` computes gamma, alpha, and beta diversity estimates
    for rarefied and extrapolated samples up to an appropriate size.
    This type of sampling curve plots the diversity estimates with
    respect to sample size.  
2.  Coverage-based R/E sampling curves: `iNEXT.beta3D` computes gamma,
    alpha, and beta diversity estimates for rarefied and extrapolated
    samples based on a standardized level of sample completeness (as
    measured by sample coverage) up to an appropriate coverage value.
    This type of sampling curve plots the diversity estimates with
    respect to sample coverage.

`iNEXT.beta3D` plots the above two types of sampling curves for gamma,
alpha, and beta diversity and four dissimilarity-turnover.

### SOFTWARE NEEDED TO RUN iNEXT.beta3D IN R

-   Required: [R](http://cran.rstudio.com/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### HOW TO RUN iNEXT.beta3D:

The `iNEXT.beta3D` package can be downloaded from Anne Chao’s
[iNEXT.beta3D_github](https://github.com/AnneChao/iNEXT.beta3D) using
the following commands. For a first-time installation, an additional
visualization extension package (`ggplot2`) and (`iNEXT.3D`) from Anne
Chao’s github must be installed and loaded.

``` r
# install_github('AnneChao/iNEXT.3D')
# library(iNEXT.3D)

## install the latest version from github
install.packages('devtools')
library(devtools)

install_github('AnneChao/iNEXT.beta3D')
## import packages
library(iNEXT.beta3D)

library(ggplot2)
```

Here are two main functions we provide in this package :

-   **iNEXTbeta3D** : Computing standardized gamma, alpha, beta
    diversity, and four dissimilarity-turnover indices for three
    dimensions: taxonomic, phylogenetic and functional diversity at
    specified sample coverage or sample size.

-   **ggiNEXTbeta3D** : Visualizing the output from the function
    `iNEXTbeta3D`

## MAIN FUNCTION: iNEXTbeta3D()

We first describe the main function `iNEXTbeta3D()` with default
arguments:

<br> iNEXTbeta3D(data, diversity = “TD”, q = c(0, 1, 2), datatype =
“abundance”, base = “coverage”, level = NULL, nboot = 20, conf = 0.95,
PDtree = NULL, PDreftime = NULL, PDtype = “meanPD”, FDdistM = NULL,
FDtype = “AUC”, FDtau = NULL, FDcut_number = 30) <br>

The arguments of this function are briefly described below, and will be
explained in more details by illustrative examples in later text. This
main function computes gamma, alpha and beta diversity estimates of
order q at specified sample coverage or sample size and measure of
diversity. By default of <code>base = “size”</code> and <code>level =
NULL</code>, then this function computes the gamma and alpha diversity
estimates up to double the reference sample size in each region. If
<code>base = “coverage”</code> and <code>level = NULL</code>, then this
function computes the gamma, alpha, beta diversity, and four
dissimilarity-turnover indices estimates up to one (for q = 1, 2) or up
to the coverage of double the reference sample size (for q = 0).

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style='text-align: left;'>

1.  For datatype = ‘abundance’, data can be input as a matrix/data.frame
    (species by assemblages), or a list of matrices/data.frames, each
    matrix represents species-by-assemblages abundance matrix.
2.  For datatype = ‘incidence_raw’ data can be input as a list (a
    region) with several lists (assemblages) of matrices/data.frames,
    each matrix represents species-by-sampling units.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    diversity
    </td>
    <td style="text-align: left;">
    selection of diversity type: ‘TD’ = Taxonomic diversity, ‘PD’ =
    Phylogenetic diversity, and ‘FD’ = Functional diversity.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    q
    </td>
    <td style="text-align: left;">
    a numerical vector specifying the diversity orders. Default is c(0,
    1, 2).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    datatype
    </td>
    <td style="text-align: left;">
    data type of input data: individual-based abundance data (datatype =
    ‘abundance’), sampling-unit-based incidence frequencies data
    (datatype = ‘incidence_freq’), or species by sampling-units
    incidence matrix (datatype = ‘incidence_raw’) with all entries being
    0 (non-detection) or 1 (detection).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    base
    </td>
    <td style="text-align: left;">
    Sample-sized-based rarefaction and extrapolation for gamma and alpha
    diversity (base = ‘size’) or coverage-based rarefaction and
    extrapolation for gamma, alpha and beta diversity (base =
    ‘coverage’). Default is base = ‘coverage’.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    level
    </td>
    <td style="text-align: left;">
    A numerical vector specifying the particular value of sample
    coverage (between 0 and 1 when base = ‘coverage’) or sample size
    (base = ‘size’). level = 1 (base = ‘coverage’) means complete
    coverage(the corresponding diversity represents asymptotic
    diversity). If base = ‘size’ and ‘level = NULL’, then this function
    computes the gamma and alpha diversity estimates up to double the
    reference sample size. If base = ‘coverage’ and ‘level = NULL’, then
    this function computes the gamma and alpha diversity estimates up to
    one (for ‘q = 1, 2’) or up to the coverage of double the reference
    sample size (for ‘q = 0’); the corresponding beta diversity is
    computed up to the same maximum coverage as the alpha diversity.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nboot
    </td>
    <td style="text-align: left;">
    a positive integer specifying the number of bootstrap replications
    when assessing sampling uncertainty and constructing confidence
    intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    conf
    </td>
    <td style="text-align: left;">
    a positive number \< 1 specifying the level of confidence interval.
    Default is ‘0.95’.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    PDtree
    </td>
    <td style="text-align: left;">
    (required only when diversity = ‘PD’), a phylogenetic tree in Newick
    format for all observed species in the pooled assemblage.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    PDreftime
    </td>
    <td style="text-align: left;">
    (required only when diversity = ‘PD’), a vector of numerical values
    specifying reference times for PD. Default is NULL (i.e., the age of
    the root of PDtree).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    PDtype
    </td>
    <td style="text-align: left;">
    (required only when diversity = ‘PD’), select PD type: PDtype = ‘PD’
    (effective total branch length) or PDtype = ‘meanPD’ (effective
    number of equally divergent lineages). Default is ‘meanPD’, where
    meanPD = PD/tree depth.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDdistM
    </td>
    <td style="text-align: left;">
    (required only when diversity = ‘FD’), select FD type: FDtype =
    ‘tau_values’ for FD under specified threshold values, or FDtype =
    ‘AUC’ (area under the curve of tau-profile) for an overall FD which
    integrates all threshold values between zero and one. Default is
    ‘AUC’.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDtype
    </td>
    <td style="text-align: left;">
    (required only when diversity = ‘FD’), select FD type: FDtype =
    ‘tau_values’ for FD under specified threshold values, or FDtype =
    ‘AUC’ (area under the curve of tau-profile) for an overall FD which
    integrates all threshold values between zero and one. Default is
    ‘AUC’.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDtau
    </td>
    <td style="text-align: left;">
    (required only when diversity = ‘FD’ and FDtype = ‘tau_values’), a
    numerical vector between 0 and 1 specifying tau values (threshold
    levels). If NULL (default), then threshold is set to be the mean
    distance between any two individuals randomly selected from the
    pooled assemblage (i.e., quadratic entropy).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    details
    </td>
    <td style="text-align: left;">
    a logical variable to decide whether do you want to print out the
    detailed value for each plots, default is `FALSE`.
    </td>
    </tr>
    <tr>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    FDcut_number
    </td>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    (required only when diversity = ‘FD’ and FDtype = ‘AUC’), a numeric
    number to split zero to one into several equal-spaced length.
    Default is ‘30’.
    </td>
    </tr>
    </tbody>
    </table>

This function returns an `"iNEXTbeta3D"` object which can be further
used to make plots using the function `ggiNEXTbeta3D()` to be described
below.

## DATA FORMAT/INFORMATION

Two types of data are supported:

1.  Individual-based abundance data (`datatype = "abundance"`): Input
    data for each region with several assemblages/sites include samples
    species abundances in an empirical sample of n individuals
    (“reference sample”). When there are N assemblages in a region,
    input data consist of a list with an S by N abundance matrix; For M
    regions consisting N assemblages, input data should be M lists of S
    by N abundance matrix.

2.  Sampling-unit-based incidence data (`datatype = "incidence_raw"`):
    Input data for a reference sample consist of a
    species-by-sampling-unit matrix; The number of sampling units and
    species for each assemblages in a region should be the same. For M
    regions consisting N assemblages, input data should be M lists of N
    lists of species-by-sampling-unit matrix.

Two data sets ( Beetles under different time periods for abundance and
incidence raw data ) are included in `iNEXT.beta3D` package for
illustration. For abundance data, the beetles data consist a list of two
matrices (“unlogged” and “logged”). Each matrix represents a
species-by-time (year 2008 and 2009) in a regions. For incidence raw
data, the beetles data consist two lists (“unlogged” and “logged”). For
each list, it has two matrices (each means year 2008 and 2009) which
represents a species-by-sampling units assemblage. For these data, the
following commands display how to compute taxonomic, phylogenetic and
functional beta decomposition and dissimilarities at specified sample
coverage for abundance and incidence raw data respectively.

Run the following code to view Beetles abundance data: (Here we only
show the first ten rows for each matrix)

``` r
data("beetle_abu")
beetle_abu
```

    $`2008 vs. 2009 unlogged`
                          2008 2009
    Abdera_flexuosa          0    2
    Onyxacalles_pyrenaeus    0    0
    Acanthocinus_griseus     2    0
    Agathidium_nigripenne    1    2
    Ampedus_aethiops        24    9
    Ampedus_auripes         90   79
    Ampedus_balteatus        2    5
    Ampedus_erythrogonus     3    0
    Ampedus_melanurus        0    0
    Ampedus_nigrinus        29   35

    $`2008 vs. 2009 logged`
                          2008 2009
    Abdera_flexuosa          0    0
    Onyxacalles_pyrenaeus    0    0
    Acanthocinus_griseus     0    0
    Agathidium_nigripenne    0    0
    Ampedus_aethiops        34   24
    Ampedus_auripes         88   40
    Ampedus_balteatus       19   11
    Ampedus_erythrogonus     0    2
    Ampedus_melanurus        0    1
    Ampedus_nigrinus       106   68

Run the following code to view Beetles incidence raw data: (Here we only
show the first ten rows and six columns for each matrix)

``` r
data("beetle_inc")
beetle_abu
```

    $`2008 vs. 2009 unlogged`
    $`2008 vs. 2009 unlogged`[[1]]
                          2008FKN_17 2008LAO_10 2008LAO_12 2008LAO_13 2008LAO_14 2008LAO_5
    Abdera_flexuosa                0          0          0          0          0         0
    Onyxacalles_pyrenaeus          0          0          0          0          0         0
    Acanthocinus_griseus           0          0          0          0          0         0
    Agathidium_nigripenne          0          0          0          0          0         0
    Ampedus_aethiops               0          1          1          1          0         0
    Ampedus_auripes                0          1          1          1          1         0
    Ampedus_balteatus              0          0          0          0          1         0
    Ampedus_erythrogonus           0          0          0          0          0         0
    Ampedus_melanurus              0          0          0          0          0         0
    Ampedus_nigrinus               0          1          1          0          1         0

    $`2008 vs. 2009 unlogged`[[2]]
                          2009FKN_17 2009LAO_10 2009LAO_12 2009LAO_13 2009LAO_14 2009LAO_5
    Abdera_flexuosa                0          0          0          0          0         0
    Onyxacalles_pyrenaeus          0          0          0          0          0         0
    Acanthocinus_griseus           0          0          0          0          0         0
    Agathidium_nigripenne          0          0          0          0          0         0
    Ampedus_aethiops               0          1          0          0          0         0
    Ampedus_auripes                0          1          1          1          1         0
    Ampedus_balteatus              0          1          0          0          0         0
    Ampedus_erythrogonus           0          0          0          0          0         0
    Ampedus_melanurus              0          0          0          0          0         0
    Ampedus_nigrinus               0          1          1          1          1         0


    $`2008 vs. 2009 logged`
    $`2008 vs. 2009 logged`[[1]]
                          2008FAE_1 2008FAE_11 2008FAE_16 2008FAE_18 2008FAE_2 2008FAE_21
    Abdera_flexuosa               0          0          0          0         0          0
    Onyxacalles_pyrenaeus         0          0          0          0         0          0
    Acanthocinus_griseus          0          0          0          0         0          0
    Agathidium_nigripenne         0          0          0          0         0          0
    Ampedus_aethiops              0          0          1          1         1          1
    Ampedus_auripes               0          0          1          1         1          1
    Ampedus_balteatus             0          1          1          1         0          1
    Ampedus_erythrogonus          0          0          0          0         0          0
    Ampedus_melanurus             0          0          0          0         0          0
    Ampedus_nigrinus              0          1          1          1         1          1

    $`2008 vs. 2009 logged`[[2]]
                          2009FAE_1 2009FAE_11 2009FAE_16 2009FAE_18 2009FAE_2 2009FAE_21
    Abdera_flexuosa               0          0          0          0         0          0
    Onyxacalles_pyrenaeus         0          0          0          0         0          0
    Acanthocinus_griseus          0          0          0          0         0          0
    Agathidium_nigripenne         0          0          0          0         0          0
    Ampedus_aethiops              0          0          1          1         0          1
    Ampedus_auripes               0          0          0          1         0          0
    Ampedus_balteatus             0          1          1          0         0          1
    Ampedus_erythrogonus          0          1          0          0         1          0
    Ampedus_melanurus             0          0          0          0         0          0
    Ampedus_nigrinus              0          1          1          1         1          1

When `base = 'coverage'`, the `iNEXTbeta3D()` function returns the
`"iNEXTbeta3D"` object including seven data frames for each regions:

-   gamma
-   alpha
-   beta
-   C ( Sorensen-type non-overlap )
-   U ( Jaccard-type non-overlap )
-   V ( Sorensen-type turnover )
-   S ( Jaccard-type turnover )

When `base = 'size'`, the `iNEXTbeta3D()` function returns the
`"iNEXTbeta3D"` object including two data frames for each regions:

-   gamma
-   alpha

## Rarefaction/Extrapolation Via Examples

## Taxonomic Diversity

First, we run the `iNEXTbeta3D()` function with beetles abundance data
to compute standardized gamma, alpha, beta diversity, and four
dissimilarities under `base = 'coverage'`. (Here we only show the first
six outputs for each output in region “unlogged”)

``` r
data("beetle_abu")

# Taxonomic diversity
Abundance_TD = iNEXTbeta3D(data = beetle_abu, diversity = 'TD', datatype = "abundance", 
                           base = 'coverage', nboot = 30)
Abundance_TD
```

    $gamma
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    2.032       0 Rarefaction 0.500 2.534 0.023 1.987 2.076 2008 vs. 2009 unlogged        TD
    2    1.909       1 Rarefaction 0.500 2.534 0.020 1.871 1.948 2008 vs. 2009 unlogged        TD
    3    1.776       2 Rarefaction 0.500 2.534 0.017 1.743 1.809 2008 vs. 2009 unlogged        TD
    4    2.174       0 Rarefaction 0.525 2.787 0.026 2.123 2.225 2008 vs. 2009 unlogged        TD
    5    2.028       1 Rarefaction 0.525 2.787 0.023 1.983 2.072 2008 vs. 2009 unlogged        TD
    6    1.871       2 Rarefaction 0.525 2.787 0.019 1.833 1.909 2008 vs. 2009 unlogged        TD

    $alpha
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    1.566       0 Rarefaction 0.500 4.141 0.022 1.522 1.610 2008 vs. 2009 unlogged        TD
    2    1.434       1 Rarefaction 0.500 4.141 0.019 1.397 1.471 2008 vs. 2009 unlogged        TD
    3    1.291       2 Rarefaction 0.500 4.141 0.016 1.260 1.322 2008 vs. 2009 unlogged        TD
    4    1.686       0 Rarefaction 0.525 4.615 0.026 1.636 1.737 2008 vs. 2009 unlogged        TD
    5    1.527       1 Rarefaction 0.525 4.615 0.022 1.484 1.570 2008 vs. 2009 unlogged        TD
    6    1.360       2 Rarefaction 0.525 4.615 0.018 1.325 1.395 2008 vs. 2009 unlogged        TD

    $beta
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    1.298       0 Rarefaction 0.500 4.141 0.005 1.288 1.307 2008 vs. 2009 unlogged        TD
    2    1.332       1 Rarefaction 0.500 4.141 0.005 1.322 1.341 2008 vs. 2009 unlogged        TD
    3    1.376       2 Rarefaction 0.500 4.141 0.005 1.367 1.385 2008 vs. 2009 unlogged        TD
    4    1.289       0 Rarefaction 0.525 4.615 0.005 1.279 1.299 2008 vs. 2009 unlogged        TD
    5    1.328       1 Rarefaction 0.525 4.615 0.005 1.318 1.338 2008 vs. 2009 unlogged        TD
    6    1.376       2 Rarefaction 0.525 4.615 0.005 1.366 1.385 2008 vs. 2009 unlogged        TD

    $C
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    0.298       0 Rarefaction 0.500 4.141 0.005 0.288 0.307 2008 vs. 2009 unlogged        TD
    2    0.413       1 Rarefaction 0.500 4.141 0.005 0.403 0.423 2008 vs. 2009 unlogged        TD
    3    0.546       2 Rarefaction 0.500 4.141 0.005 0.537 0.556 2008 vs. 2009 unlogged        TD
    4    0.289       0 Rarefaction 0.525 4.615 0.005 0.279 0.299 2008 vs. 2009 unlogged        TD
    5    0.409       1 Rarefaction 0.525 4.615 0.005 0.398 0.420 2008 vs. 2009 unlogged        TD
    6    0.546       2 Rarefaction 0.525 4.615 0.005 0.536 0.556 2008 vs. 2009 unlogged        TD

    $U
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    0.459       0 Rarefaction 0.500 4.141 0.006 0.447 0.470 2008 vs. 2009 unlogged        TD
    2    0.413       1 Rarefaction 0.500 4.141 0.005 0.403 0.423 2008 vs. 2009 unlogged        TD
    3    0.376       2 Rarefaction 0.500 4.141 0.005 0.367 0.385 2008 vs. 2009 unlogged        TD
    4    0.449       0 Rarefaction 0.525 4.615 0.006 0.437 0.461 2008 vs. 2009 unlogged        TD
    5    0.409       1 Rarefaction 0.525 4.615 0.005 0.398 0.420 2008 vs. 2009 unlogged        TD
    6    0.376       2 Rarefaction 0.525 4.615 0.005 0.366 0.385 2008 vs. 2009 unlogged        TD

    $V
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    0.298       0 Rarefaction 0.500 4.141 0.005 0.288 0.307 2008 vs. 2009 unlogged        TD
    2    0.332       1 Rarefaction 0.500 4.141 0.005 0.322 0.341 2008 vs. 2009 unlogged        TD
    3    0.376       2 Rarefaction 0.500 4.141 0.005 0.367 0.385 2008 vs. 2009 unlogged        TD
    4    0.289       0 Rarefaction 0.525 4.615 0.005 0.279 0.299 2008 vs. 2009 unlogged        TD
    5    0.328       1 Rarefaction 0.525 4.615 0.005 0.318 0.338 2008 vs. 2009 unlogged        TD
    6    0.376       2 Rarefaction 0.525 4.615 0.005 0.366 0.385 2008 vs. 2009 unlogged        TD

    $S
      Estimate Order.q      Method    SC  Size  s.e.   LCL   UCL                 Region diversity
    1    0.459       0 Rarefaction 0.500 4.141 0.006 0.447 0.470 2008 vs. 2009 unlogged        TD
    2    0.498       1 Rarefaction 0.500 4.141 0.005 0.487 0.509 2008 vs. 2009 unlogged        TD
    3    0.546       2 Rarefaction 0.500 4.141 0.005 0.537 0.556 2008 vs. 2009 unlogged        TD
    4    0.449       0 Rarefaction 0.525 4.615 0.006 0.437 0.461 2008 vs. 2009 unlogged        TD
    5    0.494       1 Rarefaction 0.525 4.615 0.006 0.483 0.505 2008 vs. 2009 unlogged        TD
    6    0.546       2 Rarefaction 0.525 4.615 0.005 0.536 0.556 2008 vs. 2009 unlogged        TD

The output contains seven data frames: `gamma`, `alpha`, `beta`, `C`,
`U`, `V`, `S`. For each data frame, it includes the diversity estimate
(`Estimate`), the diversity order (`Order.q`), `Method` (Interpolated,
Observed, or Extrapolated, depending on whether the size `m` is less
than, equal to, or greater than the reference sample size), the sample
coverage estimate (`SC`), the sample size (`Size`), the standard error
from bootstrap replications (`s.e.`), the 95% lower and upper confidence
limits of diversity (`LCL`, `UCL`), and the name of region (`Region`).
These diversity estimates with confidence intervals are used for
plotting the diversity curve.

Then, we can also use beetles incidence raw data to compute standardized
gamma, alpha, beta diversity, and four dissimilarities under
`base = 'coverage'`. (Here we don’t show the output)

``` r
data("beetle_inc")

Incidence_TD = iNEXTbeta3D(data = beetle_inc, diversity = 'TD', datatype = "incidence_raw",
                           base = 'coverage', nboot = 30)
Incidence_TD
```

## Phylogenetic Diversity

We can compute standardized gamma, alpha, beta diversity under
phylogenetic level with newick format data (`beetle_tree`). Here,
`PDtree` is the essential setting, and user can select `PDtype` and
`PDreft` alternatively. `PDtype` has two selection: `"PD"` (effective
total branch length) or `"meanPD"` (effective number of equally
divergent lineages). `PDreft` means a numerical value specifying
reference time for computing phylogenetic diversity. If `PDreft = NULL`,
then it is set the reference time to the age of the root of PDtree.
(Here we don’t show the output)

``` r
data("beetle_abu")
data("beetle_inc")
data("beetle_tree")


# Phylogenetic diversity
Abundance_PD = iNEXTbeta3D(data = beetle_abu, diversity = 'PD', datatype = "abundance", 
                           base = 'coverage', nboot = 30, PDtree = beetle_tree)
Abundance_PD

Incidence_PD = iNEXTbeta3D(data = beetle_inc, diversity = 'PD', datatype = "incidence_raw", 
                           base = 'coverage', nboot = 30, PDtree = beetle_tree)
Incidence_PD
```

## Functional Diversity

We can compute standardized gamma, alpha, beta diversity under
functional level with distance matrix (`beetle_distM`). Here, `FDdistM`
is the essential setting, and user can select `FDtype`, `FDtau`, and
`FDcut_number` alternatively. `FDtype` hase two selections: “tau_value”
means compute functional diversity under a specific threshold value and
“AUC” means compute functional diversity by integrating all threshold
values between zero and one. If user selects `FDtype = "tau_value"`, it
need to input a `FD_tau` threshold value. If `FD_tau = NULL`, then it is
set to be the mean distance between any two individuals randomly
selected from the pooled assemblage (i.e., quadratic entropy). (Here we
don’t show the output)

``` r
data("beetle_abu")
data("beetle_inc")
data("beetle_distM")


# Functional diversity - under specified threshold
Abundance_FDtau = iNEXTbeta3D(data = beetle_abu, diversity = 'FD', datatype = "abundance", 
                              base = 'coverage', nboot = 30, 
                              FDdistM = beetle_distM, FDtype = 'tau_value')
Abundance_FDtau

Incidence_FDtau = iNEXTbeta3D(data = beetle_inc, diversity = 'FD', datatype = "incidence_raw", 
                              base = 'coverage', nboot = 30, 
                              FDdistM = beetle_distM, FDtype = 'tau_value')
Incidence_FDtau


# Functional diversity - area under curve (integrates all threshold values between zero and one)
Abundance_FDAUC = iNEXTbeta3D(data = beetle_abu, diversity = 'FD', datatype = "abundance", 
                              base = 'coverage', nboot = 30, 
                              FDdistM = beetle_distM, FDtype = 'AUC')
Abundance_FDAUC

Incidence_FDAUC = iNEXTbeta3D(data = beetle_inc, diversity = 'FD', datatype = "incidence_raw", 
                              base = 'coverage', nboot = 30, 
                              FDdistM = beetle_distM, FDtype = 'AUC')
Incidence_FDAUC
```

### GRAPHIC DISPLAYS: FUNCTION ggiNEXTbeta3D()

The function `ggiNEXTbeta3D()`, which extends `ggplot2` to the
`"iNEXTbeta3D"` object with default arguments, is described as follows:

<br> ggiNEXTbeta3D(output, type = “B”, scale = “free”)  
<br>

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
output
</td>
<td style="text-align: left;">
</td>
</tr>
<tr>
<td style="text-align: left;">
type
</td>
<td style="text-align: left;">
(required only when base = ‘coverage’), selection of plot type : type =
‘B’ for plotting the gamma, alpha, and beta diversity ; type = ‘D’ for
plotting 4 turnover dissimilarities.
</td>
</tr>
<tr>
<td style="text-align: left;">
scale
</td>
<td style="text-align: left;">
Are scales shared across all facets (the default, ‘fixed’), or do they
vary across rows (‘free_x’), columns (‘free_y’), or both rows and
columns (‘free’)? Default is ‘free’.
</td>
</tr>
<tr>
<td style="border-bottom: 2px solid grey; text-align: left;">
transp
</td>
<td style="border-bottom: 2px solid grey; text-align: left;">
a value between 0 and 1 controlling transparency. ‘transp = 0’ is
completely transparent, default is ‘0.4’.
</td>
</tr>
</tbody>
</table>

The `ggiNEXTbeta3D()` function is a wrapper around the `ggplot2` package
to create a R/E curve using a single line of code. The resulting object
is of class `"ggplot"`, so it can be manipulated using the `ggplot2`
tools. Users can visualize the output of beta diversity or four
dissimilarities by setting the parameter ‘type’:

``` r
Abundance_TD = iNEXTbeta3D(data = beetle_abu, diversity = 'TD', datatype = "abundance", 
                           base = 'coverage', nboot = 30)
ggiNEXTbeta3D(Abundance_TD, type = 'B')
```

<img src="README/README-unnamed-chunk-13-1.png" width="1152" style="display: block; margin: auto;" />

``` r
ggiNEXTbeta3D(Abundance_TD, type = 'D')
```

<img src="README/README-unnamed-chunk-14-1.png" width="1152" style="display: block; margin: auto;" />

The following commands return the size-based R/E sampling curves for
gamma and alpha taxonomic diversity:

``` r
output = iNEXTbeta3D(data = beetle_abu, diversity = 'TD', datatype = "abundance", 
                     base = "size", nboot = 30)
ggiNEXTbeta3D(output)
```

<img src="README/README-unnamed-chunk-15-1.png" width="1152" style="display: block; margin: auto;" />

### How to cite

-   Chiu, C.-H., Jost, L. and Chao\*, A. (2014). Phylogenetic beta
    diversity, similarity, and differentiation measures based on Hill
    numbers. Ecological Monographs 84, 21-44.

### Referance

Chiu, C.-H., Jost, L. and Chao\*, A. (2014). Phylogenetic beta
diversity, similarity, and differentiation measures based on Hill
numbers. Ecological Monographs 84, 21-44.
