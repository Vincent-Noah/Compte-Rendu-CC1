<font size="6"> <span style="color:Red">*Script d’installation de tous
les packages.*</span> </font>
================

  - [<font size="4"> <span style="color:Red">*Mise à jour de la machine
    virtuelle.*</span> </font>](#mise-à-jour-de-la-machine-virtuelle.)
  - [<font size="4"> <span style="color:Red">*Installation de DADA2 et
    de Phyloseq avec BiocManager.*</span>
    </font>](#installation-de-dada2-et-de-phyloseq-avec-biocmanager.)
      - [<font size="3"> <span style="color:green">*Installation de
        BiocManager.*</span> </font>](#installation-de-biocmanager.)
      - [<font size="3"> <span style="color:green">*Installation de
        DADA2.*</span> </font>](#installation-de-dada2.)
      - [<font size="3"> <span style="color:green">*Installation de
        Phyloseq.*</span> </font>](#installation-de-phyloseq.)
  - [<font size="4"> <span style="color:Red">*Installation de packages
    pour phyloseq.*</span>
    </font>](#installation-de-packages-pour-phyloseq.)
      - [<font size="3"> <span style="color:green">*Installation de
        phangorn.*</span> </font>](#installation-de-phangorn.)
      - [<font size="3"> <span style="color:green">*Installation de
        DECIPHER.*</span> </font>](#installation-de-decipher.)
      - [<font size="3"> <span style="color:green">*Installation de
        gridExtra.*</span> </font>](#installation-de-gridextra.)
  - [<font size="4"> <span style="color:Red">*Installation de packages
    pour des analyses complémentaires sur phyloseq.*</span>
    </font>](#installation-de-packages-pour-des-analyses-complémentaires-sur-phyloseq.)

# <font size="4"> <span style="color:Red">*Mise à jour de la machine virtuelle.*</span> </font>

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
```

# <font size="4"> <span style="color:Red">*Installation de DADA2 et de Phyloseq avec BiocManager.*</span> </font>

A partir des instructions
<https://benjjneb.github.io/dada2/dada-installation.html>

## <font size="3"> <span style="color:green">*Installation de BiocManager.*</span> </font>

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.11')
```

## <font size="3"> <span style="color:green">*Installation de DADA2.*</span> </font>

``` r
BiocManager::install("dada2", version = "3.11")
```

## <font size="3"> <span style="color:green">*Installation de Phyloseq.*</span> </font>

``` r
BiocManager::install("phyloseq", version = "3.11")
```

# <font size="4"> <span style="color:Red">*Installation de packages pour phyloseq.*</span> </font>

## <font size="3"> <span style="color:green">*Installation de phangorn.*</span> </font>

``` r
BiocManager::install("phangorn")
```

## <font size="3"> <span style="color:green">*Installation de DECIPHER.*</span> </font>

``` r
BiocManager::install("DECIPHER")
```

## <font size="3"> <span style="color:green">*Installation de gridExtra.*</span> </font>

``` r
install.packages("gridExtra")
```

# <font size="4"> <span style="color:Red">*Installation de packages pour des analyses complémentaires sur phyloseq.*</span> </font>

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
```

``` r
install.packages(.cran_packages)
devtools::install_github(.github_packages)
BiocManager::install(.bioc_packages)                         
```
