
<!-- README.md is generated from README.Rmd. Please edit that file -->
o2mod.transphylo: an outbreaker2 module using the TransPhylo model
------------------------------------------------------------------

*o2mod.transphylo* is a module of [*outbreaker2*](https://github.com/reconhub/outbreaker2) which uses the [*TransPhylo*](https://github.com/xavierdidelot/TransPhylo) model of within-host evolution.

Installation
------------

To install the latest version from github:

``` r
devtools::install_github("xavierdidelot/o2mod.transphylo")
```

Then load the package using:

``` r
library("o2mod.transphylo")
```

Documentation
-------------

Once loaded, the package provides a function `o2mod.transphylo()` which takes the same parameters as the `outbreaker()` function of the `outbreaker2` package, so for example if `data` is an outbreaker2 input object, analysis under the TransPhylo model can be carried out using:

``` r
result <- o2mod.transphylo(data)
```

For more details, see the [introductory vignette](https://github.com/xavierdidelot/o2mod.transphylo/blob/master/vignettes/introduction.md) on o2mod.transphylo.

Authors
-------

-   [Xavier Didelot](https://github.com/xavierdidelot)
-   [Finlay Campbell](https://github.com/finlaycampbell)
-   [Thibaut Jombart](https://github.com/thibautjombart)
