
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RCTC

<!-- badges: start -->

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.4792447-0096FF.svg)](https://zenodo.org/badge/latestdoi/103871397)
[![GitHub
license](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/wltcwpf/RCTC/blob/master/LICENSE)
[![GitHub
release](https://img.shields.io/badge/Release-v1.0.0-blue.svg)](https://github.com/wltcwpf/RCTC/releases)
[![Report
Issues!](https://img.shields.io/badge/Report%20Issues-Here-1abc9c.svg)](https://github.com/wltcwpf/RCTC/issues)
[![Open Source?
Yes!](https://img.shields.io/badge/Open%20Source-Yes-green.svg)](https://github.com/wltcwpf/RCTC)
<!-- badges: end -->

The package was developed for Rotated Combination of Two-Component
(RCTC) Ground Motions calculation. For the details, please refer to the
[PEER
report](https://peer.berkeley.edu/sites/default/files/2017_09_stewart_9.10.18.pdf)

## Installation

1.  Make sure “Rcpp” and “pracma” are installed as they are
    dependencies. Otherwise, you can run the following command to
    install them.

``` r
install.packages(c('Rcpp', 'pracma'))
```

2.  For Mac users, you may need to download Xcode for command line tools
    or (by running the following command in Terminal):

<!-- -->

    xcode-select --install

For Windows users: you may need to download and install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) (and remember
to add path into the system environment path when installing):

3.  Install utility pacakge “devtools”:

``` r
install.packages("devtools")
```

4.  Lastly, install “RCTC” pacakge:

``` r
library(devtools)
install_github('wltcwpf/RCTC')
```

## Citation

If you use RCTC (directly or as a dependency of another package) for
work resulting in an academic publication, we would appreciate if the
report is cited:

    Wang, P., Stewart, J. P., Bozorgina, Y., Boore, D. M., Kishida, T. (2017). “R” Package for Computation of Earthquake Ground Motion Response Spectra. Report No. 2017/09. PEER, UC Berkeley.

## Usage

Please check out [wiki](https://github.com/wltcwpf/RCTC/wiki) for more
detailed instructions!
