## Package

The `genbankr` package is part of the Bioconductor ecosystem of R packages. As
such, it's primary splash page - including continuous integration and testing
badges - can be found [here](https://bioconductor.org/packages/release/bioc/html/genbankr.html).


## Introduction

`genbankr` is a package which provides utilities to parse GenBank and GenPept 
files into data structures which integrate with the Bioconductor ecosystem.

## Installation

This package depends on the [Bioconductor](http://bioconductor.org) ecosystem of
R packages. If you already have a version of Bioconductor (>=3.3) installed,
you can do gthe following:

```
libary(BiocManager)
BiocManager::install("genbankr")
```

If you do not currently have the Bioconductor core machinery installed, you can
get the current or release version like so:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
``` 
for release and

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
```
For devel.

After doing one of the above (once), you can install `genbankr` as above, via

```
library(BiocManager)
BiocManager::install("genbankr")
```

Note that Bioconductor is a sychronized development and release
platform, so release and development versions of Bioconductor packages
cannot be safely mixed in the same package (installation) library. Use `switchr`
or direct `.libPaths` management to maintain multiple side-by-side installations
if necessary.

To install directly from github (this will generally not be necessary
unless you intend to contribute to `genbankr`'s development), do

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("gmbecker/genbankr")
```

The Bioconductor development repository will contain
the latest development version of genbankr which has passed testing (lagged
by about a day).

## Code of Conduct

This project operates under the Contributor Covenenant Code of Coduct [see here](./CONDUCT.md)

## Package Affiliations

The `genbankr` package is a part of the Bioconductor and rOpenSci projects.

| [![bioconductor_footer](http://bioconductor.org/images/logo_bioconductor.gif)](http://bioconductor.org) | [![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org) |
|:-------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------:|
