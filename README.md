## Introduction

`genbankr` is a package which provides utilities to parse GenBank and GenPept 
files into data structures which integrate with the Bioconductor ecosystem.

## Installation

This package depends on the [Bioconductor](http://bioconductor.org) ecosystem of
R packages. To get the current release or development versions of the package, 
do

```
source("http://bioconductor.org/biocLite.R")
biocLite("genbankr")
``` 
for release and

```
source("http://bioconductor.org/biocLite.R")
try(useDevel(TRUE))
biocLite("genbankr")
```
For devel. Note that Bioconductor is a sychronized development and release
platform, so release and development versions of Bioconductor packages
cannot be safely mixed in the same package (installation) library. Use `switchr`
or direct `.libPaths` management to maintain multiple side-by-side installations
if necessary.

To install directly from github, do

```
source("http://bioconductor.org/biocLite.R")
try(useDevel(TRUE))
biocLite("gmbecker/genbankr")
```

Though this should generally not be necessary unless you are contributing to 
`genbankr`'s development. The Bioconductor development repository will contain
the latest development version of genbankr which has passed testing (lagged
by about a day).

NOTE: after the first time, sourcing `biocLite.R` can be replaced with 
`library(BiocInstaller)` and `useDevel` only need to be run once per 
package library.