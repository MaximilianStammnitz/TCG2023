![Scripting](https://img.shields.io/badge/Language-R-yellow.svg) ![Copyright](https://img.shields.io/badge/Copyright-(c)_2022_Max\_Stammnitz\_@TCG\_Cambridge-green.svg)

## Analysis scripts accompanying Stammnitz et al., 2023

_Maximilian Stammnitz, Transmissible Cancer Group, University of Cambridge (2015–2022)_

This repository contains custom R scripts which - in conjunction with the associated supplementary files - can be used to replicate the main figures presented in: **The evolution of two transmissible cancers in Tasmanian devils ([Stammnitz _et al._ 2023, Science 380:6642](https://doi.org/10.1126/science.abq6453))**.

The scripts are written to run on **[R](https://www.r-project.org/)** version 3.6.3 or later. You should be able to check your current version of R by running the command below:

```
R --version
```

The current R version is also shown when opening RStudio or the R Console.

Scripts require the following R packages: [**`readxl`**](https://cran.r-project.org/web/packages/readxl/index.html), [**`stringr`**](https://cran.r-project.org/web/packages/stringr/index.html), [**`scales`**](https://cran.r-project.org/web/packages/scales/index.html), [**`lubridate`**](https://cran.r-project.org/web/packages/lubridate/index.html), [**`data.table`**](https://cran.r-project.org/web/packages/data.table/index.html), [**`Biostrings`**](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [**`GenomicRanges`**](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [**`ggplot2`**](https://cran.r-project.org/web/packages/ggplot2/index.html), [**`ggmap`**](https://cran.r-project.org/web/packages/ggmap/index.html), [**`rstudioapi`**](https://cran.r-project.org/web/packages/rstudioapi/index.html), [**`ggsn`**](https://cran.r-project.org/web/packages/ggsn/index.html), [**`circlize`**](https://cran.r-project.org/web/packages/circlize/index.html), [**`caper`**](https://cran.r-project.org/web/packages/caper/index.html), [**`phytools`**](https://cran.r-project.org/web/packages/phytools/index.html), [**`treeio`**](https://bioconductor.org/packages/release/bioc/html/treeio.html), [**`ggtree`**](https://bioconductor.org/packages/release/bioc/html/ggtree.html), [**`rectanglePacking`**](https://github.com/TransmissibleCancerGroup/rectanglePacking)

Although care has been taken to make the code distribution-independent, it is possible that some of the scripts only work on Unix/MacOS systems, and may need to be modified in order to run on Windows systems.

For any feedback or requests, please get in touch with me via email: maxrupsta@gmail.com
