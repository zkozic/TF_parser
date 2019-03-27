### Usage:

`Rscript TF_parser.R path_to_file`

### Note:
1. File names can not have the $ sign in their names, it has to be removed from the filename. Bash is unable to properly use it as an argument for the R script.
2. 16th line in the TF_parser.R has to be hardcoded with a path to the Python binary which has pandas library installed. If deleted or commented out, it will default to `/usr/bin/python`, which is fine if this binary has the pandas library.
3. .xlsx file will be output in the directory from which the script has been run from.

### Required libraries:
#### Python
* Bio  
* pandas  

#### R
* tidyverse
* stringr
* GenomicRanges (Bioconductor)
* WriteXLS
* reticulate
