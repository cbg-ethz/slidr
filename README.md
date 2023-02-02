![Build Status](https://travis-ci.org/cbg-ethz/slidr.svg?branch=master)
[![DOI](https://zenodo.org/badge/137060992.svg)](https://zenodo.org/badge/latestdoi/137060992)

# slidr <img src="https://github.com/cbg-ethz/slidr/blob/master/slidr.png" align="right" width="155 px"/>

An R package for identifying synthetic lethal pairs from large-scale perturbation screens.


### Data

The entire dataset used in the paper is big and cannot be stored on Github. The raw shRNA data has already been published as a part of project DRIVE (https://data.mendeley.com/datasets/y3ds55n88r/4 ) and all the mutation and copynumber data from CCLE are available at  https://portals.broadinstitute.org/ccle. The MutSig 2CV v3.1 MAF file for each cancer type is available at  http://firebrowse.org/. If you wish to use the processed data, please contact us and we'd be happy to share them. The code for processing the data, running both pan-cancer and cancer-type specific analyses, performing causal inference, and plotting the results are available under Scripts. 


### SLIdR usage

You can install SLIdR using `devtools`.

```
install.packages("devtools") 
library(devtools) 
install_github("cbg-ethz/slidr")
```
To run SLIdR, specify a path to store the results and use the `identifySLHits` function. An example dataset for liver cancer is available in the package under `LiverData`. 

```
library(slidr)
library(dplyr)

data(LiverData)

# Path for results
path_results <- "~/Downloads/"
# Threshold for significance in WT cell lines
thresh <-  0.1

hits <- slidr::identifySLHits(canc_data = LiverData, 
                              path_results = path_results, 
                              WT_pval_thresh = thresh)
                      
# Filtering significant hits in WT cell lines
hits <- hits %>% 
        dplyr::filter(WT_pvalue >= thresh)

```
The resulting variable `hits` is a dataframe of all the SL pairs in liver cancer as reported in the paper. SLIdR creates two folders in the specified output directory:
* `Hit_List` - including .txt file with all the hits before filtering by `WT_pval_thresh` value
* `Plots` - Boxplots of all the significant SL pairs after filtering by `WT_pval_thresh` value

### Contributions
[Sumana Srivatsa](https://bsse.ethz.ch/department/people/detail-person.MjAyOTQw.TGlzdC8yNjY5LDEwNjI4NTM0MDk=.html) <br/>
[Hesam Montazeri](http://lcbb.ut.ac.ir/)

### Contact

If you have any questions, please contact <br/>
[Sumana Srivatsa](mailto:sumana.srivatsa@bsse.ethz.ch)
  
