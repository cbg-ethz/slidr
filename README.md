# slidr <img src="https://github.com/cbg-ethz/slidr/blob/master/slidr.png" align="right" width="155 px"/>

An R package for identification of synthetic lethal partners for mutations from large perturbation screens.


## Data

The data used is big and cannot be stored on Github. The raw shRNA data has already been published as a part of project DRIVE (https://data.mendeley.com/datasets/y3ds55n88r/4 ) and all the mutation and copynumber data from CCLE is available at  https://portals.broadinstitute.org/ccle. The MutSig 2CV v3.1 MAF file for each cancer type is available at  http://firebrowse.org/. If you wish to use the processed data, please contact us and we'd be happy to share it. 


## SLIdR usage

You can install SLIdR using devtools.

```
install.packages("devtools") 
library(devtools) 
install_github("cbg-ethz/slidr")
```

The code used to process the data and run SLIdR is available in `Scripts/slidr_fin_run.Rmd`. `cellline_annot` used in the code corresponds to the supplementary Table S2 from the project DRIVE paper. `meta_data` should be a dataframe with variables:

* Bullet list
  * `Primary_site`: Name of primary site of tumor. Example: _pancreas_
  * `Driver_gene_file`: List of all MAF files for each cancer type. If the cancer has several sub-types then the MAF file names should be separated by `;`. Example for Lung cancer: _Rank_LUSC_MutSig2CV.txt ; Rank_LUAD_MutSig2CV.txt_ 
  * `Organ`: Name of the primary organ with the tumor. Can be the same as `Primary_site`.
  * `Additional_filters`: Names of tumor sub-types separated by `;`. 
  Example: _Lung:NSCLC_Adeno;Lung:NSCLC_Squamous_
  
## Author

The package and documentation is currently being updated. Therefore, if you have any questions, please contact <br/>
__Sumana Srivatsa__: [sumana.srivatsa@bsse.ethz.ch](sumana.srivatsa@bsse.ethz.ch)
  
