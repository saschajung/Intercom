# InterCom - Reconstruction of functional cell-cell communication networks 

Reconstructs functional cell-cell communication networks for a given single cell RNA-seq dataset.

## Required software
  - R v3.5 or greater
  - textshape v1.6.0
  - ggplot2 v 3.3.0
  - stuRpkg v1.3
  - dplyr v0.8.5
  - doParallel v1.0.15
  - stringr v1.4.0
  - plyr 1.8.6
  - igraph v1.2.5
  - Matrix v1.2-18
  - reshape2 v1.4.4
  - RSpectra v0.16-0
  - snow v0.4-3
  - taRifx v1.0.6.2
  - gtools v3.8.2
  - data.table v1.12.8
  - rlist v0.4.6.1

## Running InterCom
The complete InterCom workflow can be invoked through a single command after loading the library.
```R
InterCom(data,
         anno.tbl,
         species,
         sighot.cutoff=0.1,
         sighot.percentile=70,
         consv.thrs=0.05,
         ncores=4,
         sig.cutoff=0.9,
         z.score.cutoff=2,
         tissue.name,
         temp.folder.name,
         out.path
        )
```

## InterCom parameters
The main InterCom function comes with a variety of parameters for fine-tuning the output. However, in any case, the standard parameters work well in most of the cases. Below we describe every parameter, their ranges and standard values.

  - "data": Matrix or data frame of expression values; rows have to correspond to genes and columns to cells. Row names **must** be gene symbols while column names can be any identifier for a cell. Raw counts or normalized data is permitted.
  - "anno.tbl": Data frame annotating each cell with a cluster/cell type. Cell names have to be provided in the first and cluster/cell type information in the second column. 
  - "species": The organism from which the data was obtained. Currently, only "MOUSE" and "HUMAN" are supported.
  - "sighot.cutoff": Cutoff parameter for SigHotSpotter. Can be between 0 and 1.
  - "sighot.percentile": Percentile parameter for SigHotSpotter. Can be between 0 and 100.
  - "consv.thrs": Fraction of cells per cluster/cell type that must express a ligand, receptor or TF to be considered.
  - "ncores": Number of cores to use. Depending on the machine, can be any integer value greater or equal to 1.
  - "sig.cutoff": Significance cutoff between 0 (weakest) and 1 (strictest). Default: 0.9
  - "z.score.cutoff": Cutoff parameter to determine significant associations between receptors and interface TFs. Default: 2. 
  - "tissue.name": A name of the dataset.
  - "temp.folder.name": Name of the temporary folder to be created within the output path. Default: "temp" 
  - "out.path": Path to a folder where the output should be stored. Will be created if it does not exist.

## InterCom output
The main output of InterCom is a list object containing the final interactome and auxiliary information


