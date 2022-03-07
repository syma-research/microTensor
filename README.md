# microTensor: Dimension Reduction for Longitudinal Microbiome Data

This repository contains the method implementation and analysis scripts used
for the paper "A Dimension Reduction Model For Longitudinal Microbiome Studies",
by Siyuan Ma and Hongzhe Li.
 
This is the main branch, which is used for the software and analysis used in the
paper! If you are interested in the R package simply implementing the microTensor
method, please switch to the "package" repo.

To use this repo to reproduce our paper:

1. Download the repository to a local directory.

2. In R, navigate your working directory to where the repo is stored.

3. Make sure you are under the main branch. In R, install the package with 
   `devtools::install()`
   
    * This installs the version of microTensor that is necessary to reproduce
      this analyses performed in the paper.
      
4. `rmd/` has the scripts for analyses performed in the paper. 

    * `0_Format_Data.Rmd` formats the abundance tables and metadata from the 
      FARMM and DIABIMMUNE cohorts, such that they are ready for downstreama 
      analyses. Since the formatted `.RData` files are also included in this
      repo, you can actually skip this file.
      
    * `1_Simulation_TrueMod.Rmd`, `2_Simulation_MisMod.Rmd`, and 
      `3_Simulation_SparseDOSSA.Rmd` are scripts for the three simulation 
      studies conducted in the paper. These contain large-scale simulation
      computation, and are meant to be run on a high-performance computing
      server. It uses the `batchtools` R package to schedule jobs, which
      you should configure according to your HPC workload manager.
      
    * `4_FARMM.Rmd` and `5_DIABIMMUNE.Rmd` are the scripts used to analyze
      data and visualize results for the two-real world studies. If you are
      interested in quick, reproducible analyses from the paper, these are
      the best targets, as they can be run on a local computer.

5. In addition to the packages necessary for the software implementation, 
   as specified in `DESCRIPTION`, the following R packages are also necessary
   for the analyses scripts:
   
    * [`tidyverse`](https://www.tidyverse.org/), for various data table 
      manipulation and visualization utilities.
    * As mentioned, [`batchtools`](https://cran.r-project.org/web/packages/batchtools/index.html)
      was used to submit and inspect jobs on HPC clusters for the simulation 
      studies.
    * [`SparseDOSSA`](https://huttenhower.sph.harvard.edu/sparsedossa2/) which
      is the simulation software used for simulation study 3.
    * [`lmerTest`](https://cran.r-project.org/web/packages/lmerTest/index.html)
      was used for testing on the correlated loadings identified by PCA.
    
    
Contact: [Siyuan Ma](syma.research@gmail.com)