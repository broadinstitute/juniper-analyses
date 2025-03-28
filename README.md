# juniper-analyses
Code to run analyses and generate figures for JUNIPER: Reconstructing Transmission Events from Next-Generation Sequencing Data at Scale by Specht et al. (2025), available at https://doi.org/10.1101/2025.03.02.25323192.

## Setup

### Large File Merging

Before running code in this repository, please download the file juniper-analyses-oversize.zip from the Zenodo repository https://doi.org/10.5281/zenodo.14841603 and merge the contents with the juniper-analyses directory stored here. For example, the large file juniper-analyses-oversize/h5n1/res.RData from the Zenodo repository should be moved to juniper-analyses/h5n1/res.RData in this repository. The files deposited in the Zenodo repository were too large to host on this GitHub.

### Packages

The analyses in this repository require the following R packages: `juniper0`, `simulatR`, `outbreaker2`, `Rcpp`, `ape`, `ggplot2`, `igraph`, `readxl`, `reshape2`, `lubridate`, `coda`, and `cowplot`. `juniper0` and `simulatR` are custom R packages we built for this project; to install them, first install the `devtools` package, and then within R run `devtools::install_github("broadinstitute/juniper0")` and `devtools::install_github("broadinstitute/simulatR")`.

### BEAST 2

Our comparison to the `BadTrIP` outbreak reconstruction tool requires the `BEAST 2` software, which may be downloaded and installed from https://www.beast2.org/. We use `BEAST 2` version 2.7.5 for Mac OS X.

### IQTree and TreeTime

Our comparison to the `TransPhylo` outbreak reconstruction tool requires the `IQ-TREE` and `TreeTime` software packages, which may be downloaded and installed from http://www.iqtree.org and https://treetime.readthedocs.io/en/latest/installation.html. We use `IQ-TREE` version 2.3.6 and `TreeTime` version 0.11.4 for Mac OS X.

## De Novo iSNV Frequency Model

The code used to validate our model for the frequencies of *de novo* iSNVs can be found in `isnv_frequency_model.R`, and executed using
```
Rscript isnv_frequency_model.R
```
The estimate of the parameter *r* will be printed in the terminal, and a figure will be saved to the `figs` subdirectory.

## Validation of JUNIPER on Synthetic Data

### Data Generation and Model Execution

The directories `experiment_1` through `experiment_23`, which contain the synthetic outbreaks used to validate `JUNIPER` and reconstructions thereof using `JUNIPER`, `outbreaker2`, `BadTRIP`, and `TransPhylo`, can be generated using the script `synthetic_generation.R`. This script may be run from the command line using five arguments, the first being the index of the synthetic outbreak to generate and reconstruction, and the second through fifth being the number of global iterations at which to run `JUNIPER`, `outbreaker2`, `TransPhylo`, and `BadTrip`. Note that the same number of global iterations for different methods may have drastically different runtimes and effective sample sizes. On the synthetic datasets, to achieve reasonable effective sample sizes, we recommend a default of 10,000 iterations each for `JUNIPER` and `outbreaker2`, 100,000 for `BadTrIP`, and 1,000,000 for `TransPhylo`.

As an example, to generate the directory `experiment_1` with the above defaults, first navigate to the `juniper-analyses` directory. Then, from the command line, run
```
Rscript synthetic_generation.R 1 10000 10000 100000 1000000
```
Note that if the directory `experiment_1` already exists, the script will prompt the user to delete it and try again. Additionally, running different versions of `BEAST 2`, `IQ-TREE`, and `TreeTime` and/or using different operating systems may require the user to modify the calls to said software packages, located in the functions `compare_badtrip()` and `compare_outbreaker()` within `synthetic_generation.R`.

This GitHub reflects the `juniper-analyses` repository after the `synthetic_generation.R` script has already been run for all 23 synthetic datasets, so that the user may explore the results without needing to regenerate them (which is time consuming).

### Data Processing and Visualization

All performance statistics and figures for the synthetic data validation section of the `JUNIPER` manuscript were generated using the script `synthetic_analysis.R`. It can be run from the command line as follows:
```
Rscript synthetic_analysis.R
```
Summary statistics will be printed in the terminal, and figures will be saved to the `figs` subdirectory.

## Validation of JUNIPER on Real Data

The script `brsv-and-south-africa.R` runs `JUNIPER`, `outbreaker2`, `BadTRIP`, and `TransPhylo` on the Bovine RSV and South Africa SARS-CoV-2 datasets and analyzes the results. As regenerating the model runs may be time consuming, we saved the pre-computed JUNIPER and BadTrIP runs for the bovine RSV outbreak, and pre-computed runs of all models for the SARS-CoV-2 outbreaks. To run analyses and generate figures from these pre-computed model outputs, run
```
Rscript brsv-and-south-africa.R
```
To regenerate all `JUNIPER`, `outbreaker2`, and `TransPhylo` outputs, run
```
Rscript brsv-and-south-africa.R -r
```
Note that BadTrIP is heavily time-consuming on both datasets due to the presence of ambiguous sites. To avoid re-running BadTrIP, we reused earlier BadTrIP output on both datasets from a previous paper, titled "Inferring Viral Transmission Pathways from Within-Host Variation" (Specht et al. 2023). BadTrIP scrips from that project are available at `github.com/broadinstitute/transmission-inference`.

## Application of JUNIPER to Large-Scale H5N1 and SARS-CoV-2 Datasets

### Running JUNIPER on H5N1

We obtained our dataset of H5N1 genomes (FASTA and VCF files) from the repository `github.com/andersen-lab/avian-influenza/tree/master`, filtering to sequences collected from cattle with known dates of sample collection. Data were pulled in January 2025, and hence may not reflect the latest availability of genomes. As sequences in this repository were broken down by gene, we concatenated the relevant FASTA and VCFs for each case to obtain whole genomes. The results were deposited in the `h5n1/input_data` subdirectory of this repository.

Due to the large dataset size, we ran JUNIPER on a 32-CPU core Linux virtual machine via the following R script:
```
library(juniper0)
set.seed(1)
init <- initialize(
  n_global = 50000,
  init_mu = 2e-5
)
res <- run_mcmc(init)
save(res, file = "res.RData")
```
This script produced the file `h5n1/res.RData`, which was deposited in the `juniper-analyses-oversize` directory used for analysis.

### Running JUNIPER on SARS-CoV-2

We obtained our dataset of SARS-CoV-2 genomes from the NCBI SRA, BioProject PRJNA715749, filtering to samples from November 2021. FASTA and VCF files were generated using standard bioinformatics pipelines as described in *Methods* section of our manuscript, and were deposited in the `mass-10k/input_data` subdirectory of this repository.

Due to the even larger dataset size, we used a 224-CPU core Linux virtual machine to run JUNIPER and extract relevant statistics via the following R commands:
```
library(juniper0)
set.seed(1)
init <- initialize(
  n_global = 50000,
  init_mu = 2.25e-6,
  init_N_eff = 1
)
res <- run_mcmc(init, T)

# Additional post-processing steps, described below
```
After executing the `run_mcmc` function, we proceessed the output to generate a file `llik.RData` of the log-likelihood at each iteration, `n_kids.RData` of the number of offspring of each observed host at each iteration, and `ts.RData`, the times of infection of each host at each iteration. We also saved all details of a singular posterior sample as `snapshot.RData`, and recorded the console output as `log.out`. All of these files were placed in the `mass-10k` subdirectory.

### Data Processing and Visualization

All statistics and figures for the large-scale H5N1 and SARS-CoV-2 analyses were generated using the script `h5n1-and-mass.R`. As computing the number of outgoing transmissions per person per vaccination status for SARS-CoV-2 is time consuming, we saved these pre-computed statistics within the `mass-10k` subdirectory. To generate analyses and figures based on these pre-computed statistics, run
```
Rscript brsv-and-south-africa.R
```
To regenerate these statistics before proceeding to the analysis and figure generation, run
```
Rscript brsv-and-south-africa.R -r
```
Summary statistics will be printed in the terminal, and figures will be saved to the `figs` subdirectory.


