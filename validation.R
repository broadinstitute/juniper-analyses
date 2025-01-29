### Validation of JUNIPER's ability to reconstruct epidemics

library(juniper0)
library(simulatR)
library(ape)
library(reshape2)
library(outbreaker2)

args <- as.numeric(commandArgs(TRUE))
#print(args)

# Default params
defaults <- list(
  mu_g = 5,
  var_g = 5,
  mu_s = 5,
  var_s = 5,
  mu = 2e-5,
  N_eff = log(100),
  p_sample = 0.5,
  trans_sample = NA,
  R = 2,
  psi = 0.5,
  coverage = 1
)

params <- names(defaults)

# Alternative params
alternatives <- list(
  mu_g = c(2.5, 10),
  var_g = c(2.5, 10),
  mu_s = c(2.5, 10),
  var_s = c(2.5, 10),
  mu = c(1e-5, 4e-5),
  N_eff = c(log(20), log(500)),
  p_sample = c(0.25, 0.75),
  trans_sample = c(0.25, 0.75),
  R = c(1.5, 2.5),
  psi = c(0.25, 0.75),
  coverage = c(0.8, 0.9)
)

# Matrix of parameter inputs per experiment
combos <- as.data.frame(defaults)
for (p in params) {
  for (j in 1:length(alternatives[[p]])) {
    newrow <- combos[1, ]
    newrow[[p]] <- alternatives[[p]][j]
    combos <- rbind(combos, newrow)
  }
}

compare_outbreaker <- function(i, n_global){

  fasta <- read.FASTA(paste0(
    "experiment_", i, "/input_data/aligned.fasta"
  ))

  dates <- as.Date(gsub(".*\\|", "", names(fasta)))
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), "days"))

  names(fasta) <- NULL
  dat <- list(
    dna = fasta,
    dates = dates,
    w_dens = dgamma(1:20, shape = combos$mu_g[i]^2 / combos$var_g[i], rate = combos$mu_g[i] / combos$var_g[i]) # Needs to change for different generation ints...
  )

  out <- outbreaker(
    data = dat,
    config = list(n_iter = n_global, sample_every = 100, pb = T)
  )

  save(out, file = paste0("experiment_", i, "/outbreaker.RData"))
}

### Comparison to BadTrIP
compare_badtrip <- function(i, n_global){

  ## Write input
  fasta <- read.FASTA(paste0("experiment_", i, "/input_data/aligned.fasta"))

  names <- gsub("\\|.*", "", names(fasta))

  # Composition of nucleotides by site
  comp <- list()

  # Length of genome
  N_bases <- length(fasta[[1]])

  # Loop over hosts
  for (j in 1:length(fasta)) {
    comp[[j]] <- matrix(0, nrow = N_bases, ncol = 4)

    vcf <- read.delim(
      paste0("./experiment_", i, "/input_data/vcf/", names[j], ".vcf"),
      colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character")
    )

    for(k in 1:N_bases){
      if(fasta[[j]][k] %in% c(
        as.raw(136),
        as.raw(40),
        as.raw(72),
        as.raw(24)
      )){
        if(k %in% vcf$POS){
          ind <- match(k, vcf$POS)
          alleles <- c(vcf$REF[ind], vcf$ALT[ind])
          alleles <- match(alleles, c("A", "C", "G", "T"))
          af <- vcf$INFO[ind]
          af <- gsub(".*;AF=", "", af)
          af <- gsub(";SB=.*", "", af)
          af <- as.numeric(af)
          comp[[j]][k, alleles] <- round(c(1-af, af) * 10000)

        }else{

          allele <- match(fasta[[j]][k], c(
            as.raw(136),
            as.raw(40),
            as.raw(72),
            as.raw(24)
          ))
          comp[[j]][k, allele] <- 10000

        }
      }
    }
  }



  # For each position, map A, C, G, T to a random order of the those four nucleotides,
  # since BadTRiP behaves differently under different nucleotide prevalences
  # and our simulations treat the initial genome as all A's.
  shuffles <- list()
  for (k in 1:N_bases) {
    shuffles[[k]] <- sample(1:4, 4, replace = F)
  }

  for (j in 1:length(fasta)) {

    for (k in 1:N_bases) {
      comp[[j]][k, ] <- comp[[j]][k, shuffles[[k]]]
    }

    comp[[j]] <- paste(
      comp[[j]][, 1],
      comp[[j]][, 2],
      comp[[j]][, 3],
      comp[[j]][, 4],
      sep = "-"
    )

  }

  alignment <- as.data.frame(comp, col.names = paste0("S", 1:length(fasta)))

  if(!dir.exists(
    paste0("experiment_", i, "/badtrip")
  )){
    dir.create(paste0("experiment_", i, "/badtrip"))
  }

  write.table(alignment, row.names = F, quote = F, file = paste0("experiment_", i, "/badtrip/inputAlignment.txt"), sep = "\t")




  # Write dates

  dates <- as.Date(gsub(".*\\|", "", names(fasta)))
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), "days"))
  df <- data.frame(x = paste0("H", 1:length(fasta)), y = paste0("S", 1:length(fasta)), z = dates)
  write.table(df, file = paste0("experiment_", i, "/badtrip/inputSamples.txt"), quote= F, col.names = F, row.names = F, sep = "\t")

  epi <- data.frame(x = paste0("H", 1:length(fasta)), y = dates - combos$mu_s[i] - 2 * combos$mu_g[i], z = dates - combos$mu_s[i] + 2 * combos$mu_g[i])
  write.table(epi, file = paste0("experiment_", i, "/badtrip/inputEpiData.txt"), quote= F, col.names = F, row.names = F, sep = "\t")

  system(paste0(
    "python3 badtrip_scripts/create_BADTRIP_xml.py -a experiment_",
    i,
    "/badtrip/inputAlignment.txt -e experiment_",
    i,
    "/badtrip/inputEpiData.txt -s experiment_",
    i,
    "/badtrip/inputSamples.txt -m ",
    as.integer(n_global),
    " -o experiment_",
    i,
    "/badtrip/BADTRIP_setup"
  ))

  start <- Sys.time()
  save(start, file = paste0("experiment_", i, "/badtrip/start.RData"))
  # Run badtrip
  system(
    paste0(
      "/Applications/BEAST_2.7.5/bin/beast -threads 8 experiment_", i, "/badtrip/BADTRIP_setup.xml"
    )
  )
  end <- Sys.time()
  save(end, file = paste0("experiment_", i, "/badtrip/end.RData"))

  # Extract summary tree
  system(
    paste0(
      "python3 badtrip_scripts/Make_transmission_tree_alternative.py -i experiment_", i, "/badtrip/BADTRIP_setup.trees -b 20 -o experiment_", i, "/badtrip/summary"
    )
  )

}



# Generate a simulated epidemic using parameters in combos[i, ] and reconstruct it
validate <- function(i, n_global){

  # Create a new directory, if needed, for experiment i
  dname <- paste0("experiment_", i)

  # if(!dir.exists(dname)){
  #   dir.create(dname)
  # }
  #
  # # For reproducible results that vary by experiment
  # set.seed(i)
  #
  # # Simulate epidemic by rejection sampling
  # done <- F
  # while(!done){
  #   done <- epi_sim(
  #     a_g = combos$mu_g[i]^2 / combos$var_g[i],
  #     lambda_g = combos$mu_g[i] / combos$var_g[i],
  #     a_s = combos$mu_s[i]^2 / combos$var_s[i],
  #     lambda_s = combos$mu_s[i] / combos$var_s[i],
  #     R = combos$R[i],
  #     psi = combos$psi[i],
  #     mu = combos$mu[i],
  #     N_eff = combos$N_eff[i],
  #     init_genome = rep("A", 10000), # Initialize to genome of all A's
  #     p_samp = combos$p_sample[i],
  #     trans_samp = combos$trans_sample[i],
  #     coverage = combos$coverage[i],
  #     n_obs = 50,
  #     include_root = F,
  #     outdir = paste0("./", dname, "/input_data")
  #   )
  # }

  ## Run badtrip
  #set.seed(1)
  #compare_badtrip(i, n_global)

  ## Run outbreaker
  #set.seed(1)
  #compare_outbreaker(i, n_global)

  ## Reconstruct, first with perfectly-specified inputs
  # set.seed(1)
  # init <- initialize(
  #   n_subtrees = 1, # We will parallelize over simulations, not within each one
  #   n_global = n_global,
  #   indir = paste0(dname, "/input_data"),
  #   a_g = combos$mu_g[i]^2 / combos$var_g[i],
  #   lambda_g = combos$mu_g[i] / combos$var_g[i],
  #   a_s = combos$mu_s[i]^2 / combos$var_s[i],
  #   lambda_s = combos$mu_s[i] / combos$var_s[i],
  #   psi = combos$psi[i],
  #   init_mu = combos$mu[i],
  #   init_N_eff = combos$N_eff[i]
  # )
  #
  # res <- run_mcmc(init)
  # out <- summarize(res)
  #
  # # Save results for correctly-specified model
  # save(out, file = paste0("./", dname, "/output_ideal.RData"))
  #
  # ## Reconstruct, with misspecified (i.e. default inputs)
  # set.seed(1)
  # init <- initialize(
  #   n_subtrees = 1,
  #   n_global = n_global,
  #   indir = paste0(dname, "/input_data")
  # )
  #
  # res <- run_mcmc(init)
  # out <- summarize(res)
  #
  # # Save results for missipecified model
  # save(out, file = paste0("./", dname, "/output_misspecified.RData"))

  ## Reconstruct, next assuming all cases are sampled
  if(i == 3){
    set.seed(1)
    init <- initialize(
      n_subtrees = 1, # We will parallelize over simulations, not within each one
      n_global = n_global,
      indir = paste0(dname, "/input_data"),
      a_g = combos$mu_g[i]^2 / combos$var_g[i],
      lambda_g = combos$mu_g[i] / combos$var_g[i],
      a_s = combos$mu_s[i]^2 / combos$var_s[i],
      lambda_s = combos$mu_s[i] / combos$var_s[i],
      psi = combos$psi[i],
      init_mu = combos$mu[i],
      init_N_eff = combos$N_eff[i],
      init_pi = 1,
      ongoing = F,
      fixed_pi = T,
      rooted = T
    )

    res <- run_mcmc(init)
    out <- summarize(res)

    # Save results for model assuming perfect sequencing
    save(out, file = paste0("./", dname, "/output_fully_sampled.RData"))

    ## Reconstruct using ONLY consensus genomes
    unlink(paste0("./", dname, "/input_data/vcf"), recursive=TRUE)
  }

  set.seed(1)
  init <- initialize(
    n_subtrees = 1, # We will parallelize over simulations, not within each one
    n_global = n_global,
    indir = paste0(dname, "/input_data"),
    a_g = combos$mu_g[i]^2 / combos$var_g[i],
    lambda_g = combos$mu_g[i] / combos$var_g[i],
    a_s = combos$mu_s[i]^2 / combos$var_s[i],
    lambda_s = combos$mu_s[i] / combos$var_s[i],
    psi = combos$psi[i],
    init_mu = combos$mu[i],
    init_N_eff = combos$N_eff[i]
  )

  res <- run_mcmc(init)
  out <- summarize(res)

  # Save results for consensus JUNIPER
  save(out, file = paste0("./", dname, "/output_consensus.RData"))

}

validate(args[1], args[2])

# Create command
# cmd <- paste0("Rscript validation.R ", 1:12, " 100000 & ")
# writeLines(cmd, con = "badtrip_script.txt")



