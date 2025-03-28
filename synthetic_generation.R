### Validation of JUNIPER's ability to reconstruct epidemics
library(juniper0)
library(simulatR)
library(ape)
library(reshape2)
library(outbreaker2)
library(lubridate)
library(TransPhylo)

args <- as.numeric(commandArgs(TRUE))

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
    "synthetic-outbreaks/experiment_", i, "/input_data/aligned.fasta"
  ))

  dates <- as.Date(gsub(".*\\|", "", names(fasta)))
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), "days"))

  names(fasta) <- NULL
  dat <- list(
    dna = fasta,
    dates = dates,
    w_dens = dgamma(1:20, shape = combos$mu_g[i]^2 / combos$var_g[i], rate = combos$mu_g[i] / combos$var_g[i]) # Needs to change for different generation ints...
  )

  start <- Sys.time()
  out <- outbreaker(
    data = dat,
    config = list(n_iter = n_global, sample_every = 100, pb = T)
  )
  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/o2.RData"))

  save(out, file = paste0("synthetic-outbreaks/experiment_", i, "/outbreaker.RData"))
}

### Comparison to BadTrIP
compare_badtrip <- function(i, n_global){

  ## Write input
  fasta <- read.FASTA(paste0("synthetic-outbreaks/experiment_", i, "/input_data/aligned.fasta"))

  names <- gsub("\\|.*", "", names(fasta))

  # Composition of nucleotides by site
  comp <- list()

  # Length of genome
  N_bases <- length(fasta[[1]])

  # Loop over hosts
  for (j in 1:length(fasta)) {
    comp[[j]] <- matrix(0, nrow = N_bases, ncol = 4)

    vcf <- read.delim(
      paste0("./synthetic-outbreaks/experiment_", i, "/input_data/vcf/", names[j], ".vcf"),
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
    paste0("synthetic-outbreaks/experiment_", i, "/badtrip")
  )){
    dir.create(paste0("synthetic-outbreaks/experiment_", i, "/badtrip"))
  }

  write.table(alignment, row.names = F, quote = F, file = paste0("synthetic-outbreaks/experiment_", i, "/badtrip/inputAlignment.txt"), sep = "\t")

  # Write dates
  dates <- as.Date(gsub(".*\\|", "", names(fasta)))
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), "days"))
  df <- data.frame(x = paste0("H", 1:length(fasta)), y = paste0("S", 1:length(fasta)), z = dates)
  write.table(df, file = paste0("synthetic-outbreaks/experiment_", i, "/badtrip/inputSamples.txt"), quote= F, col.names = F, row.names = F, sep = "\t")

  epi <- data.frame(x = paste0("H", 1:length(fasta)), y = dates - combos$mu_s[i] - 2 * combos$mu_g[i], z = dates - combos$mu_s[i] + 2 * combos$mu_g[i])
  write.table(epi, file = paste0("synthetic-outbreaks/experiment_", i, "/badtrip/inputEpiData.txt"), quote= F, col.names = F, row.names = F, sep = "\t")

  system(paste0(
    "python3 badtrip_scripts/create_BADTRIP_xml.py -a synthetic-outbreaks/experiment_",
    i,
    "/badtrip/inputAlignment.txt -e synthetic-outbreaks/experiment_",
    i,
    "/badtrip/inputEpiData.txt -s synthetic-outbreaks/experiment_",
    i,
    "/badtrip/inputSamples.txt -m ",
    as.integer(n_global),
    " -o synthetic-outbreaks/experiment_",
    i,
    "/badtrip/BADTRIP_setup"
  ))

  start <- Sys.time()
  # Run badtrip
  system(
    paste0(
      "/Applications/BEAST_2.7.5/bin/beast -threads 8 synthetic-outbreaks/experiment_", i, "/badtrip/BADTRIP_setup.xml"
    )
  )
  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/bt.RData"))

  # Extract summary tree
  system(
    paste0(
      "python3 badtrip_scripts/Make_transmission_tree_alternative.py -i synthetic-outbreaks/experiment_", i, "/badtrip/BADTRIP_setup.trees -b 20 -o synthetic-outbreaks/experiment_", i, "/badtrip/summary"
    )
  )

}

### Comparison to IQTree + TransPhylo
compare_transphylo <- function(i, n_global){

  # Directory for iqtree files
  dir.create(
    paste0(
      "synthetic-outbreaks/experiment_", i, "/iqtree"
    )
  )

  # Run IQTree
  system(
    paste0(
      "iqtree2 -s synthetic-outbreaks/experiment_", i, "/input_data/aligned.fasta -m JC --prefix synthetic-outbreaks/experiment_", i, "/iqtree/iqtree-results"
    )
  )

  # Prepare date csv
  fasta <- read.FASTA(paste0("synthetic-outbreaks/experiment_", i, "/input_data/aligned.fasta"))
  names <- names(fasta)
  dates <- as.Date(gsub(".*\\|", "", names))

  # Convert dates to years, transphylo doesn't like small date spans in years
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01")))
  max_year <- 2000 + max(dates)
  dates <- as.Date("2000-01-01") + years(dates)

  date_csv <- data.frame(name = names, date = dates)
  write.csv(date_csv, file = paste0("synthetic-outbreaks/experiment_", i, "/iqtree/date.csv"), row.names = F, quote = F)

  # Run TreeTime
  system(
    paste0(
      "treetime --tree synthetic-outbreaks/experiment_",
      i,
      "/iqtree/iqtree-results.treefile --dates synthetic-outbreaks/experiment_",
      i,
      "/iqtree/date.csv --aln synthetic-outbreaks/experiment_",
      i,
      "/input_data/aligned.fasta --outdir synthetic-outbreaks/experiment_",
      i,
      "/treetime --coalescent skyline --n-skyline 2 --greedy-resolve"
    )
  )

  # Load timed tree
  phy <- read.nexus(paste0("synthetic-outbreaks/experiment_", i, "/treetime/timetree.nexus"))

  # Convert to binary tree
  phy <- multi2di(phy)

  # Make branch lengths nonzero
  phy$edge.length <- phy$edge.length + 0.01

  # Add max date
  ptree<-ptreeFromPhylo(phy,dateLastSample= max_year)
  #plot(ptree)

  start <- Sys.time()


  res_TransPhylo <- inferTTree(
    ptree,
    mcmcIterations=n_global,
    w.mean = combos$mu_g[i],
    w.std = sqrt(combos$var_g[i]),
    ws.mean = combos$mu_s[i],
    ws.std = sqrt(combos$var_s[i]),
    startOff.r = combos$R[i] * combos$psi[i] / (1 - combos$psi[i]),
    startOff.p = 1 - combos$psi[i],
    startPi = combos$p_sample[i],
    updateOff.r = T,
    updateOff.p = F,
    updatePi = T,
    dateT=max_year + 0.01
  )

  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/tp.RData"))

  # We typically run TransPhylo at a very high number of iterations yet achieve relatively low ESSs. Here slim down the posterior samples (without typically altering ESS)
  if(n_global > 1e5){
    res_TransPhylo <- res_TransPhylo[seq(1000, n_global, 1000)]
  }

  save(res_TransPhylo, file = paste0("synthetic-outbreaks/experiment_", i, "/res_TransPhylo_slim.RData"))

}



# Generate a simulated epidemic using parameters in combos[i, ] and reconstruct it using each method
validate <- function(
    i, # Which experiment number, corresponding to a row in the combos data frame
    n_global = 10000, # Number of global iterations for each Juniper method
    n_global_outbreaker = 10000, # Number of iterations for outbresker2
    n_global_transphylo = 1e6, # Number of iterations for TransPhylo (very slow convergence observed in practice)
    n_global_badtrip = 1e5 # Number of iterations for BadTrIP
){

  # Create new directory for synthetic outbreaks
  if(!dir.exists("synthetic-outbreaks")){
    dir.create("synthetic-outbreaks")
  }

  # Create a new directory, if needed, for experiment i
  dname <- paste0("synthetic-outbreaks/experiment_", i)


  if(!dir.exists(dname)){
    dir.create(dname)
  }else{
    stop(paste0("A directory called ", dname, " already exists. Please delete it and try again."))
  }

  dir.create(paste0(dname, "/runtime"))

  # For reproducible results that vary by experiment
  set.seed(i)

  # Simulate epidemic by rejection sampling
  done <- F
  while(!done){
    done <- epi_sim(
      a_g = combos$mu_g[i]^2 / combos$var_g[i],
      lambda_g = combos$mu_g[i] / combos$var_g[i],
      a_s = combos$mu_s[i]^2 / combos$var_s[i],
      lambda_s = combos$mu_s[i] / combos$var_s[i],
      R = combos$R[i],
      psi = combos$psi[i],
      mu = combos$mu[i],
      N_eff = combos$N_eff[i],
      init_genome = rep("A", 10000), # Initialize to genome of all A's
      p_samp = combos$p_sample[i],
      trans_samp = combos$trans_sample[i],
      coverage = combos$coverage[i],
      n_obs = 50,
      include_root = F,
      outdir = paste0("./", dname, "/input_data")
    )
  }

  ## Write metadata
  fasta <- ape::read.FASTA(paste0("./", dname, "/input_data/aligned.fasta"))
  samples <- gsub("\\|.*", "", names(fasta))
  dates <- gsub(".*\\|", "", names(fasta))
  meta <- data.frame(sample = samples, date = dates)
  big_meta <- data.frame(sample = c(samples, "ref_genome"), date = c(dates, "2000-01-01"))
  write.csv(meta, file = paste0("./", dname, "/input_data/metadata.csv"), row.names = F, quote = F)

  # Run badtrip
  set.seed(1)
  compare_badtrip(i, n_global_badtrip)

  # Run outbreaker
  set.seed(1)
  compare_outbreaker(i, n_global_outbreaker)

  # Run TransPhylo
  set.seed(1)
  compare_transphylo(i, n_global_transphylo)

  # Reconstruct with Juniper, first with perfectly-specified inputs
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

  start <- Sys.time()

  res <- run_mcmc(init)

  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/juniper_ideal.RData"))

  out <- summarize(res)

  # Save results for correctly-specified model
  save(out, file = paste0("./", dname, "/output_ideal.RData"))

  ## Reconstruct, with misspecified (i.e. default inputs)
  set.seed(1)
  init <- initialize(
    n_subtrees = 1,
    n_global = n_global,
    indir = paste0(dname, "/input_data")
  )

  start <- Sys.time()

  res <- run_mcmc(init)

  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  #save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/juniper_default.RData"))


  out <- summarize(res)

  # Save results for missipecified model
  save(out, file = paste0("./", dname, "/output_misspecified.RData"))

  ## Reconstruct, next assuming all cases are sampled
  ## In this case we have to supply the root of the cluster
  combined_fasta <- c(
    read.FASTA(paste0(dname, "/input_data/ref.fasta")),
    read.FASTA(paste0(dname, "/input_data/aligned.fasta"))
  )
  old_fasta <- read.FASTA(paste0(dname, "/input_data/aligned.fasta"))
  # Overwrite the old fasta
  write.FASTA(combined_fasta, file = paste0(dname, "/input_data/aligned.fasta"))
  # Overwrite metadata
  write.csv(big_meta, file = paste0("./", dname, "/input_data/metadata.csv"), row.names = F, quote = F)

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
    root = "ref_genome"
  )

  start <- Sys.time()

  res <- run_mcmc(init)

  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  #save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/juniper_fully_sampled.RData"))


  out <- summarize(res)

  # Save results for model assuming perfect sequencing
  save(out, file = paste0("./", dname, "/output_fully_sampled.RData"))

  ## Reconstruct using ONLY consensus genomes
  unlink(paste0("./", dname, "/input_data/vcf"), recursive=TRUE)

  # Overwrite the combined fasta to revert to what we had before
  write.FASTA(old_fasta, file = paste0(dname, "/input_data/aligned.fasta"))
  # Overwrite metadata
  write.csv(meta, file = paste0("./", dname, "/input_data/metadata.csv"), row.names = F, quote = F)

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

  start <- Sys.time()

  res <- run_mcmc(init)

  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "hours"))
  #save(runtime, file = paste0("synthetic-outbreaks/experiment_", i, "/runtime/juniper_consensus.RData"))


  out <- summarize(res)

  # Save results for consensus JUNIPER
  save(out, file = paste0("./", dname, "/output_consensus.RData"))

}

validate(args[1], args[2], args[3], args[4], args[5])



