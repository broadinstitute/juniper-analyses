## See how well other methods do reconstructing the same epidemics
library(ape)
library(outbreaker2)
library(lubridate)
library(TransPhylo)

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


### Comparison to IQTree + TransPhylo
compare_transphylo <- function(i, n_global){

  # Directory for iqtree files
  # dir.create(
  #   paste0(
  #     "experiment_", i, "/iqtree"
  #   )
  # )
  #
  # # Run IQTree
  # system(
  #   paste0(
  #     "iqtree2 -s experiment_", i, "/input_data/aligned.fasta -m JC --prefix experiment_", i, "/iqtree/iqtree-results"
  #   )
  # )

  # Prepare date csv
  fasta <- read.FASTA(paste0("experiment_", i, "/input_data/aligned.fasta"))
  names <- names(fasta)
  dates <- as.Date(gsub(".*\\|", "", names))

  # Convert dates to years, transphylo doesn't like small date spans in years
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01")))
  max_year <- 2000 + max(dates)
  dates <- as.Date("2000-01-01") + years(dates)

  date_csv <- data.frame(name = names, date = dates)
  #write.csv(date_csv, file = paste0("experiment_", i, "/iqtree/date.csv"), row.names = F, quote = F)

  # Run TreeTime
  # system(
  #   paste0(
  #     "treetime --tree experiment_",
  #     i,
  #     "/iqtree/iqtree-results.treefile --dates experiment_",
  #     i,
  #     "/iqtree/date.csv --aln experiment_",
  #     i,
  #     "/input_data/aligned.fasta --outdir experiment_",
  #     i,
  #     "/treetime --coalescent skyline --n-skyline 2 --greedy-resolve"
  #   )
  # )

  # Load timed tree
  phy <- read.nexus(paste0("experiment_", i, "/treetime/timetree.nexus"))

  # Convert to binary tree
  phy <- multi2di(phy)

  # Make branch lengths nonzero
  phy$edge.length <- phy$edge.length + 0.01

  # Add max date
  ptree<-ptreeFromPhylo(phy,dateLastSample= max_year)
  #plot(ptree)

  start <- Sys.time()
  save(start, file = paste0("experiment_", i, "/start_tp.RData"))

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
  save(end, file = paste0("experiment_", i, "/end_tp.RData"))

  save(res_TransPhylo, file = paste0("experiment_", i, "/res_TransPhylo.RData"))

}

### Comparison to BadTrIP
compare_badtrip <- function(i, n_global){

  ## Write input
  fasta <- read.FASTA(paste0("experiment_", i, "/input_data/aligned.fasta"))

  names <- gsub("\\|.*", "", names(fasta))

  # All positions with mutations
  all_pos <- integer(0)

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

          # A mutation occurred at this site
          all_pos <- c(all_pos, k)
        }else{

          allele <- match(fasta[[j]][k], c(
            as.raw(136),
            as.raw(40),
            as.raw(72),
            as.raw(24)
          ))
          comp[[j]][k, allele] <- 10000

          # If the allele is not A, a mutation occurred at this site
          if(allele != 1){
            all_pos <- c(all_pos, k)
          }
        }
      }
    }
  }

  all_pos <- sort(unique(all_pos))

  ## Extract positions with mutations for efficiency

  # Also, for each position, map A, C, G, T to a random order of the those four nucleotides,
  # since BadTRiP behaves differently under different nucleotide prevalences
  # and our simulations treat the initial genome as all A's.
  shuffles <- list()
  for (k in 1:length(all_pos)) {
    shuffles[[k]] <- sample(1:4, 4, replace = F)
  }

  for (j in 1:length(fasta)) {
    comp[[j]] <- comp[[j]][all_pos, ]

    for (k in 1:length(all_pos)) {
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


  system(
    paste0(
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
    )
  )

  # Run badtrip
  system(
    paste0(
      "/Applications/BEAST_2.7.7/bin/beast -threads 8 experiment_", i, "/badtrip/BADTRIP_setup.xml"
    )
  )


  # Extract summary tree
  system(
    paste0(
      "python3 badtrip_scripts/Make_transmission_tree_alternative.py -i experiment_", i, "/badtrip/BADTRIP_setup.trees -b 20 -o experiment_", i, "/badtrip/summary"
    )
  )

}

compare <- function(i, n_global, badtrip){
  if(badtrip == 1){
    compare_badtrip(i, n_global)
  }else{
    compare_transphylo(i, n_global)
  }
}

compare(args[1], args[2], args[3])
