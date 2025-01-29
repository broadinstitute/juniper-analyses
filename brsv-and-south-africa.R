### Script to generate figures for BRSV and South Africa COVID-19 clusters

library(juniper0)
library(ggplot2)
library(cowplot)
library(readxl)
library(outbreaker2)
library(TransPhylo)
library(lubridate)
library(igraph)

### Bovine RSV Cluster

# Set to TRUE to re-generate model runs; otherwise use saved model runs
regenerate <- FALSE

if(regenerate){
  # Run JUNIPER
  set.seed(0)
  init <- initialize(
    indir = "brsv-outbreak/input_data/",
    n_subtrees = 1,
    init_pi = 1,
    init_mu = 3e-6,
    n_global = 1000,
    #split_bottlenecks = T,
    rooted = T,
    a_g = 7,
    ongoing = F #,
    #fixed_pi = T
  )

  res <- run_mcmc(init, T)
  save(res, file = "brsv-outbreak/output.RData")
  # Run JUNIPER with split bottlenecks
  set.seed(0)
  init <- initialize(
    indir = "brsv-outbreak/input_data/",
    n_subtrees = 1,
    #init_pi = 1,
    init_mu = 3e-6,
    n_global = 1000,
    split_bottlenecks = T,
    rooted = T,
    a_g = 7,
    ongoing = F #,
    #fixed_pi = T
  )

  res_unsampled <- run_mcmc(init, T)
  save(res_unsampled, file = "brsv-outbreak/output_unsampled.RData")
  # Run JUNIPER with NO unsampled intermediates
  set.seed(0)
  init <- initialize(
    indir = "brsv-outbreak/input_data/",
    n_subtrees = 1,
    init_pi = 1,
    init_mu = 3e-6,
    n_global = 1000,
    #split_bottlenecks = T,
    rooted = T,
    a_g = 7,
    ongoing = F,
    fixed_pi = T
  )

  res_complete_bot <- run_mcmc(init, T)
  save(res_complete_bot, file = "brsv-outbreak/output_complete_bot.RData")
}else{
  load("brsv-outbreak/output.RData")
  load("brsv-outbreak/output_unsampled.RData")
  load("brsv-outbreak/output_complete_bot.RData")
}

## Get adjacency matrices of transmissions
adj_juniper_brsv <- summarize(res)$direct_transmission
colnames(adj_juniper_brsv)[1] <- gsub(" \\(root\\)", "", colnames(adj_juniper_brsv)[1])
rownames(adj_juniper_brsv)[1] <- gsub(" \\(root\\)", "", rownames(adj_juniper_brsv)[1])

adj_juniper_unsampled_brsv <- summarize(res_unsampled)$direct_transmission
colnames(adj_juniper_unsampled_brsv)[1] <- gsub(" \\(root\\)", "", colnames(adj_juniper_unsampled_brsv)[1])
rownames(adj_juniper_unsampled_brsv)[1] <- gsub(" \\(root\\)", "", rownames(adj_juniper_unsampled_brsv)[1])

adj_juniper_complete_bot_brsv <- summarize(res_complete_bot)$direct_transmission
colnames(adj_juniper_complete_bot_brsv)[1] <- gsub(" \\(root\\)", "", colnames(adj_juniper_complete_bot_brsv)[1])
rownames(adj_juniper_complete_bot_brsv)[1] <- gsub(" \\(root\\)", "", rownames(adj_juniper_complete_bot_brsv)[1])

## Run outbreaker2
fasta <- read.FASTA("brsv-outbreak/input_data/aligned.fasta")
ref <- read.FASTA("brsv-outbreak/input_data/ref.fasta")
fasta <- c(ref, fasta)

dates <- as.Date(gsub(".*\\|", "", names(fasta)))
dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), units = "days"))

names(fasta) <- NULL
dat <- list(
  dna = fasta,
  dates = dates,
  w_dens = dgamma(1:30, shape = 7, rate = 1)
)

out <- outbreaker(
  data = dat,
  config = list(n_iter = 10000, sample_every = 100, move_kappa = F, pb = T)
)

# Filter to columns relating to ancestry
out <- out[, 9:17]
out <- as.matrix(out)

adj_outbreaker_brsv <- matrix(0, nrow = ncol(out), ncol = ncol(out))

# Burnin
out <- out[round(0.2 * nrow(out)):nrow(out), ]

for (r in 1:nrow(out)) {

  # Which transmissions do we update?
  inds <- cbind(out[r, ], 1:ncol(out))
  inds <- inds[!is.na(inds[, 1]), ]

  adj_outbreaker_brsv[inds] <- adj_outbreaker_brsv[inds] + 1
  #print(r)
}
adj_outbreaker_brsv <- adj_outbreaker_brsv / nrow(out)
colnames(adj_outbreaker_brsv) <- colnames(adj_juniper_brsv)
rownames(adj_outbreaker_brsv) <- rownames(adj_juniper_brsv)

## Run TransPhylo

# Read in FASTA
fasta <- read.FASTA("brsv-outbreak/input_data/aligned.fasta")

names <- gsub("\\|.*", "", names(fasta))
dates <- gsub(".*\\|", "", names(fasta))

# No mutations at consensus level in entire dataset, so let each generation of cases simply coalesce halfway back to the previous one

edges <- matrix(
  c(
    9, 1,
    9, 5,
    9, 10,
    10, 6,
    10, 7,
    10, 11,
    11, 2,
    11, 3,
    11, 4,
    11, 8
  ),
  byrow = T,
  ncol = 2
)
times <- c(7, 21, 21, 21, 7, 14, 14, 21, 3.5, 10.5, 17.5)
lengths <- c(3.5, 3.5, 7, 3.5, 3.5, 7, 3.5, 3.5, 3.5, 3.5)

phy0 <- list(edge = edges, edge.length = lengths, Nnode = 3, node.label = c("NODE_0000001", "NODE_0000002", "NODE_0000003"), tip.label = names, root.edge = 0.1)
class(phy0) <- "phylo"

# Convert to binary tree
phy0 <- multi2di(phy0)

# Make branch lengths nonzero
phy0$edge.length <- phy0$edge.length + 0.01

# Add max date
ptree<-ptreeFromPhylo(phy0,dateLastSample= as.Date("2021-01-01"))

res_TransPhylo <- inferTTree(
  ptree,
  mcmcIterations=10000,
  w.mean = 7,
  w.std = sqrt(7),
  ws.mean = 5,
  ws.std = sqrt(5),
  startOff.r = 1,
  startOff.p = 0.5,
  startPi = 0.9999, # Error out if set to 1
  updateOff.r = T,
  updateOff.p = F,
  updatePi = F,
  dateT= Inf
)

# Names according to transphylo
ref <- read.FASTA("brsv-outbreak/input_data/ref.fasta")
tp_names <- c(names(ref), res_TransPhylo[[1]]$ctree$nam)
tp_names <- gsub("\\|.*", "", tp_names)
n <- length(tp_names)

adj_transphylo_brsv <- matrix(0, ncol = n, nrow = n)
colnames(adj_transphylo_brsv) <- tp_names
rownames(adj_transphylo_brsv) <- tp_names


# Knock out burnin
res_TransPhylo <- res_TransPhylo[round(0.2 * length(res_TransPhylo)):length(res_TransPhylo)]

for (r in 1:length(res_TransPhylo)) {

  ttree=extractTTree(res_TransPhylo[[r]]$ctree)$ttree

  infectors=ttree[1:(n-1),3]
  infecteds=1:(n-1)

  # Also row for 0 -> someone
  infectors <- c(infectors, 0)
  infecteds <- c(infecteds, which(ttree[,3] == 0))

  # Reindex
  infectors <- infectors + 1
  infecteds <- infecteds + 1

  # Subset to which transmissions are between sampled individuals
  keep <- which(infectors <= n & infecteds <= n)

  infectors <- infectors[keep]
  infecteds <- infecteds[keep]

  # Reindexing and update adjacencey matrix
  adj_transphylo_brsv[cbind(infectors,infecteds)] <- adj_transphylo_brsv[cbind(infectors,infecteds)] + 1/length(res_TransPhylo)
}

# Reorder columns and rows per Juniper output
adj_transphylo_brsv <- adj_transphylo_brsv[, match(colnames(adj_juniper_brsv), colnames(adj_transphylo_brsv))]
adj_transphylo_brsv <- adj_transphylo_brsv[match(rownames(adj_juniper_brsv), rownames(adj_transphylo_brsv)), ]

## BadTrIP
fasta <- read.FASTA("brsv-outbreak/input_data/aligned.fasta")
names <- gsub("\\|.*", "", names(fasta))

net <- readLines("brsv-outbreak/badtrip/summary_network.txt")

# Get number of hosts
N <- length(fasta) + 1
start <- which(net == "Probabilities of direct transmittor to each sampled host: ")

adj_badtrip_brsv <- matrix(0, nrow = N, ncol = N)
for (r in (start+1):length(net)) {
  if(grepl("To host", net[r])){
    to <- sub("To host HH", "", net[r])
    to <- sub(" from : ", "", to)
    to <- as.numeric(to)

    from <- unlist(strsplit(net[r+1], ", "))
    probs <- sub(".* ", "", from)
    probs <- as.numeric(probs)

    from <- sub(" .*", "", from)
    from <- sub("HH", "", from)
    from[from == "Unsampled"] <- "0" # then adding 1
    from <- as.numeric(from)
    probs <- probs[from != 0]
    from <- from[from != 0]

    adj_badtrip_brsv[from, to] <- probs

  }
}

# BadTrIP was run with the ref seq at the end
adj_badtrip_brsv <- adj_badtrip_brsv[c(9, 1:8), ]
adj_badtrip_brsv <- adj_badtrip_brsv[, c(9, 1:8)]

colnames(adj_badtrip_brsv) <- colnames(adj_juniper_brsv)
rownames(adj_badtrip_brsv) <- rownames(adj_juniper_brsv)

# Manually input true links
truth <- cbind(
  c(1,8,8,8,1,6,2,8),
  2:9
)

# Rename rows and columns for simplicity
colnames(adj_juniper_brsv) <- paste("Cow", 1:9)
rownames(adj_juniper_brsv) <- paste("Cow", 1:9)
colnames(adj_juniper_unsampled_brsv) <- paste("Cow", 1:9)
rownames(adj_juniper_unsampled_brsv) <- paste("Cow", 1:9)
colnames(adj_juniper_complete_bot_brsv) <- paste("Cow", 1:9)
rownames(adj_juniper_complete_bot_brsv) <- paste("Cow", 1:9)
colnames(adj_badtrip_brsv) <- paste("Cow", 1:9)
rownames(adj_badtrip_brsv) <- paste("Cow", 1:9)
colnames(adj_outbreaker_brsv) <- paste("Cow", 1:9)
rownames(adj_outbreaker_brsv) <- paste("Cow", 1:9)
colnames(adj_transphylo_brsv) <- paste("Cow", 1:9)
rownames(adj_transphylo_brsv) <- paste("Cow", 1:9)

### Visualizations

ig <- igraph::graph_from_adjacency_matrix(adj_juniper_unsampled_brsv, mode = "directed", weighted = T)

# Extract edges in from and to
edges <- as_edgelist(ig)
status <- factor(rep("Incorrect", nrow(edges)), levels = c("Correct", "Incorrect"))

# True links

for (i in 1:nrow(truth)) {
  row <- which(edges[, 1] == paste("Cow", truth[i, 1]) & edges[, 2] == paste("Cow", truth[i, 2]))
  if(length(row) == 1){
    status[row] <- "Correct"
  }
}

E(ig)$Inference <- status
E(ig)$Probability <- E(ig)$weight


set.seed(1)
brsv_network <- ggraph(ig, layout = "fr") +
  geom_edge_link(
    aes(width = Probability, alpha = Probability, color = Inference),
    arrow = arrow(type = "open", length = unit(4, "mm")),
    #start_cap = circle(5, 'mm'),
    end_cap = circle(8, 'mm')
  ) +
  geom_node_point(size = 14.3, shape = 21, fill = 'white', stroke = 1.2) +
  geom_node_text(aes(label = name), size = 3.5) +
  scale_edge_width_continuous(range = c(0.5, 4)) +
  scale_edge_alpha_continuous(range = c(0.25, 1)) +
  scale_edge_color_manual(values = c("#11CC55", "darkgray")) +
  theme_graph(background = "white", base_family = "sans") +
  theme(legend.position = c(0.85, 0.75))

## Epi link posterior probabilities by method
models <- c("JUNIPER, standard", "JUNIPER, incomplete bottlenecks", "JUNIPER, no intermediates", "BadTrIP", "outbreaker2", "TransPhylo")
brsv_link_df <- data.frame(
  Transmission = paste("Cow", truth[, 1], "to", truth[, 2]),
  Probability = c(
    adj_juniper_brsv[truth],
    adj_juniper_unsampled_brsv[truth],
    adj_juniper_complete_bot_brsv[truth],
    adj_badtrip_brsv[truth],
    adj_outbreaker_brsv[truth],
    adj_transphylo_brsv[truth]
  ),
  Model = factor(rep(models, each = nrow(truth)), levels = models[c(2,1,3:6)])
)

# Mean of posterior probabilities assigned to epi links
means_brsv <- c(
  mean(adj_juniper_brsv[truth]),
  mean(adj_juniper_unsampled_brsv[truth]),
  mean(adj_juniper_complete_bot_brsv[truth]),
  mean(adj_outbreaker_brsv[truth]),
  mean(adj_badtrip_brsv[truth]),
  mean(adj_transphylo_brsv[truth])
)

brsv_links <- ggplot(data = brsv_link_df) +
  geom_tile(mapping = aes(x = Model, y = Transmission, alpha = Probability), fill = "#11CC55", color = "black", linewidth = 0.4) +
  geom_text(
    data = data.frame(
      Model = factor(models, levels = models),
      Value = paste0("µ = ", round(means_brsv, 3))
    ),
    aes(
      x = Model,
      y = 8.7,
      label = Value
    ),
    size = 3
  ) +
  scale_alpha_continuous(range = c(0,1)) +
  scale_x_discrete(labels=function(x){gsub("\\s", "\n", x)}) +
  ylim(c(unique(brsv_link_df$Transmission), "")) +
  theme_minimal() +
  theme(plot.margin = margin(-1,0,0,0, "cm"), panel.grid = element_blank())

### South Africa cluster
metadata <- read_xlsx("south-africa-outbreak/raw/leger.xlsx", skip = 1)

if(regenerate){
  # Run JUNIPER
  set.seed(0)
  init <- initialize(
    filters = list(af = 0.05, sb = 10, dp = 100), # Because iSNVs were pre-filtered to 5% AF
    indir = "south-africa-outbreak/input_data/",
    n_subtrees = 1,
    init_mu = 1e-5,
    init_N_eff = 0.1,
    init_pi = 1,
    init_R = 0.5,
    n_global = 1000,
    split_bottlenecks = F,
    rooted = T,
    a_g = 10,
    ongoing = F,
    fixed_pi = T
  )

  res <- run_mcmc(init, T)
  save(res, file = "south-africa-outbreak/output.RData")
  # Run JUNIPER without assuming perfect sampling, i.e. normal way
  set.seed(0)
  init <- initialize(
    filters = list(af = 0.05, sb = 10, dp = 100), # Because iSNVs were pre-filtered to 5% AF
    indir = "south-africa-outbreak/input_data/",
    n_subtrees = 1,
    init_mu = 1e-5,
    init_N_eff = 0.1,
    #init_pi = 1,
    init_R = 0.5,
    n_global = 1000,
    split_bottlenecks = F,
    rooted = T,
    a_g = 10,
    ongoing = F #,
    #fixed_pi = T
  )

  res_unsampled <- run_mcmc(init, T)
  save(res_unsampled, file = "south-africa-outbreak/output_unsampled.RData")
  # Run JUNIPER without assuming perfect sampling, i.e. normal way, plus incomplete bottlenecks
  set.seed(0)
  init <- initialize(
    filters = list(af = 0.05, sb = 10, dp = 100), # Because iSNVs were pre-filtered to 5% AF
    indir = "south-africa-outbreak/input_data/",
    n_subtrees = 1,
    init_mu = 1e-5,
    init_N_eff = 0.1,
    #init_pi = 1,
    init_R = 0.5,
    n_global = 1000,
    split_bottlenecks = T,
    rooted = T,
    a_g = 10,
    ongoing = F #,
    #fixed_pi = T
  )

  res_split <- run_mcmc(init, T)
  save(res_split, file = "south-africa-outbreak/output_split.RData")
}else{
  load("south-africa-outbreak/output.RData")
  load("south-africa-outbreak/output_unsampled.RData")
  load("south-africa-outbreak/output_split.RData")
}


## Run outbreaker2
fasta <- read.FASTA("south-africa-outbreak/input_data/aligned.fasta")
ref <- read.FASTA("south-africa-outbreak/input_data/ref.fasta")
fasta <- c(ref, fasta)

dates <- as.Date(gsub(".*\\|", "", names(fasta)))
dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), units = "days"))

names(fasta) <- NULL
dat <- list(
  dna = fasta,
  dates = dates,
  w_dens = dgamma(1:30, shape = 10, rate = 1)
)

out <- outbreaker(
  data = dat,
  config = list(n_iter = 10000, sample_every = 100, move_kappa = F, pb = T)
)

#save(out, file = "south-africa-outbreak/outbreaker.RData")
#load("south-africa-outbreak/outbreaker.RData")



if(regenerate){
  ## Run IQTree + TransPhylo

  # Build a combined fasta/ref
  fasta <- read.FASTA("south-africa-outbreak/input_data/aligned.fasta")
  ref <- read.FASTA("south-africa-outbreak/input_data/ref.fasta")
  # combined <- c(ref,fasta) -- TreeTime unable to handle this; returns error
  combined <- fasta
  dates <- as.Date(gsub(".*\\|", "", names(combined)))
  names <- gsub("\\|.*", "", names(combined))
  names(combined) <- names
  write.FASTA(combined, "south-africa-outbreak/input_data/combined.fasta")

  # Run IQTree
  system(
    "iqtree2 -s south-africa-outbreak/input_data/combined.fasta -m JC --prefix south-africa-outbreak/iqtree/iqtree-results"
  )



  # Convert dates to years, transphylo doesn't like small date spans in years
  dates <- as.numeric(difftime(dates, as.Date("2000-01-01"), units = 'days'))
  max_year <- 2000 + max(dates)
  dates <- as.Date("2000-01-01") + years(dates)

  date_csv <- data.frame(name = names, date = dates)
  write.csv(date_csv, file = paste0("south-africa-outbreak/iqtree/date.csv"), row.names = F, quote = F)

  # Run TreeTime
  system(
    "treetime --tree south-africa-outbreak/iqtree/iqtree-results.treefile --dates south-africa-outbreak/iqtree/date.csv --aln south-africa-outbreak/input_data/combined.fasta --outdir south-africa-outbreak/treetime --coalescent skyline --n-skyline 2 --greedy-resolve"
  )



  # Load timed tree
  phy <- read.nexus("south-africa-outbreak/treetime/timetree.nexus")

  # Convert to binary tree
  phy <- multi2di(phy)

  # Make branch lengths nonzero
  phy$edge.length <- phy$edge.length + 0.01

  # Add max date
  ptree<-ptreeFromPhylo(phy,dateLastSample= max_year)
  res_TransPhylo <- inferTTree(
    ptree,
    mcmcIterations=10000,
    w.mean = 10,
    w.std = sqrt(10),
    ws.mean = 5,
    ws.std = sqrt(5),
    startOff.r = 0.5,
    startOff.p = 0.5,
    startPi = 0.9999, # Error out if set to 1
    updateOff.r = T,
    updateOff.p = F,
    updatePi = F,
    dateT= Inf
  )
  save(res_TransPhylo, file = "south-africa-outbreak/res_TransPhylo.RData")
}else{
  load("south-africa-outbreak/res_TransPhylo.RData")
}






## Obtain adjacency matrices

# For Juniper
s <- summarize(res)
adj_juniper <- s$direct_transmissions
colnames(adj_juniper)[1] <- gsub(" \\(root\\)", "", colnames(adj_juniper)[1])
rownames(adj_juniper)[1] <- gsub(" \\(root\\)", "", rownames(adj_juniper)[1])
colnames(adj_juniper) <- metadata$`Outbreak ID`[match(colnames(adj_juniper), metadata$`Gisaid Accession`)]
rownames(adj_juniper) <- metadata$`Outbreak ID`[match(rownames(adj_juniper), metadata$`Gisaid Accession`)]

ss <- summarize(res_unsampled)
adj_juniper_unsampled <- ss$indirect_transmissions
colnames(adj_juniper_unsampled)[1] <- gsub(" \\(root\\)", "", colnames(adj_juniper_unsampled)[1])
rownames(adj_juniper_unsampled)[1] <- gsub(" \\(root\\)", "", rownames(adj_juniper_unsampled)[1])
colnames(adj_juniper_unsampled) <- metadata$`Outbreak ID`[match(colnames(adj_juniper_unsampled), metadata$`Gisaid Accession`)]
rownames(adj_juniper_unsampled) <- metadata$`Outbreak ID`[match(rownames(adj_juniper_unsampled), metadata$`Gisaid Accession`)]

sss <- summarize(res_split)
adj_juniper_split <- sss$indirect_transmissions
colnames(adj_juniper_split)[1] <- gsub(" \\(root\\)", "", colnames(adj_juniper_split)[1])
rownames(adj_juniper_split)[1] <- gsub(" \\(root\\)", "", rownames(adj_juniper_split)[1])
colnames(adj_juniper_split) <- metadata$`Outbreak ID`[match(colnames(adj_juniper_split), metadata$`Gisaid Accession`)]
rownames(adj_juniper_split) <- metadata$`Outbreak ID`[match(rownames(adj_juniper_split), metadata$`Gisaid Accession`)]


# For BadTrIP
fasta <- read.FASTA("south-africa-outbreak/badtrip/old_input_data/aligned.fasta")
names <- gsub("\\|.*", "", names(fasta))

net <- readLines("south-africa-outbreak/badtrip/summary_network.txt")

# Get number of hosts
N <- length(fasta)
start <- which(net == "Probabilities of direct transmittor to each sampled host: ")

adj_badtrip <- matrix(0, nrow = N, ncol = N)
for (r in (start+1):length(net)) {
  if(grepl("To host", net[r])){
    to <- sub("To host HH", "", net[r])
    to <- sub(" from : ", "", to)
    to <- as.numeric(to)

    from <- unlist(strsplit(net[r+1], ", "))
    probs <- sub(".* ", "", from)
    probs <- as.numeric(probs)

    from <- sub(" .*", "", from)
    from <- sub("HH", "", from)
    from[from == "Unsampled"] <- "0" # then adding 1
    from <- as.numeric(from)
    probs <- probs[from != 0]
    from <- from[from != 0]

    adj_badtrip[from, to] <- probs

  }
}


colnames(adj_badtrip) <- metadata$`Outbreak ID`[match(names, metadata$Name)]
rownames(adj_badtrip) <- metadata$`Outbreak ID`[match(names, metadata$Name)]

colnames(adj_badtrip)[22]
colnames(adj_badtrip)[22] <- "P11B"
rownames(adj_badtrip)[22] <- "P11B"

# Reorder columns and rows per Juniper output
adj_badtrip <- adj_badtrip[, match(colnames(adj_juniper), colnames(adj_badtrip))]
adj_badtrip <- adj_badtrip[match(rownames(adj_juniper), rownames(adj_badtrip)), ]

## For outbreaker2

# Filter to columns relating to ancestry
out <- out[, 9:32]
out <- as.matrix(out)

adj_outbreaker <- matrix(0, nrow = ncol(out), ncol = ncol(out))

# Burnin
out <- out[round(0.2 * nrow(out)):nrow(out), ]

for (r in 1:nrow(out)) {

  # Which transmissions do we update?
  inds <- cbind(out[r, ], 1:ncol(out))
  inds <- inds[!is.na(inds[, 1]), ]

  adj_outbreaker[inds] <- adj_outbreaker[inds] + 1
  #print(r)
}

adj_outbreaker <- adj_outbreaker / nrow(out)

fasta <- read.FASTA("south-africa-outbreak/input_data/aligned.fasta")
ref <- read.FASTA("south-africa-outbreak/input_data/ref.fasta")
fasta <- c(ref, fasta)
names <- gsub("\\|.*", "", names(fasta))
colnames(adj_outbreaker) <- names
rownames(adj_outbreaker) <- names
colnames(adj_outbreaker) <- metadata$`Outbreak ID`[match(colnames(adj_outbreaker), metadata$`Gisaid Accession`)]
rownames(adj_outbreaker) <- metadata$`Outbreak ID`[match(rownames(adj_outbreaker), metadata$`Gisaid Accession`)]


## For IQTree + TransPhylo

# Names according to transphylo
ref <- read.FASTA("south-africa-outbreak/input_data/ref.fasta")
tp_names <- c(names(ref), res_TransPhylo[[1]]$ctree$nam)
tp_names <- gsub("\\|.*", "", tp_names)
n <- length(tp_names)

adj_transphylo <- matrix(0, ncol = n, nrow = n)
colnames(adj_transphylo) <- tp_names
rownames(adj_transphylo) <- tp_names


# Knock out burnin
res_TransPhylo <- res_TransPhylo[round(0.2 * length(res_TransPhylo)):length(res_TransPhylo)]

for (r in 1:length(res_TransPhylo)) {

  ttree=extractTTree(res_TransPhylo[[r]]$ctree)$ttree

  infectors=ttree[1:(n-1),3]
  infecteds=1:(n-1)

  # Also row for 0 -> someone
  infectors <- c(infectors, 0)
  infecteds <- c(infecteds, which(ttree[,3] == 0))

  # Reindex
  infectors <- infectors + 1
  infecteds <- infecteds + 1

  # Subset to which transmissions are between sampled individuals
  keep <- which(infectors <= n & infecteds <= n)

  infectors <- infectors[keep]
  infecteds <- infecteds[keep]

  # Reindexing and update adjacencey matrix
  adj_transphylo[cbind(infectors,infecteds)] <- adj_transphylo[cbind(infectors,infecteds)] + 1/length(res_TransPhylo)
}

# Match column and row names
colnames(adj_transphylo) <- metadata$`Outbreak ID`[match(colnames(adj_transphylo), metadata$`Gisaid Accession`)]
rownames(adj_transphylo) <- metadata$`Outbreak ID`[match(rownames(adj_transphylo), metadata$`Gisaid Accession`)]

# Reorder columns and rows per Juniper output
adj_transphylo <- adj_transphylo[, match(colnames(adj_juniper), colnames(adj_transphylo))]
adj_transphylo <- adj_transphylo[match(rownames(adj_juniper), rownames(adj_transphylo)), ]



### Visualizations

# List facilities by patient
facilities <- metadata$Facility[match(colnames(adj_juniper), metadata$`Outbreak ID`)]
# Correct P11 vs P11B
facilities[10] <- "NW"
# Split multiple facilities
facilities <- sapply(facilities, function(s){strsplit(s, ",")}, USE.NAMES = F)

ig <- igraph::graph_from_adjacency_matrix(adj_juniper, mode = "directed", weighted = T)

# Extract edges in from and to
edges <- as_edgelist(ig)
status <- factor(rep("Neither", nrow(edges)), levels = c("Epidemiologically Confirmed", "Shared Facility", "Neither"))

links <- read_xlsx("south-africa-outbreak/raw/links.xlsx")

for (i in 1:nrow(edges)) {
  fac1 <- facilities[[match(edges[i, 1], colnames(adj_juniper))]]
  fac2 <- facilities[[match(edges[i, 2], colnames(adj_juniper))]]
  if(length(intersect(fac1, fac2)) > 0){
    status[i] = "Shared Facility"
  }
}

# Epi links

for (i in 1:nrow(links)) {
  row <- which(edges[, 1] == links$from[i] & edges[, 2] == links$to[i])
  if(length(row) == 1){
    status[row] <- "Epidemiologically Confirmed"
  }
}

E(ig)$Status <- status
E(ig)$Probability <- E(ig)$weight


short_names <- gsub("B", "", V(ig)$name)

set.seed(36)
sa_network <- ggraph(ig, layout = "fr") +
  geom_edge_link(
    aes(width = Probability, alpha = Probability, color = Status),
    arrow = arrow(type = "open", length = unit(3, "mm")),
    #start_cap = circle(5, 'mm'),
    end_cap = circle(6, 'mm')
  ) +
  geom_node_point(size = 10, shape = 21, fill = 'white', stroke = 1.2) +
  geom_node_text(aes(label = short_names), size = 3) +
  scale_edge_width_continuous(range = c(0.5, 2.5)) +
  scale_edge_alpha_continuous(range = c(0.25, 1)) +
  scale_edge_color_manual(values = c("#11CC55", "#4477EE", "darkgray")) +
  theme_graph(background = "white", base_family = "sans") +
  theme(legend.position = c(0.85, 0.25))

#print(sa_network)





## Epi link posterior probabilities by method
models <- c("JUNIPER, no intermediates", "JUNIPER, standard", "JUNIPER, incomplete bottlenecks", "outbreaker2", "BadTrIP", "TransPhylo")

sa_link_df <- data.frame(
  Transmission = paste(links$from, "to", links$to),
  Probability = c(
    adj_juniper[as.matrix(links)],
    adj_juniper_unsampled[as.matrix(links)],
    adj_juniper_split[as.matrix(links)],
    adj_outbreaker[as.matrix(links)],
    adj_badtrip[as.matrix(links)],
    adj_transphylo[as.matrix(links)]
  ),
  Model = rep(factor(models, levels = models[c(1,4,5,3,6,2)]), each = nrow(links))
)


# Mean of posterior probabilities assigned to epi links
means_sa <- c(
  mean(adj_juniper[as.matrix(links)]),
  mean(adj_juniper_unsampled[as.matrix(links)]),
  mean(adj_juniper_split[as.matrix(links)]),
  mean(adj_outbreaker[as.matrix(links)]),
  mean(adj_badtrip[as.matrix(links)]),
  mean(adj_transphylo[as.matrix(links)])
)

sa_links <- ggplot(data = sa_link_df) +
  geom_tile(mapping = aes(x = Model, y = Transmission, alpha = Probability), fill = "#11CC55", color = "black", linewidth = 0.4) +
  geom_text(
    data = data.frame(
      Model = factor(models, levels = models),
      Value = paste0("µ = ", signif(means_sa, 3))
    ),
    aes(
      x = Model,
      y = 9.7,
      label = Value
    ),
    size = 3
  ) +
  scale_alpha_continuous(range = c(0,1)) +
  scale_x_discrete(labels=function(x){gsub("\\s", "\n", x)}) +
  ylim(c(unique(sa_link_df$Transmission), "")) +
  theme_minimal() +
  theme(plot.margin = margin(-1,0,0,0, "cm"), panel.grid = element_blank())


cowplot::plot_grid(

  cowplot::plot_grid(
    brsv_links,
    sa_links,
    align = "hv",
    labels = c("A", "C"),
    ncol = 1
  ),
  cowplot::plot_grid(
    brsv_network,
    sa_network,
    labels = c("B", "D"),
    ncol = 1
  ),
  ncol = 2,
  #rel_widths = c(2,2),
  labels = ""
)

ggsave("./figs/brsv-and-south-africa.pdf", width = 12, height = 12, limitsize = F)
ggsave("./figs/brsv-and-south-africa.png", width = 12, height = 12, limitsize = F)


