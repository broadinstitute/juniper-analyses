### Visualize MCMC output

library(ggplot2)
library(cowplot)
library(ape)
library(juniper0)

Rcpp::sourceCpp("cpp_subroutines.cpp")


## Lineage by case
# Read in sequence names and metadata
cons <- read.FASTA("mass-10k/input_data/aligned.fasta")
n_obs <- length(cons) + 1
names <- names(cons)
names <- gsub("\\|.*", "", names)
meta <- read.csv("mass-10k/metadata.txt", sep = "\t")

# Vaccination status
vaxed <- rep(NA, n_obs)
vaxed[2:n_obs] <- meta$Vaccination.Status[match(names, meta$Sequence.Name)]

# Load in preprocessed file of number of offspring at each MCMC iteration, which was run separately on a VM
load("mass-10k/n_kids.RData")
n_iters <- length(n_kids)

n_kids <- n_kids[25001:50000]

# Indices of each category
boosted <- which(vaxed == "boosted")
fully <- which(vaxed == "fully vaccinated")
partial <- which(vaxed == "partially vaccinated")
unvax <- which(vaxed == "unvaccinated")
vax <- which(vaxed %in% c("boosted", "fully vaccinated", "partially vaccinated"))



### Correct individual-level R calculation
load("mass-10k/ts.RData")

# Get R for each iter
logging <- readLines("mass-10k/log.out")
Rs <- logging[which(grepl("Reproductive number:", logging))]
Rs <- gsub(".*\\: ", "", Rs)
Rs <- gsub("\\ .*", "", Rs)
Rs <- gsub("\\\"", "", Rs)
Rs <- as.numeric(Rs)
Rs <- Rs[25001:50000]

pis <- logging[which(grepl("Sampling rate:", logging))]
pis <- gsub(".*\\: ", "", pis)
pis <- gsub("\\ .*", "", pis)
pis <- gsub("\\\"", "", pis)
pis <- as.numeric(pis)
pis <- pis[25001:50000]

for (i in 1:(n_iters/2)) {
  n_kids[[i]] <- n_kids[[i]][1:n_obs]
}

# Mean number of transmissions per category
mean_trans_boosted <- rep(0, n_iters/2)
mean_trans_fully <- rep(0, n_iters/2)
mean_trans_partial <- rep(0, n_iters/2)
mean_trans_unvax <- rep(0, n_iters/2)
mean_trans_vax <- rep(0, n_iters/2)

# PDF of transmissions by category
p_trans_boosted <- rep(list(rep(0, 7)), n_iters/2)
p_trans_fully <- rep(list(rep(0, 7)), n_iters/2)
p_trans_partial <- rep(list(rep(0, 7)), n_iters/2)
p_trans_unvax <- rep(list(rep(0, 7)), n_iters/2)
p_trans_vax <- rep(list(rep(0, 7)), n_iters/2)

if(FALSE){
  for (i in 1:(n_iters/2)) {
    mean_trans <- rep(0, n_obs)
    min_t <- min(ts[[i]])
    wbar0 <- rev(wbar(min_t - 1, 0, Rs[i] * 0.5 / (1 - 0.5), 1 - 0.5, pis[i], 5, 1, 5, 1, 0.1))

    rho <- Rs[i]
    psi <- 0.5

    # Get correct wbar for each j
    ws <- wbar0[round(-ts[[i]]/0.1)] # Still log scale

    # Normalizing constant
    norms <- rep(0, n_obs)
    for (j in 1:n_obs) {
      norms[j] <- exp(alpha(n_kids[[i]][j], psi, rho, ws[j])) # Takes in wbar on log scale
    }

    # P(n_kids = k) for each host, for k = 0, 1, ..., 6
    for (k in 1:7) {
      ps <- choose(k-1, n_kids[[i]]) * dnbinom(k-1, rho, psi) * exp(ws)^(k - 1 - n_kids[[i]]) / norms

      p_trans_boosted[[i]][k] <- mean(ps[boosted])
      #p_trans_fully[[i]][k] <- mean(ps[fully])
      #p_trans_partial[[i]][k] <- mean(ps[partial])
      p_trans_vax[[i]][k] <- mean(ps[vax])
      p_trans_unvax[[i]][k] <- mean(ps[unvax])
    }


    # Sum_{k=d}^\infty (k choose d) alpha(k) wbar^(k-d) * k / normalizing constant, majorly simplifies to
    means <- rep(0, n_obs)
    for (j in 1:n_obs) {
      means[j] <- exp(alpha2(n_kids[[i]][j], psi, rho, ws[j])) / norms[j] # Takes in wbar on log scale
    }

    mean_trans_boosted[i] <- mean(means[boosted])
    #mean_trans_fully[i] <- mean(means[fully])
    #mean_trans_partial[i] <- mean(means[partial])
    mean_trans_vax[i] <- mean(means[vax])
    mean_trans_unvax[i] <- mean(means[unvax])

    if(i %% 100 == 0){
      print(i)
    }
  }

  save(mean_trans_boosted, file = "mass-10k/mean_trans_boosted.RData")
  save(mean_trans_vax, file = "mass-10k/mean_trans_vax.RData")
  save(mean_trans_unvax, file = "mass-10k/mean_trans_unvax.RData")
}

load("mass-10k/mean_trans_boosted.RData")
load("mass-10k/mean_trans_vax.RData")
load("mass-10k/mean_trans_unvax.RData")





# Vaccination
vax_df <- data.frame(
  Transmissions = c(mean_trans_boosted, mean_trans_vax, mean_trans_unvax),
  Status = factor(c(rep("Boosted", 0.5*n_iters), rep("Vaccinated", 0.5*n_iters), rep("Unvaccinated", 0.5*n_iters)), levels = c("Unvaccinated", "Vaccinated", "Boosted"))
)

trans_plot <- ggplot2::ggplot(data = vax_df, mapping = aes(x = Transmissions, fill = Status, color = Status)) +
  geom_density(color = NA, alpha = 0.9) +
  xlab("Transmissions per Case") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position= c(0.15, 0.75))

mean(vax_df$Transmissions[vax_df$Status == "Boosted"])
mean(vax_df$Transmissions[vax_df$Status == "Vaccinated"])
mean(vax_df$Transmissions[vax_df$Status == "Unvaccinated"])

# Relative change in transmission rate
(mean(vax_df$Transmissions[vax_df$Status == "Boosted"]) - mean(vax_df$Transmissions[vax_df$Status == "Unvaccinated"])) / mean(vax_df$Transmissions[vax_df$Status == "Unvaccinated"])
(mean(vax_df$Transmissions[vax_df$Status == "Vaccinated"]) - mean(vax_df$Transmissions[vax_df$Status == "Unvaccinated"])) / mean(vax_df$Transmissions[vax_df$Status == "Unvaccinated"])



### CDF Plot

cdf_trans_boosted <- lapply(p_trans_boosted, cumsum)
cdf_trans_vax <- lapply(p_trans_vax, cumsum)
cdf_trans_unvax <- lapply(p_trans_unvax, cumsum)



cdf_df <- data.frame(
  Transmissions = c(unlist(cdf_trans_boosted), unlist(cdf_trans_vax), unlist(cdf_trans_unvax)),
  Status = factor(rep(c("Boosted", "Vaccinated", "Unvaccinated"), each = 7*n_iters/2), levels = c("Unvaccinated", "Vaccinated", "Boosted")),
  k = as.factor(0:6)
)

cdf_plot <- ggplot(cdf_df, aes(group = interaction(Status, k), x = k, y = Transmissions, fill = Status)) +
  geom_violin(scale = "width", color = NA, adjust = 5, position = position_dodge(width = 0), alpha = 0.5) +
  stat_summary(aes(color = Status), fun = "mean", geom = "point", alpha = 1) +
  xlab("Transmissions per Case, k") +
  ylab("Pr(X <= k)") +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.35))


### Giant radial plot

## Helper function
# DFS order, expanding nodes by lineage
dfs_lineage <- function(kids, lineage){
  stack <- 1
  explored <- rep(0, length(unlist(kids)) + 1)
  i <- 1 # counter
  while (length(stack) > 0) {
    who <- kids[[stack[1]]]
    who <- who[sort.int(lineage[who], index.return = T)$ix]
    explored[i] <- stack[1]
    stack <- stack[-1]
    stack <- c(who, stack)
    i <- i+1
  }
  return(explored)
}

# BFS traversal of tree
bfs <- function(h){
  n <- length(h)
  kids <- rep(list(integer(0)), n)

  for (i in 2:n) {
    kids[[h[i]]] <- c(kids[[h[i]]], i)
  }

  out <- rep(1, n)
  frontier <- kids[[1]]
  counter <- 2
  while (length(frontier) > 0) {
    out[counter:(counter + length(frontier) - 1)] <- frontier
    counter <- counter + length(frontier)
    frontier <- unlist(kids[frontier])
  }
  return(out)
}

# Load in sample for visualization
load("mass-10k/snapshot.RData")
n <- snapshot$n
h <- snapshot$h

# Get kids of each node, to speed up DFS
kids <- rep(list(integer(0)), n)
for (i in 2:n) {
  kids[[h[i]]] <- c(kids[[h[i]]], i)
}

# Number of kids
n_kids <- sapply(kids, length)

# Lineage
lin <- meta$Sequence.Pango.Lineage[match(names, meta$Sequence.Name)]
lineage <- rep("Unsampled", n)
lineage[2:(length(lin) + 1)] <- lin



ord <- rev(dfs_lineage(kids, lineage))

# Leaves
leaves <- which(!(1:n %in% h))
n_leaves <- length(leaves)

# Angle of each node
thetas <- rep(0, length(ord))
leaf_count <- 0

for (i in ord[ord %in% leaves]) {
  thetas[i] <- leaf_count * 2 * pi / n_leaves
  leaf_count <- leaf_count + 1
}

for (i in ord[!(ord %in% leaves)]) {
    kids_i <- kids[[i]]
    thetas[i] <- mean(thetas[kids_i])


  #print(i)
}

# Reset unsampled lineages to unsampled
lineage[(length(lin) + 2):length(lineage)] <- "Unsampled"
lineage[1] <- "Unsampled"
#lineage[lineage == "recombinant"] <- "Recombinant"

t <- sapply(snapshot$seq, function(v){v[1]})
t <- t - min(t)
t <- log(1 + t)


offset <- 4.7
rotation <- 2.9 #radians
df <- data.frame(x = ((t-offset)*cos(thetas + rotation))[t>offset], y = ((t-offset)*sin(thetas + rotation))[t>offset], lineage = lineage[t>offset])
df_standard <- data.frame(x = t, y = thetas, lineage = lineage)

df <- df[nrow(df):1, ]

df_sampled <- df[df$lineage != "Unsampled", ]

# Set lineages with <100 cases to "other"
minor <- names(table(df$lineage))[table(df$lineage) < 100]
df_sampled$lineage[df_sampled$lineage %in% minor] <- "Other"
df_standard$lineage[df_standard$lineage %in% minor] <- "Other"
df_standard <- df_standard[df_standard$lineage != "Unsampled", ]



colors <- colorRampPalette(c('red', 'green', 'blue', 'red'))(length(unique(df_sampled$lineage)))[1:(length(unique(df_sampled$lineage)) - 1)]
colors <- c(colors, "gray")

colnames(df_sampled[3]) <- "Lineage"

radial <- ggplot(df_sampled, aes(x=x,y=y,color=lineage)) +
  geom_point(size = 1, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = colors, name = "Lineage") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.direction = "vertical",
    legend.position=c(0.5, -0.15),
    legend.text = element_text(size=8),
    legend.key.width= unit(0.1, 'cm'),
    legend.key.height= unit(0.5, 'cm'),
    plot.margin = margin(1,0,3,0, "cm"),
    legend.title = element_text(hjust = 0.5, vjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 8))


#ggsave("./figs/100k.pdf", width=10, height = 6)
#ggsave("./figs/100k.png", width=10, height = 6)

### H5N1

# Juniper output
load("h5n1/res.RData")

ss <- juniper0::summarize(res, burnin = 0.5)

N_iters <- length(res[[1]])

# Extract a single posterior tree, for visualization
snapshot <- res[[2]][[N_iters]]

n <- snapshot$n
h <- snapshot$h

# Get kids of each node, to speed up DFS
kids <- rep(list(integer(0)), n)
for (i in 2:n) {
  kids[[h[i]]] <- c(kids[[h[i]]], i)
}

# Number of kids
n_kids <- sapply(kids, length)

load("h5n1/raw/meta_subsampled.RData")
cons <- read.FASTA("h5n1/input_data/aligned.fasta")
names <- gsub("\\|.*", "", names(cons))

# State of each case
case_state <- sub("^[^/]*/[^/]*/", "", meta$strain)
case_state <- gsub("\\/.*", "", case_state)
abbreviated <- which(case_state %in% state.abb)
unspaced_state <- gsub(" ", "", state.name)
case_state[abbreviated] <- unspaced_state[match(case_state[abbreviated], state.abb)]
case_state <- state.name[match(case_state, unspaced_state)]

state <- rep("ZZZ", n)
state[2:(length(cons) + 1)] <- case_state
state[is.na(state)] <- "Unknown"

ord <- rev(dfs_lineage(kids, state))

# Leaves
leaves <- which(!(1:n %in% h))
n_leaves <- length(leaves)

# Angle of each node
thetas <- rep(0, length(ord))
leaf_count <- 0

for (i in ord[ord %in% leaves]) {
  thetas[i] <- leaf_count * 2 * pi / n_leaves
  leaf_count <- leaf_count + 1
}

for (i in ord[!(ord %in% leaves)]) {
  kids_i <- kids[[i]]
  thetas[i] <- mean(thetas[kids_i])
}

t <- sapply(snapshot$seq, function(v){v[1]})

t0 <- gsub(".*\\|", "", names(cons))
t0 <- max(as.Date(t0))
t <- t + t0


df_standard <- data.frame(x = t, y = thetas, state = state)
df_standard <- df_standard[df_standard$state != "ZZZ", ]


small <- names(which(table(df_standard$state) < 50))
df_standard$state[df_standard$state %in% small] <- "Other"
df_standard$state <- factor(df_standard$state, levels = c(sort(setdiff(unique(df_standard$state), "Other")), "Other"))


## Phylo tree plot

# vertical segments
xs <- as.Date(integer(0))
ystart <- c()
yend <- c()
for (i in 1:n) {
  kids <- which(h == i)
  if(length(kids) > 0){
    xs <- c(xs, t[i])
    ystart <- c(ystart, min(thetas[kids]))
    yend <- c(yend, max(thetas[kids]))
  }
}



h5n1_phylo <- ggplot() +
  geom_segment(mapping = aes(x = (t[h])[-1], xend = t[-1], y = thetas[-1], yend = thetas[-1]), linewidth = 0.2) +
  geom_segment(mapping = aes(x = xs, xend = xs, y = ystart, yend = yend), linewidth = 0.2) +
  geom_point(mapping= aes(x = df_standard$x, y = df_standard$y, color = df_standard$state), size = 0.5) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(length(unique(df_standard$state)) - 1, "Dark2"), "grey"), name = "State") +
  scale_shape(name = "Host") +
  xlab("Date") +
  scale_x_date(date_labels = "%b %Y") +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    legend.position = c(0.85, 0.85)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))

h5n1_phylo

## Subset to California cases and plot transmission network
ca_direct <- ss$direct_transmissions[which(state == "California") - 1, which(state == "California") - 1]
ca_direct[ca_direct < 0.5] <- 0
keep <- which(colSums(ca_direct) > 0 | rowSums(ca_direct) > 0)
ca_direct <- ca_direct[keep, keep]
ig <- igraph::graph_from_adjacency_matrix(ca_direct, mode = "directed", weighted = T)

E(ig)$Probability <- E(ig)$weight


state_count <- table(state[state %in% state.name])


big_states <- table(case_state)
big_states <- names(big_states[which(big_states >= 50)])


set.seed(0)
ca_network <- ggraph(ig, layout = "fr") +
  geom_edge_link(
    aes(width = Probability, alpha = Probability),
    arrow = arrow(type = "open", length = unit(1, "mm")),
    #start_cap = circle(5, 'mm'),
    end_cap = circle(1.5, 'mm')
  ) +
  geom_node_point(size = 2.5, color = RColorBrewer::brewer.pal(sum(state_count >= 50), "Dark2")[1]) +
  scale_edge_width_continuous(range = c(0.25, 0.75)) +
  scale_edge_alpha_continuous(range = c(0.25, 1)) +
  theme_graph(background = "white", base_family = "sans") +
  theme(legend.position = "top")

# Posterior density of various params



### Transmission analysis of h5n1
n_iters <- length(res[[1]])


## Calculate mean offspring by state
n_obs <- 1519

load("h5n1/mean_trans_state.RData")

h5n1_trans_df <- data.frame(
  Offspring = unlist(mean_trans_state),
  State = factor(big_states, levels = rev(big_states))
)


h5n1_trans <- ggplot(h5n1_trans_df, aes(x = Offspring, y = State, fill = State)) +
  geom_violin(adjust = 2, color = NA) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(length(big_states), "Dark2"))) +
  theme_minimal() +
  xlab("Transmissions per Case") +
  theme(legend.position = "none")

cowplot::plot_grid(
  cowplot::plot_grid(
    h5n1_phylo,
    cowplot::plot_grid(
      h5n1_trans,
      #mu_dens,
      #R_dens,
      #pi_dens,
      #N_eff_dens + xlab("Intrahost doubling time (days)"),
      ca_network,
      ncol = 1,
      labels = c("B", "C"),
      rel_heights = c(3,4)
    ),
    labels = c("A",""),
    rel_widths = c(3,2)
  ),
  cowplot::plot_grid(
    radial,
    cowplot::plot_grid(
      trans_plot,
      cdf_plot,
      ncol = 1,
      labels = c("E", "F")
    ),
    labels = c("D","")
  ),
  ncol = 1,
  rel_heights = c(3,2)
)

ggsave("./figs/h5n1-and-mass-10k.pdf", width = 9, height = 12, limitsize = F)
ggsave("./figs/h5n1-and-mass-10k.png", width = 9, height = 12, limitsize = F)


### Performance / convergence of mass 10k
logging <- readLines("mass-10k/log.out")

param_names <- c(
  "Evolution rate",
  "Within-host effective population size",
  "Sampling rate",
  "Reproductive number"
)

params <- list()

for (i in 1:length(param_names)) {
  params[[i]] <- logging[which(grepl(param_names[i], logging))]
  params[[i]] <- gsub(".*\\: ", "", params[[i]])
  params[[i]] <- gsub("\\ .*", "", params[[i]])
  params[[i]] <- gsub("\\\"", "", params[[i]])
  if(i==2){
    params[[i]] <- gsub("exp\\(", "", params[[i]])
  }
  params[[i]] <- as.numeric(params[[i]])
}

params[[1]] <- params[[1]] * 365.25

params[[2]] <- log(2) / params[[2]]
param_names[2] <- "Intrahost doubling time"

param_names[4] <- "Reproductive rate"

mass_convergence_plots <- list()
mass_convergence_plots_zoomed <- list()

# Load in log-likelihood from run on VM
load("mass-10k/llik.RData")
n_iters <- length(llik)

xlabs <- c(
  "Evolution Rate (subs/site/year)",
  "Intrahost Doubling Time (days)",
  "Sampling Probability",
  "Reproductive Number"
)

for (i in 1:length(param_names)) {
  mass_convergence_plots[[i]] <- ggplot(data.frame(x = 1:n_iters, y = params[[i]]), aes(x=x,y=y)) +
    geom_line() +
    xlab("Global Iteration") +
    ggtitle(xlabs[i]) +
    ylab(element_blank()) +
    theme_minimal()

  mass_convergence_plots_zoomed[[i]] <- ggplot(data.frame(x = (n_iters / 2 + 1):n_iters, y = params[[i]][(n_iters / 2 + 1):n_iters]), aes(x=x,y=y)) +
    geom_line() +
    xlab("Global Iteration") +
    ggtitle(xlabs[i]) +
    ylab(element_blank()) +
    theme_minimal()
}




mass_convergence_plots[[5]] <- ggplot2::ggplot() +
  geom_line(mapping = aes(x = 1:n_iters, y = llik)) +
  xlab("Global Iteration") +
  ggtitle("Log Posterior") +
  ylab(element_blank()) +
  theme_minimal()

mass_convergence_plots_zoomed[[5]] <- ggplot2::ggplot() +
  geom_line(mapping = aes(x = (n_iters / 2 + 1):n_iters, y = llik[(n_iters / 2 + 1):n_iters])) +
  xlab("Global Iteration") +
  ggtitle("Log Posterior") +
  ylab(element_blank()) +
  theme_minimal()


cowplot::plot_grid(
  mass_convergence_plots[[5]],
  mass_convergence_plots_zoomed[[5]],
  mass_convergence_plots[[1]],
  mass_convergence_plots_zoomed[[1]],
  mass_convergence_plots[[2]],
  mass_convergence_plots_zoomed[[2]],
  mass_convergence_plots[[3]],
  mass_convergence_plots_zoomed[[3]],
  mass_convergence_plots[[4]],
  mass_convergence_plots_zoomed[[4]],
  ncol = 2,
  labels = "AUTO"
)

ggsave("figs/mass-convergence.pdf", width = 10, height = 12)
ggsave("figs/mass-convergence.png", width = 10, height = 12)

print(coda::effectiveSize(llik[(n_iters / 2 + 1):n_iters]))
for (i in 1:length(param_names)) {
  print(coda::effectiveSize(params[[i]][(n_iters / 2 + 1):n_iters]))
}


### Performance / convergence of h5n1
ss <- juniper0::summarize(res, burnin= 0)
h5n1_params <- list()
h5n1_params[[1]] <- ss$mu * 365.25
h5n1_params[[2]] <- log(2) / ss$N_effs
h5n1_params[[3]] <- ss$pi
h5n1_params[[4]] <- ss$R


h5n1_convergence_plots <- list()
h5n1_convergence_plots_zoomed <- list()

# Load in log-likelihood from run on VM
llik <- res[[1]]
n_iters <- length(llik)



for (i in 1:length(param_names)) {
  h5n1_convergence_plots[[i]] <- ggplot(data.frame(x = 1:n_iters, y = h5n1_params[[i]]), aes(x=x,y=y)) +
    geom_line() +
    xlab("Global Iteration") +
    ggtitle(xlabs[i]) +
    ylab(element_blank()) +
    theme_minimal()

  h5n1_convergence_plots_zoomed[[i]] <- ggplot(data.frame(x = (n_iters * 0.5 + 1):n_iters, y = h5n1_params[[i]][(n_iters * 0.5 + 1):n_iters]), aes(x=x,y=y)) +
    geom_line() +
    xlab("Global Iteration") +
    ggtitle(xlabs[i]) +
    ylab(element_blank()) +
    theme_minimal()
}




h5n1_convergence_plots[[5]] <- ggplot2::ggplot() +
  geom_line(mapping = aes(x = 1:n_iters, y = llik)) +
  xlab("Global Iteration") +
  ggtitle("Log Posterior") +
  ylab(element_blank()) +
  theme_minimal()

h5n1_convergence_plots_zoomed[[5]] <- ggplot2::ggplot() +
  geom_line(mapping = aes(x = (n_iters * 0.5 + 1):n_iters, y = llik[(n_iters * 0.5 + 1):n_iters])) +
  xlab("Global Iteration") +
  ggtitle("Log Posterior") +
  ylab(element_blank()) +
  theme_minimal()


cowplot::plot_grid(
  h5n1_convergence_plots[[5]],
  h5n1_convergence_plots_zoomed[[5]],
  h5n1_convergence_plots[[1]],
  h5n1_convergence_plots_zoomed[[1]],
  h5n1_convergence_plots[[2]],
  h5n1_convergence_plots_zoomed[[2]],
  h5n1_convergence_plots[[3]],
  h5n1_convergence_plots_zoomed[[3]],
  h5n1_convergence_plots[[4]],
  h5n1_convergence_plots_zoomed[[4]],
  ncol = 2,
  labels = "AUTO"
)

ggsave("figs/h5n1-convergence.pdf", width = 10, height = 12)
ggsave("figs/h5n1-convergence.png", width = 10, height = 12)

print(coda::effectiveSize(llik[(n_iters * 0.5 + 1):n_iters]))
for (i in 1:length(param_names)) {
  print(coda::effectiveSize(h5n1_params[[i]][(n_iters * 0.5 + 1):n_iters]))
}

## Parameter postrior densities for h5n1 and sars-cov-2

h5n1_dens_plots <- list()
mass_dens_plots <- list()


param_names <- c(
  "Evolution Rate",
  "Intrahost Doubling Time",
  "Sampling Probability",
  "Reproductive Number"
)

for (i in 1:4) {
  h5n1_dens_plots[[i]] <- ggplot2::ggplot(data.frame(x = h5n1_params[[i]][(n_iters * 0.5 + 1):n_iters])) +
    geom_density(aes(x=x), fill = rgb(0,0,0,0.25), adjust = 2) +
    xlab(xlabs[i]) +
    ylab("Density") +
    ggtitle(paste0(param_names[i], ", H5N1")) +
    theme_minimal()
  mass_dens_plots[[i]] <- ggplot2::ggplot(data.frame(x = params[[i]][(n_iters * 0.5 + 1):n_iters])) +
    geom_density(aes(x=x), fill = rgb(0,0,0,0.25), adjust = 2) +
    xlab(xlabs[i]) +
    ylab("Density") +
    ggtitle(paste0(param_names[i], ", SARS-CoV-2")) +
    theme_minimal()
}

cowplot::plot_grid(
  h5n1_dens_plots[[1]],
  h5n1_dens_plots[[2]],
  h5n1_dens_plots[[3]],
  h5n1_dens_plots[[4]],
  mass_dens_plots[[1]],
  mass_dens_plots[[2]],
  mass_dens_plots[[3]],
  mass_dens_plots[[4]],
  ncol = 2,
  labels = "AUTO",
  align = "hv"
)

ggsave("figs/h5n1-mass-density.pdf", width = 8, height = 10)
ggsave("figs/h5n1-mass-density.png", width = 8, height = 10)

## Inferences

# Evolution rate
mean(h5n1_params[[1]][(n_iters * 0.5 + 1):n_iters])
quantile(h5n1_params[[1]][(n_iters * 0.5 + 1):n_iters], c(0.025, 0.975))

# Reproductive number
mean(h5n1_params[[4]][(n_iters * 0.5 + 1):n_iters])
quantile(h5n1_params[[4]][(n_iters * 0.5 + 1):n_iters], c(0.025, 0.975))

# Transmissions per sampled case, California
mean(h5n1_trans_df$Offspring[h5n1_trans_df$State == "California"])
quantile(h5n1_trans_df$Offspring[h5n1_trans_df$State == "California"], c(0.025, 0.975))

# Number of links w/ posterior probability > 50%
sum(ca_direct > 0.5)






