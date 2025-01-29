## Turn h5n1 run output into reproductive number by state

library(ape)
Rcpp::sourceCpp("cpp_subroutines.cpp")

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

# Juniper output
load("res.RData")
n_iters <- length(res[[1]])
n_obs <- 1519

load("meta_subsampled.RData")
cons <- read.FASTA("aligned.fasta")
names <- gsub("\\|.*", "", names(cons))

# State of each case
case_state <- sub("^[^/]*/[^/]*/", "", meta$strain)
case_state <- gsub("\\/.*", "", case_state)
abbreviated <- which(case_state %in% state.abb)
unspaced_state <- gsub(" ", "", state.name)
case_state[abbreviated] <- unspaced_state[match(case_state[abbreviated], state.abb)]
case_state <- state.name[match(case_state, unspaced_state)]


# Objective function for finding R by state
obj <- function(R, ts, pii, n_kids){
  min_t <- min(ts)
  wbar0 <- rev(wbar(min_t - 1, 0, R * 0.5 / (1 - 0.5), 1 - 0.5, pii, 5, 1, 5, 1, 0.1))
  ws <- wbar0[round(-ts/0.1)] # Still log scale
  out <- 0
  for (j in 1:length(ts)) {
    out <- out + alpha(n_kids[j], 0.5, R, ws[j]) # Takes in wbar on log scale, spits out also log scale
  }
  return(out)
}

big_states <- table(case_state)
big_states <- names(big_states[which(big_states >= 50)])

seqs <- list()
Rs <- c()
piis <- c()
ns <- c()
hs <- list()

for (i in seq(n_iters / 500, n_iters / 2, n_iters / 500)) {
  seqs[[length(seqs) + 1]] <- res[[2]][[i+n_iters/2]]$seq
  Rs[length(Rs) + 1] <- res[[2]][[i+n_iters/2]]$R
  piis[length(piis) + 1] <- res[[2]][[i+n_iters/2]]$pi
  ns[length(ns) + 1] <- res[[2]][[i+n_iters/2]]$n
  hs[[length(hs) + 1]] <- res[[2]][[i+n_iters/2]]$h
}

## Get mean trans by state for ith iteration
get_trans <- function(seq, R, pii, n, h){
  # Unpack results

  ts <- unlist(seq)

  #min_t <- min(ts)
  #wbar0 <- rev(wbar(min_t - 1, 0, R * 0.5 / (1 - 0.5), 1 - 0.5, pii, 5, 1, 5, 1, 0.1))

  #rho <- R
  psi <- 0.5

  # Get state by case through parsimony ancestral state reconstruction
  state <- rep(list(character(0)), n)
  for (j in 2:(n_obs + 1)) {
    state[[j]] <- case_state[j-1]
  }

  # Set of kids for branch point / sampled cases
  kids <- rep(list(integer(0)), n)
  for (j in 2:n) {
    kids[[h[j]]] <- c(kids[[h[j]]], j)
  }

  # FORWARD BFS order
  ord <- bfs(h)

  # Loop over nodes in reverse BFS order
  for (j in rev(ord)) {
    # If you have kids...
    if(length(kids[[j]]) > 0){
      # Table of states of all of the kids, flattened
      kid_states <- table(unlist(state[kids[[j]]]))

      if(length(kid_states) == 0){
        state[[j]] <- NA
      }else{
        # The valid states at j are the ones that maximize the number of kid states
        state[[j]] <- names(kid_states)[kid_states == max(kid_states)]
      }
    }
  }

  # Now, pick the state
  if(!all(is.na(state[[1]]))){
    state[[1]] <- sample(state[[1]], 1)
  }

  for (j in ord[-1]) {

    if(all(is.na(state[[j]]))){
      if(!is.na(state[[h[j]]])){
        state[[j]] <- state[[h[j]]]
      }
    }else{

      if(state[[h[j]]] %in% state[[j]]){
        state[[j]] <- state[[h[j]]]
      }else{
        state[[j]] <- sample(state[[j]], 1)
      }
    }
  }

  # Flatten resulting list
  state <- unlist(state)

  # Number of hosts along each edge
  n_hosts_edge <- sapply(seq, length)

  # Update state for ALL hosts
  state <- unlist(mapply(rep, state, n_hosts_edge, USE.NAMES = F, SIMPLIFY = F))

  # Number of kids for ALL hosts
  n_kids <- sapply(kids, length)
  n_kids <- unlist(mapply(function(d, n){c(d, rep(1, n-1))}, n_kids, n_hosts_edge, USE.NAMES = F, SIMPLIFY = F))

  n_tot <- length(n_kids)

  # Get correct wbar for each j
  #ws <- wbar0[round(-ts/0.1)] # Still log scale

  # Total number of hosts


  # # Normalizing constant
  # norms <- rep(0, n_tot)
  # for (j in 1:n_tot) {
  #   norms[j] <- exp(alpha(n_kids[j], psi, rho, ws[j])) # Takes in wbar on log scale
  # }
  #
  # # Sum_{k=d}^\infty (k choose d) alpha(k) wbar^(k-d) * k
  # means <- rep(0, n_tot)
  # for (j in 1:n_tot) {
  #   means[j] <- exp(alpha2(n_kids[j], psi, rho, ws[j])) / norms[j] # Takes in wbar on log scale
  # }

  ## Estimation approach to finding rho for each state

  out <- rep(0, length(big_states))

  for (k in 1:length(big_states)) {
    ts_state <- ts[state == big_states[k]]
    n_kids_state <- n_kids[state == big_states[k]]

    optimize(function(R){-obj(R, ts_state, pii, n_kids_state)}, lower = 0, upper = 5)



    out[k] <- optimize(function(R){-obj(R, ts_state, pii, n_kids_state)}, lower = 0, upper = 5)$minimum
  }

  return(out)
}

mean_trans_state <- list()

for (i in 1:length(Rs)) {
  mean_trans_state[[i]] <- get_trans(seqs[[i]], Rs[i], piis[i], ns[i], hs[[i]])
  print(i)
}
save(mean_trans_state, file = "mean_trans_state.RData")


