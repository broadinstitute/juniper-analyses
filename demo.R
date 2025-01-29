## Comparison and visualization of outbreak reconstruction on simulated datasets
library(devtools)
library(juniper0)
library(simulatR)
library(ggplot2)
library(cowplot)
library(TransPhylo)
library(outbreaker2)

set.seed(4)

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

## Summarize results from online run

plot_params <- function(ideal = T, output_plot = T){

  ## Visualization of each parameter for each experiment
  viz_params <- list()
  params <- c("mu", "pi", "R", "N_eff", "time_of_MRCA")
  param_names <- c("Evolution Rt.", "Sampling Prob.", "Reproductive Nr.", "Intrahost Doubling Time", "Start Date")

  # True parameter values by experiment
  truths <- list()
  truths[["mu"]] <- combos$mu
  truths[["pi"]] <- combos$p_sample
  truths[["R"]] <- combos$R
  truths[["N_eff"]] <- log(2) / combos$N_eff
  truths[["time_of_MRCA"]] <- rep(as.Date("2000-01-01"), nrow(combos))

  errors <- list()

  for (p in params) {
    df <- data.frame(Experiment = integer(0), Values = numeric(0))
    for (i in 1:nrow(combos)) {
      dir <- paste0("experiment_", i)
      load(paste0(dir, "/output_ideal.RData"))
      if(p == "N_eff"){
        df <- rbind(
          df,
          data.frame(Experiment = factor(i), Values = log(2) / out[["N_effs"]])
        )
      }else{
        df <- rbind(
          df,
          data.frame(Experiment = factor(i), Values = out[[p]])
        )
      }

    }

    df2 <- data.frame(Experiment = factor(1:nrow(combos)), Values = truths[[p]])

    err <- c()
    for (j in 1:nrow(combos)) {
      err <- c(
        err,
        as.numeric(df$Values[df$Experiment == j] - df2$Values[j])
      )
    }
    errors[[match(p, params)]] <- err


    viz_params[[match(p, params)]] <- ggplot() +
      geom_boxplot(data = df, mapping = aes(x = Experiment, y = Values, group = Experiment), fill = "green", color = "#008800", outliers = F) +
      geom_point(data = df2, mapping = aes(x = Experiment, y = Values), size = 2, shape = 23, fill = "#1188FF") +
      ylab(param_names[match(p, params)]) +
      theme_minimal()
  }

  if(!output_plot){
    return(errors)
  }

  list(
    viz_params[[1]],
    viz_params[[2]],
    viz_params[[3]],
    viz_params[[4]],
    viz_params[[5]]
  )

}

# Models we're considering
models <- c("ideal", "misspecified", "consensus", "fully_sampled", "transphylo", "outbreaker2", "badtrip")
light_colors <- colorRampPalette(c("green", "red"))(length(models))
dark_colors <- colorRampPalette(c("#008800", "#880000"))(length(models))
pale_colors <- colorRampPalette(c("#BBFFBB", "#FFBBBB"))(length(models))
mid_colors <- colorRampPalette(c("#77FF77", "#FF7777"))(length(models))


# for (i in 1:23) {
#   dir <- paste0("experiment_", i)
#   load(paste0(dir, "/output_ideal.RData"))
#   oo <- out
#   load(paste0(dir, "/output_misspecified.RData"))
#   print(identical(oo$indirect_transmissions, out$indirect_transmissions))
# }


## Visualization of direct transmission reconstruction accuracy
plot_trans <- function(run = "ideal", plot = T){

  df <- data.frame(Experiment = integer(0), Probability = numeric(0))
  df2 <- data.frame(Experiment = integer(0), Accuracy = numeric(0), Coverage = numeric(0), Hit = numeric(0))
  for (i in 1:nrow(combos)) {

    print(i)

    dir <- paste0("experiment_", i)

    if(run == "ideal"){
      load(paste0(dir, "/output_ideal.RData"))
    }else if(run == "misspecified"){
      load(paste0(dir, "/output_misspecified.RData"))
    }else if(run == "consensus"){
      load(paste0(dir, "/output_consensus.RData"))
    }else if(run == "fully_sampled"){
      load(paste0(dir, "/output_fully_sampled.RData"))
      out$direct_transmissions <- out$direct_transmissions[2:51, 2:51]
      out$indirect_transmissions <- out$indirect_transmissions[2:51, 2:51]
    }else if(run == "outbreaker2"){
      load(paste0(dir, "/outbreaker.RData"))

      # Kill burn-in
      out <- out[(nrow(out) * 0.2 + 1):nrow(out), ]

      # Filter to columns relating to ancestry
      kappas <- out[, 109:158]
      out <- out[, 9:58]
      out <- as.matrix(out)

      # NA means root; replace with 0 for now
      out[is.na(out)] <- 0

      # Shift by 1 since 1st column is going to represent root
      out <- out + 1

      adj <- matrix(0, nrow = ncol(out) + 1, ncol = ncol(out) + 1)
      adj_indirect <- matrix(0, nrow = ncol(out) + 1, ncol = ncol(out) + 1)

      out <- out[round(0.2 * nrow(out)):nrow(out), ]

      for (r in 1:nrow(out)) {

        # Which transmissions do we update?
        inds <- cbind(out[r, ], 2:(ncol(out) + 1))

        adj_indirect[inds] <- adj_indirect[inds] + 1

        direct <- kappas[r, ] == 1
        direct[is.na(direct)] <- TRUE

        # Only update direct transmissions
        inds <- inds[direct, ]

        adj[inds] <- adj[inds] + 1
        #print(r)
      }

      adj <- adj / nrow(out)
      adj_indirect <- adj_indirect / nrow(out)



      fasta <- read.FASTA(paste0(dir, "/input_data/aligned.fasta"))
      names <- gsub("\\|.*", "", names(fasta))

      colnames(adj) <- c("reference_genome (root)", names)
      rownames(adj) <- c("reference_genome (root)", names)

      colnames(adj_indirect) <- c("reference_genome (root)", names)
      rownames(adj_indirect) <- c("reference_genome (root)", names)

      # For sake of comparison: delete 1st row and column
      # May bring this back later
      adj <- adj[2:(ncol(out) + 1), 2:(ncol(out) + 1)]
      adj_indirect <- adj_indirect[2:(ncol(out) + 1), 2:(ncol(out) + 1)]

      out <- list()
      out$direct_transmissions <- adj
      out$indirect_transmissions <- adj_indirect

    }else if(run == "badtrip"){

      if(i <= 21){
        net <- readLines(paste0(
          "experiment_", i, "/badtrip/summary_network.txt"
        ))
        # Get number of hosts
        N <- length(unlist(strsplit(net[1], ","))) - 1
        start <- which(net == "Probabilities of direct transmittor to each sampled host: ")

        adj <- matrix(0, nrow = N+1, ncol = N+1)
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

            adj[from + 1, to + 1] <- probs

          }
        }

        if(all(abs(colSums(adj)[2:(N+1)] - 1) < 1e-10)){
          print(paste("No indirect transmissions for experiment", i))
        }

        adj_indirect <- adj

        fasta <- read.FASTA(paste0(dir, "/input_data/aligned.fasta"))
        names <- gsub("\\|.*", "", names(fasta))

        colnames(adj) <- c("reference_genome (root)", names)
        rownames(adj) <- c("reference_genome (root)", names)

        colnames(adj_indirect) <- c("reference_genome (root)", names)
        rownames(adj_indirect) <- c("reference_genome (root)", names)

        # For sake of comparison: delete 1st row and column
        # May bring this back later
        adj <- adj[2:51, 2:51]
        adj_indirect <- adj_indirect[2:51, 2:51]

        out <- list()
        out$direct_transmissions <- adj
        out$indirect_transmissions <- adj_indirect
      }


    }else if(run == "transphylo"){

      fasta <- read.FASTA(paste0(dir, "/input_data/aligned.fasta"))
      names <- gsub("\\|.*", "", names(fasta))
      names <- c("reference_genome (root)", names)

      # Read in transmission network
      load(paste0(
        "experiment_",
        i,
        "/res_TransPhylo_slim.RData"
      ))

      res_TransPhylo <- res_TransPhylo[(0.2 * length(res_TransPhylo) + 1):length(res_TransPhylo)]

      # Names according to transphylo
      tp_names <- res_TransPhylo[[1]]$ctree$nam
      tp_names <- gsub("\\|.*", "", tp_names)
      n <- length(tp_names)
      tp_names <- c("reference_genome (root)", tp_names)

      adj <- matrix(0, ncol = n+1, nrow = n+1)
      colnames(adj) <- tp_names
      rownames(adj) <- tp_names

      adj_indirect <- adj

      # Knock out burnin
      res_TransPhylo <- res_TransPhylo[round(0.2 * length(res_TransPhylo)):length(res_TransPhylo)]

      for (r in 1:length(res_TransPhylo)) {

        ttree=extractTTree(res_TransPhylo[[r]]$ctree)$ttree

        infectors=ttree[1:n,3]
        infecteds=1:n

        infectors_indirect <- infectors
        while (any(infectors_indirect > n)) {
          who <- which(infectors_indirect > n)
          anc <- ttree[infectors_indirect[who], 3]
          infectors_indirect[who] <- anc
        }

        # Offset by 1 to account for root column in adjacency matrix
        adj_indirect[cbind(infectors_indirect+1,infecteds+1)] <- adj_indirect[cbind(infectors_indirect+1,infecteds+1)] + 1/length(res_TransPhylo)

        w=which(infectors==0|infectors>n)
        infectors=infectors[setdiff(1:n,w)]
        infecteds=infecteds[setdiff(1:n,w)]

        adj[cbind(infectors+1,infecteds+1)] <- adj[cbind(infectors+1,infecteds+1)] + 1/length(res_TransPhylo)
      }

      # Match column and row names
      adj <- adj[match(names, tp_names), match(names, tp_names)]
      adj_indirect <- adj_indirect[match(names, tp_names), match(names, tp_names)]

      names <- names[-1]

      adj <- adj[2:51, 2:51]
      adj_indirect <- adj_indirect[2:51, 2:51]

      out <- list()
      out$direct_transmissions <- adj
      out$indirect_transmissions <- adj_indirect

    }

    if(i > 21 & run == "badtrip"){
      df <- rbind(
        df,
        data.frame(
          Experiment = factor(i),
          Probability = NA
        )
      )

      df2 <- rbind(
        df2,
        data.frame(
          Experiment = factor(i),
          Accuracy = NA,
          Coverage = NA,
          Hit = NA,
          Accuracy_Indirect = NA,
          Coverage_Indirect = NA,
          Hit_Indirect = NA
        )
      )
    }else{
      # Read true transmission network
      truth <- read.csv(paste0(dir, "/input_data/transmission.csv"))
      truth$from[truth$from == "person_1"] <- "reference_genome (root)"
      truth$to[truth$to == "person_1"] <- "reference_genome (root)"

      # FOR NOW: deduplicate truth$to
      truth <- truth[!duplicated(truth$to), ]

      # Names of sequences reconstructed
      names <- colnames(out$direct_transmissions)

      # Indirect transmissions
      truth_indirect <- truth
      for (n in names) {
        while(!(truth_indirect$from[which(truth_indirect$to == n)] %in% c(names, "reference_genome (root)"))){
          anc <- truth_indirect$from[which(truth_indirect$to == n)]
          truth_indirect$from[which(truth_indirect$to == n)] <- truth_indirect$from[which(truth_indirect$to == anc)]
        }
      }

      truth <- truth[truth$from %in% names & truth$to %in% names, ]
      truth_indirect <- truth_indirect[truth_indirect$to %in% names, ]

      truth <- as.matrix(truth[, 1:2])
      truth_indirect <- as.matrix(truth_indirect[, 1:2])

      # Bind extra column designating root
      out$indirect_transmissions <- rbind(1 - colSums(out$indirect_transmissions), out$indirect_transmissions)
      out$indirect_transmissions <- cbind(rep(0, nrow(out$indirect_transmissions)), out$indirect_transmissions)
      rownames(out$indirect_transmissions)[1] <- "reference_genome (root)"
      colnames(out$indirect_transmissions)[1] <- "reference_genome (root)"

      df <- rbind(
        df,
        data.frame(
          Experiment = factor(i),
          Probability = out$indirect_transmissions[truth_indirect] #Change this if adjusting stat
        )
      )

      # Accuracy
      # true_adj <- matrix(0, ncol = ncol(out$direct_transmissions), nrow = nrow(out$direct_transmissions))
      # colnames(true_adj) <- colnames(out$direct_transmissions)
      # rownames(true_adj) <- rownames(out$direct_transmissions)
      # true_adj[truth] <- 1
      #
      # true_adj_indirect <- matrix(0, ncol = ncol(out$indirect_transmissions), nrow = nrow(out$indirect_transmissions))
      # colnames(true_adj_indirect) <- colnames(out$indirect_transmissions)
      # rownames(true_adj_indirect) <- rownames(out$indirect_transmissions)
      # true_adj_indirect[truth_indirect] <- 1
      #
      # #accuracy <- sum(out$direct_transmissions[truth]) / sum(out$direct_transmissions)
      # #accuracy_indirect <- sum(out$indirect_transmissions[truth_indirect]) / sum(out$indirect_transmissions)
      # accuracy <- mean(1 - 0.5*colSums(abs(true_adj - out$direct_transmissions)))
      # accuracy_indirect <- mean(1 - 0.5*colSums(abs(true_adj_indirect - out$indirect_transmissions)))


      # Probability that 95% credible set for ancestor contains the true ancestor, or if no such interval exists, probabiltiy of being infected by unobserved ancestor
      accuracy <- numeric(0)
      accuracy_indirect <- numeric(0)

      coverage <- integer(0)
      coverage_indirect <- integer(0)
      hit <- integer(0)
      hit_indirect <- integer(0)
      for (n in names) {
        froms <- out$direct_transmissions[, n]
        froms_indirect <- out$indirect_transmissions[, n]

        froms <- c(froms, 1 - sum(froms))

        names(froms)[length(froms)] <- "unobserved"

        froms <- sort(froms, decreasing = T)
        froms_indirect <- sort(froms_indirect, decreasing = T)

        # Cumulative probability
        prob <- 0
        # Credible set
        set <- character(0)
        k <- 1
        while (prob < 0.95 & k <= length(froms)) {
          set <- c(set, names(froms)[k])
          prob <- prob + froms[k]
          k <- k + 1
        }

        # Cumulative probability - indirect ancestor
        prob <- 0
        # Credible set
        set_indirect <- character(0)
        k <- 1
        while (prob < 0.95 & k <= length(froms_indirect)) {
          set_indirect <- c(set_indirect, names(froms_indirect)[k])
          prob <- prob + froms_indirect[k]
          k <- k + 1
        }

        #print(truth)
        #stop("hheh")

        if(n %in% truth[,2]){

          accuracy <- c(accuracy, froms[which(names(froms) == truth[,1][which(truth[,2] == n)])])

          # The immediate ancestor of n was sampled
          if(truth[,1][which(truth[,2] == n)] %in% set){
            coverage <- c(coverage, 1)
          }else{
            coverage <- c(coverage, 0)
          }

          if(names(froms)[which.max(froms)] == truth[,1][which(truth[,2] == n)]){
            hit <- c(hit, 1)
          }else{
            hit <- c(hit, 0)
          }

        }else{

          accuracy <- c(accuracy, froms[which(names(froms) == "unobserved")])

          # The immediate ancestor of n was not sampled
          if("unobserved" %in% set){
            coverage <- c(coverage, 1)
          }else{
            coverage <- c(coverage, 0)
          }

          if(names(froms)[which.max(froms)] == "unobserved"){
            hit <- c(hit, 1)
          }else{
            hit <- c(hit, 0)
          }

        }

        # Indirect accuracy
        accuracy_indirect <- c(accuracy_indirect, froms_indirect[which(names(froms_indirect) == truth_indirect[,1][which(truth_indirect[,2] == n)])])

        # Indirect coverage
        if(truth_indirect[,1][which(truth_indirect[,2] == n)] %in% set_indirect){
          coverage_indirect <- c(coverage_indirect, 1)
        }else{
          coverage_indirect <- c(coverage_indirect, 0)
        }

        # Indirect hit rate
        if(truth_indirect[,1][which(truth_indirect[,2] == n)] == names(froms_indirect)[which.max(froms_indirect)]){
          hit_indirect <- c(hit_indirect, 1)
        }else{
          hit_indirect <- c(hit_indirect, 0)
        }
      }

      accuracy <- mean(accuracy)
      accuracy_indirect <- mean(accuracy_indirect)

      coverage <- mean(coverage)
      hit <- mean(hit)

      coverage_indirect <- mean(coverage_indirect)
      hit_indirect <- mean(hit_indirect)

      df2 <- rbind(
        df2,
        data.frame(
          Experiment = factor(i),
          Accuracy = accuracy,
          Coverage = coverage,
          Hit = hit,
          Accuracy_Indirect = accuracy_indirect,
          Coverage_Indirect = coverage_indirect,
          Hit_Indirect = hit_indirect
        )
      )
    }
  }

  if(!plot){
    return(df2)
  }

  out <- list(ggplot() +
    geom_boxplot(
      data = df,
      mapping = aes(x = Experiment, y = Probability, group = Experiment),
      outliers = F,
      fill = pale_colors[match(run, models)],
      color = mid_colors[match(run, models)]
    ) +
    geom_point(
      data = df2,
      mapping = aes(x = Experiment, y = Accuracy_Indirect, group = Experiment),
      shape = 1,
      size = 2,
      stroke = 0.75,
      color = "black"
    ) +
    geom_point(
      data = df2,
      mapping = aes(x = Experiment, y = Coverage_Indirect, group = Experiment),
      shape = 0,
      size = 2,
      stroke = 0.75,
      color = "black"
    ) +
    geom_point(
      data = df2,
      mapping = aes(x = Experiment, y = Hit_Indirect, group = Experiment),
      shape = 2,
      size = 2,
      stroke = 0.75,
      color = "black"
    ) +
    theme_minimal())
}



plots <- c(
  plot_params(),
  plot_trans(models[1]),
  plot_trans(models[2]),
  plot_trans(models[3]),
  plot_trans(models[4]),
  plot_trans(models[5]),
  plot_trans(models[6]),
  plot_trans(models[7])
)

comp_df <- data.frame()
for (m in models) {
  tmp <- cbind(plot_trans(run = m, plot = F), m)
  colnames(tmp)[ncol(tmp)] <- "Model"
  comp_df <- rbind(comp_df, tmp)
}

comp_df$Model <- factor(comp_df$Model, levels = models)

# Direct transmission accuracy
for (m in models) {
  print(mean(comp_df$Accuracy[comp_df$Model == m], na.rm = T))
}

# Indirect transmission accuracy
for (m in models) {
  print(mean(comp_df$Accuracy_Indirect[comp_df$Model == m], na.rm = T))
}

# Direct transmission percent change
for (m in models[2:length(models)]) {
  print(
    100 * (mean(comp_df$Accuracy[comp_df$Model == "ideal"], na.rm = T) - mean(comp_df$Accuracy[comp_df$Model == m], na.rm = T)) /
      mean(comp_df$Accuracy[comp_df$Model == m], na.rm = T)
  )
}

# Indirect transmission percent change
for (m in models[2:length(models)]) {
  print(
    100 * (mean(comp_df$Accuracy_Indirect[comp_df$Model == "ideal"], na.rm = T) - mean(comp_df$Accuracy_Indirect[comp_df$Model == m], na.rm = T)) /
      mean(comp_df$Accuracy_Indirect[comp_df$Model == m], na.rm = T)
  )
}

# Direct donor set accuracy
for (m in models) {
  print(mean(comp_df$Coverage[comp_df$Model == m], na.rm = T))
}

# Indirect donor set accuracy
for (m in models) {
  print(mean(comp_df$Coverage_Indirect[comp_df$Model == m], na.rm = T))
}

# Direct DSA percent change
for (m in models[2:length(models)]) {
  print(
    100 * (mean(comp_df$Coverage[comp_df$Model == "ideal"], na.rm = T) - mean(comp_df$Coverage[comp_df$Model == m], na.rm = T)) /
      mean(comp_df$Coverage[comp_df$Model == m], na.rm = T)
  )
}

# Indirect DSA percent change
for (m in models[2:length(models)]) {
  print(
    100 * (mean(comp_df$Coverage_Indirect[comp_df$Model == "ideal"], na.rm = T) - mean(comp_df$Coverage_Indirect[comp_df$Model == m], na.rm = T)) /
      mean(comp_df$Coverage_Indirect[comp_df$Model == m], na.rm = T)
  )
}

# Direct most likely donor accuracy
for (m in models) {
  print(mean(comp_df$Hit[comp_df$Model == m], na.rm = T))
}

# Indirect most likely donor set accuracy
for (m in models) {
  print(mean(comp_df$Hit_Indirect[comp_df$Model == m], na.rm = T))
}



comp_df1 <- comp_df[, c(1:4, 8)]
comp_df2 <- comp_df[, c(1, 5:8)]

colnames(comp_df1)[2:4] <- c("Donor-Recipient Accuracy", "Donor Set Accuracy", "Most-Likely Donor Accuracy")
colnames(comp_df2)[2:4] <- c("Donor-Recipient Accuracy", "Donor Set Accuracy", "Most-Likely Donor Accuracy")

comp_df1 <- reshape2::melt(comp_df1)
comp_df2 <- reshape2::melt(comp_df2)

model_names <- c(
  "Targeted JUNIPER",
  "Default JUNIPER",
  "Consensus JUNIPER",
  "No-Intermediates JUNIPER",
  "TransPhylo + IQTree",
  "outbreaker2",
  "BadTrIP"
)

plots[[length(plots) + 1]] <- ggplot(data = comp_df2, mapping = aes(x = variable, y = value)) +
  geom_boxplot(aes(x = variable, group = interaction(variable, Model), color = Model, fill = Model)) +
  xlab("Metric") +
  ylab("Value") +
  scale_color_manual(values = dark_colors, labels = model_names) +
  scale_fill_manual(values = light_colors, labels = model_names) +
  ggtitle("Indirect Transmissions") +
  xlab(element_blank()) +
  theme_minimal() +
  guides(color=guide_legend(ncol=length(models))) +
  theme(legend.title=element_blank())

legend1 <- cowplot::get_legend(plots[[length(plots)]])
plots[[length(plots)]] <- plots[[length(plots)]] + theme(legend.position = "none")

plots[[length(plots) + 1]] <- ggplot(data = comp_df1, mapping = aes(x = variable, y = value)) +
  geom_boxplot(aes(x = variable, group = interaction(variable, Model), color = Model, fill = Model)) +
  xlab("Metric") +
  ylab("Value") +
  scale_color_manual(values = dark_colors) +
  scale_fill_manual(values = light_colors) +
  ggtitle("Direct Transmissions") +
  xlab(element_blank()) +
  theme_minimal() +
  theme(legend.position="none")

o2_dens <- function(param, plot = TRUE){
  df <- data.frame()
  for (i in 1:nrow(combos)) {
    dir <- paste0("experiment_", i)
    load(paste0(dir, "/outbreaker.RData"))

    # Cut out burnin
    out <- out[(0.2*nrow(out) + 1):nrow(out), ]

    df <- rbind(df, data.frame(
      Experiment = factor(i),
      Values = out[, param]
    ))

  }

  truths <- list()
  truths[["mu"]] <- combos$mu
  truths[["pi"]] <- combos$p_sample


  df2 <- data.frame(Experiment = factor(1:nrow(combos)), Values = truths[[param]])

  if(!plot){
    err <- c()
    for (j in 1:nrow(combos)) {
      err <- c(
        err,
        as.numeric(df$Values[df$Experiment == j] - df2$Values[j])
      )
    }
    return(err)
  }


  ggplot() +
    geom_boxplot(data = df, mapping = aes(x = Experiment, y = Values, group = Experiment), fill = light_colors[match("outbreaker2", models)], color = dark_colors[match("outbreaker2", models)], outliers = F) +
    geom_point(data = df2, mapping = aes(x = Experiment, y = Values), size = 2, shape = 23, fill = "#1188FF") +
    ylab(ifelse(param == "mu", "Mutation Rate", "Sampling Rate")) +
    theme_minimal()



}

# Parameter density for transphylo
tp_dens <- function(param, plot = TRUE){
  df <- data.frame()
  for (i in 1:nrow(combos)) {
    dir <- paste0("experiment_", i)
    load(paste0(dir, "/res_TransPhylo_slim.RData"))

    # Cut out burnin
    res_TransPhylo <- res_TransPhylo[(0.2*length(res_TransPhylo) + 1):length(res_TransPhylo)]

    vals <- rep(0, length(res_TransPhylo))
    for (r in 1:length(res_TransPhylo)) {
      if(param == "pi"){
        vals[r] <- res_TransPhylo[[r]]$pi
      }
      if(param == "R"){
        vals[r] <- res_TransPhylo[[r]]$off.r * res_TransPhylo[[r]]$off.p / (1 - res_TransPhylo[[r]]$off.p)
      }
    }

    df <- rbind(df, data.frame(
      Experiment = factor(i),
      Values = vals
    ))

  }

  truths <- list()
  truths[["R"]] <- combos$R
  truths[["pi"]] <- combos$p_sample


  df2 <- data.frame(Experiment = factor(1:nrow(combos)), Values = truths[[param]])

  if(!plot){
    err <- c()
    for (j in 1:nrow(combos)) {
      err <- c(
        err,
        as.numeric(df$Values[df$Experiment == j] - df2$Values[j])
      )
    }
    return(err)
  }


  ggplot() +
    geom_boxplot(data = df, mapping = aes(x = Experiment, y = Values, group = Experiment), fill = light_colors[match("transphylo", models)], color = dark_colors[match("transphylo", models)], outliers = F) +
    geom_point(data = df2, mapping = aes(x = Experiment, y = Values), size = 2, shape = 23, fill = "#1188FF") +
    ylab(ifelse(param == "R", "Reproductive Number", "Sampling Rate")) +
    theme_minimal()



}

# Parameter density for badtrip
bt_dens <- function(param, plot = TRUE){
  df <- data.frame()
  for (i in 1:nrow(combos)) {
    if(i <= 21){
      bt_log <- read.table(paste0("experiment_", i, "/badtrip/BADTRIP_setup.log"), sep = "\t", header = T)
      bt_log <- bt_log[(0.2*nrow(bt_log) + 1):nrow(bt_log), ]

      vals <- bt_log$PoMo.mutRates

      df <- rbind(df, data.frame(
        Experiment = factor(i),
        Values = vals
      ))
    }else{
      df <- rbind(df, data.frame(
        Experiment = factor(i),
        Values = NA
      ))
    }
  }




  df2 <- data.frame(Experiment = factor(1:nrow(combos)), Values = combos$mu)

  if(!plot){
    err <- c()
    for (j in 1:nrow(combos)) {
      err <- c(
        err,
        as.numeric(df$Values[df$Experiment == j] - df2$Values[j])
      )
    }
    return(err)
  }


  ggplot() +
    geom_boxplot(data = df, mapping = aes(x = Experiment, y = Values, group = Experiment), fill = light_colors[match("badtrip", models)], color = dark_colors[match("badtrip", models)], outliers = F) +
    geom_point(data = df2, mapping = aes(x = Experiment, y = Values), size = 2, shape = 23, fill = "#1188FF") +
    ylab("Mutation Rate") +
    theme_minimal()



}


plots[[length(plots) + 1]] <- o2_dens("mu")
plots[[length(plots) + 1]] <- o2_dens("pi")

plots[[length(plots) + 1]] <- tp_dens("R")
plots[[length(plots) + 1]] <- tp_dens("pi")

plots[[length(plots) + 1]] <- bt_dens("mu")

### Reliability of parameter estimates aggregated across experiments
err_juniper <- plot_params(T, F)

err_mu_o2 <- o2_dens("mu", F)
err_pi_o2 <- o2_dens("pi", F)

err_pi_tp <- tp_dens("pi", F)
err_R_tp <- tp_dens("R", F)

err_mu_bt <- bt_dens("mu", F)


mu_err <- data.frame(
  Model = factor(c(
    rep("JUNIPER", length(err_juniper[[1]])),
    rep("outbreaker2", length(err_mu_o2)),
    rep("BadTrIP", length(err_mu_bt))
  ), levels = c("JUNIPER", "BadTrIP", "outbreaker2")),
  Error = c(err_juniper[[1]], err_mu_o2, err_mu_bt)
)

pi_err <- data.frame(
  Model = factor(
    c(rep("JUNIPER", length(err_juniper[[2]])), rep("outbreaker2", length(err_pi_o2)), rep("TransPhylo", length(err_pi_tp))),
    levels = c("JUNIPER", "TransPhylo", "outbreaker2")
  ),
  Error = c(err_juniper[[2]], err_pi_o2, err_pi_tp)
)

R_err <- data.frame(
  Model = c(rep("JUNIPER", length(err_juniper[[3]])), rep("TransPhylo", length(err_R_tp))),
  Error = c(err_juniper[[3]], err_R_tp)
)

plots[[length(plots) + 1]] <- ggplot(mu_err, aes(x = Model, y= log10(abs(Error)), fill = Model, color = Model)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = dark_colors[match(c("ideal", "badtrip", "outbreaker2"), models)]) +
  scale_fill_manual(values = light_colors[match(c("ideal", "badtrip", "outbreaker2"), models)]) +
  ggtitle("Mutation Rate") +
  xlab(element_blank()) +
  ylab("Log10 Absolute Error") +
  theme_minimal() +
  theme(legend.position = "none")

mu_rmse_juniper <- sqrt(mean(mu_err$Error[mu_err$Model == "JUNIPER"]^2))
print(mu_rmse_juniper)
print(
  100 * (sqrt(mean(mu_err$Error[mu_err$Model == "BadTrIP"]^2, na.rm = T)) - mu_rmse_juniper) / mu_rmse_juniper
)
print(
  100 * (sqrt(mean(mu_err$Error[mu_err$Model == "outbreaker2"]^2)) - mu_rmse_juniper) / mu_rmse_juniper
)

plots[[length(plots) + 1]] <- ggplot(pi_err, aes(x = Model,y = log10(abs(Error)), fill = Model, color = Model)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = dark_colors[match(c("ideal", "transphylo", "outbreaker2"), models)]) +
  scale_fill_manual(values = light_colors[match(c("ideal", "transphylo", "outbreaker2"), models)]) +
  ggtitle("Sampling Probability") +
  xlab(element_blank()) +
  ylab("Log10 Absolute Error") +
  theme_minimal() +
  theme(legend.position = "none")

pi_rmse_juniper <- sqrt(mean(pi_err$Error[mu_err$Model == "JUNIPER"]^2))
print(pi_rmse_juniper)
print(
  100 * (sqrt(mean(pi_err$Error[pi_err$Model == "TransPhylo"]^2)) - pi_rmse_juniper) / pi_rmse_juniper
)
print(
  100 * (sqrt(mean(pi_err$Error[pi_err$Model == "outbreaker2"]^2, na.rm = T)) - pi_rmse_juniper) / pi_rmse_juniper
)

plots[[length(plots) + 1]] <- ggplot(R_err, aes(x = Model, y= log10(abs(Error)), fill = Model, color = Model)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = dark_colors[match(c("ideal", "transphylo"), models)]) +
  scale_fill_manual(values = light_colors[match(c("ideal", "transphylo"), models)]) +
  ggtitle("Reproductive Number") +
  xlab(element_blank()) +
  ylab("Log10 Absolute Error") +
  theme_minimal() +
  theme(legend.position = "none")

R_rmse_juniper <- sqrt(mean(R_err$Error[R_err$Model == "JUNIPER"]^2))
print(R_rmse_juniper)
print(
  100 * (sqrt(mean(R_err$Error[R_err$Model == "TransPhylo"]^2, na.rm = T)) - R_rmse_juniper) / R_rmse_juniper
)

## Plot of runtime
starts <- list()
ends <- list()

if(FALSE){
  for (i in 1:23) {
    set.seed(1)
    init <- initialize(
      n_subtrees = 1, # We will parallelize over simulations, not within each one
      n_global = 100, # Fewer reps here, can multiply
      a_g = combos$mu_g[i]^2 / combos$var_g[i],
      lambda_g = combos$mu_g[i] / combos$var_g[i],
      a_s = combos$mu_s[i]^2 / combos$var_s[i],
      lambda_s = combos$mu_s[i] / combos$var_s[i],
      psi = combos$psi[i],
      init_mu = combos$mu[i],
      init_N_eff = combos$N_eff[i],
      indir = paste0("experiment_", i, "/input_data")
    )
    starts[[i]] <- Sys.time()
    res <- run_mcmc(init)
    ends[[i]] <- Sys.time()
    print(i)
  }

  # Get duration of each run. Going to run for 1/100 of the time and then multiply through
  durs <- c()
  for (i in 1:23) {
    durs[i] <- difftime(ends[[i]], starts[[i]], units = "hours") * 100
  }

  #save(durs, file = "~/Desktop/durs.RData")
}
load("runtime/durs.RData")



# Get ESS with 20% burnin
juniper_ess <- c()
for (i in 1:23) {

  load(paste0("experiment_", i, "/output_ideal.RData"))

  llik <- out[[1]]
  llik <- llik[(0.2 * length(llik) + 1):length(llik)]

  juniper_ess[i] <- coda::effectiveSize(llik)
}

# Get BadTrIP ESS
badtrip_ess <- c()
badtrip_dur <- c()

for (i in 1:23) {

  out <- read.table(paste0("experiment_", i, "/badtrip/BADTRIP_setup.log"), sep = "\t", header = T)

  llik <- out$posterior
  llik <- llik[(0.2 * length(llik) + 1):length(llik)]

  badtrip_ess[i] <- coda::effectiveSize(llik)

  load(paste0("experiment_", i, "/badtrip/start.RData"))
  load(paste0("experiment_", i, "/badtrip/end.RData"))
  badtrip_dur[i] <- as.numeric(difftime(end, start, units = "hours"))

}

# Get outbreaker2 ESS and runtime
if(FALSE){
  o2_ess <- c()
  o2_dur <- c()

  for (i in 1:23) {

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

    start <- Sys.time()
    set.seed(1)
    out <- outbreaker(
      data = dat,
      config = list(n_iter = 1000, sample_every = 10, pb = T)
    )
    end <- Sys.time()

    llik <- out$post
    llik <- llik[(0.2 * length(llik) + 1):length(llik)]

    o2_ess[i] <- coda::effectiveSize(llik)
    o2_dur[i] <- as.numeric(difftime(end, start, units = "hours"))

    save(o2_dur, file = "runtime/o2_dur.RData")
    save(o2_ess, file = "runtime/o2_ess.RData")

    print(i)

  }
}
if(TRUE){
  # Get TransPhylo ESS and runtime
  tp_ess <- c()
  tp_dur <- c()

  for (i in 1:23) {

    load(paste0("experiment_", i, "/start_tp.RData"))
    load(paste0("experiment_", i, "/end_tp.RData"))

    load(paste0("experiment_", i, "/res_TransPhylo_slim.RData"))

    llik <- c()
    for (j in 1:1000) {
      llik <- c(llik, res_TransPhylo[[j]]$pTTree + res_TransPhylo[[j]]$pPTree)
    }

    llik <- llik[(0.2 * length(llik) + 1):length(llik)]

    tp_ess[i] <- coda::effectiveSize(llik)
    tp_dur[i] <- as.numeric(difftime(end, start, units = "hours"))

    print(i)

  }

  save(tp_dur, file = "runtime/tp_dur.RData")
  save(tp_ess, file = "runtime/tp_ess.RData")


}

load("runtime/tp_dur.RData")
load("runtime/tp_ess.RData")

load("runtime/o2_dur.RData")
load("runtime/o2_ess.RData")




# 22 and 23 we just want to report expected runtime per ESS
badtrip_ess[22:23] <- NA

# 22 and 23 were run to 750 iters; rest to 100000
e_dur <- badtrip_dur[22:23] * 100000/750
e_dur_per_ess <- mean(e_dur) / mean(badtrip_ess[1:21])
print(e_dur_per_ess * 100 / 24) # Convert to days per 100 ESS

# Mean time per 100 ESS
mean(durs * 100 / juniper_ess)
mean(badtrip_dur * 100 / badtrip_ess, na.rm = T)
mean(o2_dur * 100 / o2_ess)
mean(tp_dur * 100 / tp_ess)


# Plot runtime per 100 ESS
dur_df <- data.frame(
  Experiment = as.factor(1:23),
  Duration = c(o2_dur * 100 / o2_ess, tp_dur * 100 / tp_ess, durs * 100 / juniper_ess, badtrip_dur * 100 / badtrip_ess),
  Method = factor(rep(c("outbreaker2", "TransPhylo", "JUNIPER", "BadTrIP"), each = 23), levels = c("outbreaker2", "JUNIPER", "BadTrIP", "TransPhylo"))
)

plots[[length(plots) + 1]] <- ggplot(dur_df, aes(x = Experiment, y = Duration, group = Method, color = Method, fill = Method)) +
  geom_point(shape = 21, size = 2) +
  ylab("Duration (hours)") +
  ggtitle("Runtime for 100 ESS") +
  scale_color_manual(values = dark_colors[match(c("outbreaker2", "ideal", "badtrip", "transphylo"), models)]) +
  scale_fill_manual(values = light_colors[match(c("outbreaker2", "ideal", "badtrip", "transphylo"), models)]) +
  theme_minimal() +
  scale_y_continuous(trans = "log10") +
  theme(legend.position = "none")

## Make a dummy plot to extract legend
dummy <- ggplot(data.frame(x = 1:3, y = 1:3, Metric = factor(c("Donor-Recipient Accuracy", "Donor Set Accuracy", "Most-Likely Donor Accuracy"), levels = c("Donor-Recipient Accuracy", "Donor Set Accuracy", "Most-Likely Donor Accuracy")))) +
  geom_point(aes(x = 1:3, y = 1:3, shape = Metric)) +
  scale_shape_manual(values = c(1, 0, 2)) +
  theme_minimal() +
  guides(shape=guide_legend(ncol=3)) +
  theme(legend.title=element_blank())

legend2 <- cowplot::get_legend(dummy)

cowplot::plot_grid(
  cowplot::plot_grid(
    plots[[1]] + ggtitle("JUNIPER") + ylab("Evolution Rt. (subs/site/day)"),
    plots[[19]] + ggtitle("BadTrIP") + ylab("Evolution Rt. (subs/site/day)"),
    plots[[15]] + ggtitle("outbreaker2") + ylab("Evolution Rt. (subs/site/day)"),
    plots[[2]] + ggtitle("JUNIPER") + ylab("Sampling Prob."),
    plots[[18]] + ggtitle("TransPhylo") + ylab("Sampling Prob."),
    plots[[16]] + ggtitle("outbreaker2") + ylab("Sampling Prob."),
    plots[[3]] + ggtitle("JUNIPER") + ylab("Reproductive Nr."),
    plots[[17]] + ggtitle("TransPhylo") + ylab("Reproductive Nr."),
    plots[[4]] + ggtitle("JUNIPER") + ylab("Doubling Time (days)"),
    plots[[5]] + ggtitle("JUNIPER"),
    plots[[23]] + ggtitle("Runtime for 100 ESS, All Methods"),
    plots[[6]] + ylab("Accuracy") + ggtitle("Targeted JUNIPER"),
    plots[[7]] + ylab("Accuracy") + ggtitle("Default JUNIPER"),
    plots[[8]] + ylab("Accuracy") + ggtitle("Consensus JUNIPER"),
    plots[[9]] + ylab("Accuracy") + ggtitle("No-Intermediates JUNIPER"),
    plots[[10]] + ylab("Accuracy") + ggtitle("TransPhylo"),
    plots[[11]] + ylab("Accuracy") + ggtitle("outbreaker2"),
    plots[[12]] + ylab("Accuracy") + ggtitle("BadTrIP"),
    align = "hv",
    labels = "AUTO",
    ncol = 3
  ),
  legend1,
  legend2,
  rel_heights = c(30, 1, 1),
  ncol = 1
)

ggsave("figs/individual.pdf", width = 12, height = 16)
ggsave("figs/individual.png", width = 12, height = 16)

## Boxplot of runtimes over methods
plots[[length(plots) + 1]] <- ggplot(dur_df, aes(x = Method, y = Duration, group = Method, color = Method, fill = Method)) +
  geom_boxplot() +
  scale_color_manual(values = dark_colors[match(c("outbreaker2", "ideal", "badtrip", "transphylo"), models)]) +
  scale_fill_manual(values = light_colors[match(c("outbreaker2", "ideal", "badtrip", "transphylo"), models)]) +
  scale_y_continuous(trans = "log10") +
  ylab("Duration (hours)") +
  xlab(element_blank()) +
  ggtitle("Runtime, 100 ESS") +
  theme_minimal() +
  theme(legend.position = "none")



cowplot::plot_grid(
  plots[[14]] + ylab("Accuracy"),
  plots[[13]] + ylab("Accuracy"),

  cowplot::plot_grid(
    plots[[20]] + ggtitle("Evolution Rt. (subs/site/day)"), plots[[21]] + ggtitle("Sampling Prob."), plots[[22]] + ggtitle("Reproductive Nr."), plots[[24]],
    nrow = 1,
    labels = c("C", "D", "E", "F"),
    rel_widths = c(3,3,2,4)
  ),
  legend1,
  ncol = 1,
  labels = c("A", "B", "", ""),
  rel_heights = c(5,5,5, 1)
)


ggsave("figs/aggregate.pdf", width = 11, height = 11)
ggsave("figs/aggregate.png", width = 11, height = 11)






