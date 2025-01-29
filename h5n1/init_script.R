library(juniper0)
set.seed(1)
init <- initialize(
  n_global = 50000,
  init_mu = 2e-5
)
save(init, file = "init.RData")
