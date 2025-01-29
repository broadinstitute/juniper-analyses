## Fit parameters of denovo iSNV model

vars <- read.csv("isnv-frequency-model/allvariants.tsv", sep = "\t")

assembly <- read.csv("isnv-frequency-model/assembly.txt", sep = "\t")

# Subset to relevant BioProject
assembly <- assembly[assembly$bioproject_accession == "PRJNA715749", ]

vars$sample_name <- gsub("\\..*", "", vars$sample_name)
vars <- vars[vars$sample_name %in% assembly$sample_sanitized, ]

af <- vars$INFO
af <- gsub(".*;AF=", "", af)
af <- gsub(";SB=.*", "", af)
af <- as.numeric(af)

dp <- vars$INFO
dp <- gsub(".*DP=", "", dp)
dp <- gsub(";AF=.*", "", dp)
dp <- as.numeric(dp)

sb <- vars$INFO
sb <- gsub(".*;SB=", "", sb)
sb <- gsub(";DP4=.*", "", sb)
sb <- as.numeric(sb)

pos <- vars$POS
problematic <- read.csv("mass-10k/input_data/problematic.csv")

filtered_af <- af[af > 0.03 & af < 0.97 & dp > 100 & sb <= 10 & !(pos %in% problematic$x)]
filtered_af <- pmin(filtered_af, 1- filtered_af)

maf <- filtered_af
save(maf, file = "isnv-frequency-model/maf.RData")

n_seqs <- length(unique(vars$sample_name))

# Total number of observed positions, approx.
n_obs <- n_seqs * 29903

# Number of observations below threshhold
n_below <- n_obs - length(filtered_af)

# Fraction of obsrvatoions below threshhold
p_below <- n_below / n_obs
p_above <- length(filtered_af) / n_obs

# Estimate of r based on number of sites with no observed iSNV
r <- ((-1+p_below) *0.03)/(p_below* (-1+0.03))

# Generate plots


min_af_filter <- 0.03
max_af_filter <- 0.5

# PDF, CDF, and quantile funtions, using approximation in Methods
approx_pdf <- function(x){
  max_af_filter * min_af_filter / ((max_af_filter - min_af_filter) * x^2)
}

approx_cdf <- function(x){
  (max_af_filter * (min_af_filter - x))/((min_af_filter - max_af_filter) * x)
}

approx_qtile <- function(p){
  min_af_filter * max_af_filter / (max_af_filter + p*(min_af_filter - max_af_filter))
}

df <- data.frame(log(maf))

p1 <- ggplot(df, aes(x = (maf))) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.01, boundary = 0.01, color = "white", fill = "#CCCCCC") +
  geom_vline(xintercept=min_af_filter + 0.00001, linetype=3, color = "#BB5522") +
  geom_function(fun = approx_pdf, color = "#2255BB", linewidth = 1.5, linetype = "dashed") +
  xlab("Minor Allele Frequency") +
  ylab("Probability Density") +
  theme_minimal()

print(p1)


### Histogram with log-transform

round_label <- function(x){
  round(x, 2)
}

# Empirical log density
bw <- 0.01
rights <- seq(min_af_filter + bw, max_af_filter, bw)
lefts <- rights - bw
counts <- c()
for (i in 1:length(rights)) {
  counts[i] <- sum(maf < rights[i] & maf >= lefts[i])
}
counts <- counts/sum(counts)/bw

p2 <- ggplot(data.frame(x = lefts), aes(x=x)) +
  geom_line(aes(y = counts), stat = 'identity', linewidth = 1.5, color = "#CCCCCC") +
  geom_vline(xintercept=min_af_filter + 0.00001, linetype=3, color = "#BB5522") +
  geom_function(fun = approx_pdf, color = "#2255BB", linewidth = 1.5, linetype = "dashed") +
  xlab("Log Minor Allele Frequency") +
  scale_x_continuous(trans = 'log', labels = round_label) +
  scale_y_continuous(trans = 'log', labels = round_label) +
  ylab("Log Probability Density") +
  theme_minimal()

print(p2)

### Q-Q Plot

breaks <- approx_cdf(
  seq(min_af_filter, max_af_filter, 0.01)
)

# Rewrite CDF to correctly compute values below min_af_filter and above max_af_filter

full_cdf <- function(x){
  below <- which(x <= min_af_filter)
  above <- which(x  > max_af_filter)
  other <- which(x > min_af_filter & x <= max_af_filter)

  x[below] <- 0
  x[above] <- 1
  x[other] <- approx_cdf(x[other])

  x
}



### CDF comparison
p3 <- ggplot(data.frame(x=maf), aes(x=x)) +
  stat_ecdf(geom = "line", color = "#BB5522", linewidth = 1.5, pad = F) +
  geom_function(fun = full_cdf, color = "#2255BB", linewidth = 1.5, linetype = "dashed", xlim = c(min_af_filter, max_af_filter)) +
  xlab("Minor Allele Frequency") +
  ylab("Cumulative Probability Density") +
  # xlim(c(0,0.01)) +
  # ylim(c(0,0.01)) +
  # scale_x_continuous(trans='log') +
  # scale_y_continuous(trans='log') +
  theme_minimal()

print(p3)

fig2 <- plot_grid(p1, p2, p3, labels = "AUTO", ncol = 3)
fig2
ggsave("figs/allele-frequency.pdf", width = 12, height = 4)
ggsave("figs/allele-frequency.png", width = 12, height = 4)



