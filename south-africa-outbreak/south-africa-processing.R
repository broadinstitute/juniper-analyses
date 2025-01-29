## Validation on South Africa cluster

library(juniper0)
library(ape)
library(readxl)
library(lubridate)

metadata <- read_xlsx("south-africa-outbreak/raw/leger.xlsx", skip = 1)

# Read in FASTA, simplify names and dates
big_fasta <- read.FASTA("south-africa-outbreak/raw/nextclade.aligned.fasta")
big_names <- names(big_fasta)
big_dates <- as.Date(gsub(".*\\|", "", big_names))
big_names <- gsub("^.*?\\|", "", big_names)

# Read in amalgamated VCF
big_vcf <- read_xlsx("south-africa-outbreak/raw/vcf.xlsx", skip = 1)

# Add useless columns
big_vcf$ID <- "."
big_vcf$QUAL <- "."
big_vcf$FILTER <- "PASS"

# Create INFO column
big_vcf$INFO <- paste0(
  "DP=",
  big_vcf$Depth,
  ";AF=",
  big_vcf$AF,
  ";SB=",
  big_vcf$SB
)


# Get relevant columns
big_vcf <- big_vcf[, c(1:3, 27, 4:5, 28:30)]

# Reformat sample names
change <- which(grepl("KPCOVID-", big_vcf$Sample))
names <- big_vcf$Sample[change]
names <- gsub("REP", "", names)
names <- gsub("-R", "", names)
number <- gsub(".*-", "", names)
zeros <- sapply(4 - nchar(number), function(n){paste(rep(0,n), collapse = "")})
names <- paste0("KPCOVID_", zeros, number)
big_vcf$Sample[change] <- names

# Unique samples
unq_samples <- unique(big_vcf$Sample)

# Write VCF for each sample
for (s in unq_samples) {
  vcf <- big_vcf[big_vcf$Sample == s, ]
  vcf <- vcf[, -1]

  colnames(vcf)[1] <- "#CHROM"

  # Figure out name via metadata sheet
  id <- metadata$`Gisaid Accession`[match(s, metadata$Name)]

  write.table(vcf, file = paste0("south-africa-outbreak/input_data/vcf/", id, ".vcf"), quote = F, row.names = F, sep = "\t")
}

names(big_fasta) <- big_names
fasta <- big_fasta[metadata$Outbreak[match(gsub("\\|.*", "", big_names), metadata$`Gisaid Accession`)] == "CH1"]

# Fix dates
dates <- read.csv("south-africa-outbreak/raw/date.csv")
dates$V1 <- metadata$`Gisaid Accession`[match(dates$V1, metadata$Name)]
dates <- dates$V2[match(gsub("\\|.*", "", names(fasta)), dates$V1)]
dates <- as.Date("2000-01-01") + days(dates)
names(fasta) <- paste0(gsub("\\|.*", "", names(fasta)), "|", as.character(dates))

ref <- fasta[20]
fasta <- fasta[-20]

write.FASTA(fasta, "south-africa-outbreak/input_data/aligned.fasta")
write.FASTA(ref, "south-africa-outbreak/input_data/ref.fasta")
