## Process files for H5N1 run

meta <- read.csv("h5n1/raw/metadata.cattle.tsv", sep = "\t")

meta[which(grepl("SRR29694483", meta$sra_accessions)), ]

fasta_sra <- list.files("h5n1/avian-influenza/fasta/")
fasta_sra <- gsub("\\_.*", "", fasta_sra)
fasta_sra <- unique(fasta_sra)

vcf_sra <- list.files("h5n1/avian-influenza/variants/")
vcf_sra <- gsub("\\_.*", "", vcf_sra)
vcf_sra <- unique(vcf_sra)

# SRA numbers
sra <- meta$sra_accessions
sra <- strsplit(sra, ",")

# Genes
genes <- c("HA","MP","NA","NP","NS","PA","PB1","PB2")

# Figure out which one has a matching fasta
first <- T
for(i in 1:length(sra)){
  sra[[i]] <- sra[[i]][sra[[i]] %in% fasta_sra]


  # Update metadata
  if(length(sra[[i]]) > 0){
    meta$sra_accessions[i] <- sra[[i]]

    if(first){
      gene_fastas <- list()
      for (j in 1:8) {
        gene_fastas[[j]] <- read.FASTA(paste0("h5n1/avian-influenza/fasta/", sra[[i]], "_", genes[j], "_cns.fa"))
      }
      first <- F
    }else{
      for (j in 1:8) {
        gene_fastas[[j]] <- c(gene_fastas[[j]], read.FASTA(paste0("h5n1/avian-influenza/fasta/", sra[[i]], "_", genes[j], "_cns.fa")))
      }
    }

  }else{
    meta$sra_accessions[i] <- ""
  }

  print(i)
}

# Subset metadata to sampled cases
meta <- meta[meta$sra_accessions != "", ]
n <- nrow(meta)

save(meta, file = "h5n1/raw/meta_subsampled.RData")


for (j in 1:8) {
  write.FASTA(gene_fastas[[j]], paste0("h5n1/alignment/", genes[j], "_unaligned.fasta"))
}

ref <- read.FASTA("h5n1/avian-influenza/reference/reference.fasta")
for (j in 1:8) {
  write.FASTA(ref[j], paste0("h5n1/alignment/", genes[j], "_ref.fasta"))
}

# Length of each segment
lengths <- unname(sapply(ref, length))
offset <- cumsum(lengths)

## Turns out the only genes that needed alignement were PA, PB1, and PB2. We update those here
gene_fastas[[6]] <- read.FASTA("h5n1/alignment/PA_aligned.fasta")
gene_fastas[[7]] <- read.FASTA("h5n1/alignment/PB1_aligned.fasta")
gene_fastas[[8]] <- read.FASTA("h5n1/alignment/PB2_aligned.fasta")

# Concatenate into one fasta
fasta <- gene_fastas[[1]]

for (j in 2:8) {
  for (i in 1:n){
    fasta[[i]] <- c(fasta[[i]], gene_fastas[[j]][[i]])
  }
}

names(fasta) <- paste0(meta$sra_accessions, "|", meta$date)
srrs <- gsub("\\|.*", "", names(fasta))

# Delete duplicated sample
fasta <- fasta[-which(duplicated(srrs))]

write.FASTA(fasta, "h5n1/input_data/aligned.fasta")

names(fasta) <- gsub("\\|.*", "", names(fasta))



## Compile VCFs
for (i in 1:n) {
  vcf <- read.csv(
    paste0("h5n1/avian-influenza/variants/", names(fasta)[i], "_", genes[1], "_variants.tsv"),
    sep = "\t",
    colClasses = c(
      "character",
      "integer",
      "character",
      "character",
      "integer",
      "integer",
      "integer",
      "integer",
      "integer",
      "integer",
      "numeric",
      "integer",
      "numeric",
      "logical",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character"
    )
  )

  for (j in 2:8) {
    new <- read.csv(
      paste0("h5n1/avian-influenza/variants/", names(fasta)[i], "_", genes[j], "_variants.tsv"),
      sep = "\t",
      colClasses = c(
        "character",
        "integer",
        "character",
        "character",
        "integer",
        "integer",
        "integer",
        "integer",
        "integer",
        "integer",
        "numeric",
        "integer",
        "numeric",
        "logical",
        "character",
        "character",
        "character",
        "character",
        "character",
        "character"
      )
    )
    new$POS <- new$POS + offset[j-1]
    vcf <- rbind(vcf, new)
  }

  ## Match expected VCF format

  # Filter out non-substitution positions
  vcf <- vcf[vcf$ALT %in% c("A", "C", "G", "T"), ]
  vcf <- vcf[vcf$REF %in% c("A", "C", "G", "T"), ]

  # Get strand bias
  sb <- rep(0, nrow(vcf))

  for (k in 1:nrow(vcf)) {
    m <- matrix(c(vcf$REF_DP[k] - vcf$REF_RV[k], vcf$REF_RV[k], vcf$ALT_DP[k] - vcf$ALT_RV[k], vcf$ALT_RV[k]), ncol = 2)
    sb[k] <- round(-10 * log10(fisher.test(m)$p))
  }

  info <- paste0(
    "DP=",
    vcf$TOTAL_DP,
    ";AF=",
    vcf$ALT_FREQ,
    ";SB=",
    sb
  )

  vcf$ID <- "."
  vcf$QUAL <- "."
  vcf$INFO <- info

  vcf <- vcf[, c("REGION", "POS", "ID", "REF", "ALT", "QUAL", "PASS", "INFO")]
  vcf$REGION <- "PP755589.1"
  colnames(vcf)[1] <- "#CHROM"
  colnames(vcf)[7] <- "FILTER"

  ## De-duplicate rows for same iSNV
  delete <- which(duplicated(paste0(vcf$REF, vcf$POS, vcf$ALT)))
  vcf <- vcf[-delete, ]

  write.table(vcf, paste0("h5n1/input_data/vcf/", names(fasta)[i], ".vcf"), quote = F, row.names = F, sep = "\t")


  print(i)

}

length(list.files("h5n1/input_data/vcf/"))





