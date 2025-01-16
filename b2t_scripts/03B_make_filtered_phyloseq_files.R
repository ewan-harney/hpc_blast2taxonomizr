#load msa and phangorn libraries
library(msa)
library(phangorn)
library(optparse)

# parse command line options
option_list = list(
  make_option(c("-E", "--email"), type="character", default=NULL,
              help="Provide an email address to receive an email notification when the job has finished.", metavar="character"))

# Parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# get path
path<-getwd()

## read in the file
results <- read.table(file = paste(path, "/working_data/ASV_taxa_seq_counts_filtered.tsv", sep=""), header = TRUE)
row.names(results) <- results$asv

##  Redo phyloseq file creation after filtering ASV_taxa_seq_counts.tsv
##  a. taxmat
taxmat <- results[,c(4:10)]
write.table(taxmat, "working_data/ps_taxamat_filt.tsv", sep = "\t", quote=F, col.names=T, row.names = T)

##  b. otumat (counts)
countmat <- results[,-c(1:11)]
write.table(countmat, "working_data/ps_countmat_filt.tsv", sep = "\t", quote=F, col.names=T, row.names = T)

##  c. phylogeny
##   i. setup
seqs <- results$sequence
names(seqs) <- row.names(results)
##   ii. alignment and tree creation
mult <- msa(seqs, method="ClustalW", type="dna", order="input")
phang.align <- as.phyDat(mult, type="DNA", names=seqs)
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR$tree, file = "working_data/ps_phylogeny_filt.rds")

## Email user
email_plot_command <- paste("echo \"Phyloseq results have been created following filtration of summary filesbatch.\" | mail -s \"Filtered phyloseq objects\" -a working_data/ps_taxamat_filt.tsv -a working_data/ps_countmat_filt.tsv -a working_data/ps_phylogeny_filt.rds", opt$email, sep=" ")
system(email_plot_command)
