# Load libraries
library(taxonomizr)
library(dplyr)
library(optparse)

# parse command line options
option_list = list(
  make_option(c("-E", "--email"), type="character", default=NULL,
              help="Provide an email address to receive an email notification when the job has finished.", metavar="character"),
  make_option(c("-T", "--topperc"), type="numeric", default=NULL,
              help="Top percent value for taxonomizr", metavar="character"))

# Parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Get the directory containing the raw files
path<-getwd()

# Read in the filtered_blast.out.tab file
blastResults<-read.table(file = paste(path, "/filtered_blast.out.tab", sep=""),header=FALSE,stringsAsFactors=FALSE)

# Filter file by (user-supplied) top percent value
# Note that the top percent value is converted to a fraction and subtracted from 1:
# thus a top percent value of 2 becomes 0.98 
topPercent <- data.frame(blastResults %>% group_by(V1) %>% 
  arrange(V1, desc(V12)) %>% 
  filter(V12 > max(V12*(1 - (opt$topperc/100)))))

# Split the top percent values in to asvs and accessions
asvs<-topPercent[,1]
accessions<-topPercent[,2]

# Get the taxonomic information for top percent accessions from accessionTaxa.sql
taxaId<-accessionToTaxa(accessions,"./../accessionTaxa.sql")
taxa<-getTaxonomy(taxaId,'./../accessionTaxa.sql')

# LCA algorithm and merge back to top percent asvs
taxon_path <- condenseTaxa(taxa,asvs)

# Write output
write.table(taxon_path, "taxonomizr_taxon_path.tsv", sep="\t", quote = FALSE)

# Email user
email_plot_command <- paste("echo \"Taxonomizr is complete.\" | mail -s \"Taxonomizr table\" -a taxonomizr_taxon_path.tsv", opt$email, sep=" ")
system(email_plot_command)
