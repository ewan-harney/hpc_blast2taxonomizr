#!/bin/bash

#SBATCH --job-name=03_make_summary
#SBATCH -o 03_make_summary_o%j.out # Standard output
#SBATCH -e 03_make_summary_e%j.err # Standard error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=END

## load profile. Conda environment required to run Rscript command
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

USAGE=printf "
###################################################
## Using the 03_run_make_summary_files.sh script ##
###################################################
\n03_run_make_summary_files.sh -E [email address] \n
The script takes the standard taxonomizr output and generates two further taxonomizr output files. It then calls an R script that uses these 
new files and the '06' files from dada2 to make a final summary file combining the LCA, taxon path, sequence information and count data. This
file is saved to working_data and also sent to the user's email address (-E). The R script also creates three additional files that can be 
used as input for the popular R package phyloseq, and saves them to working_data (but these files are not sent to the user).\n
These are the files required for the script to run: \n
- blast_out/taxonomizr_taxon_path.tsv 
- working_data/06_ASV_seqs.fasta 
- working_data/06_ASV_counts.tsv \n
These are the files created by the script (with brief description): \n
- blast_out/taxonomizr_taxon_path_us.tsv (all spaces replaced with underscores in taxon path)
- blast_out/taxonomizr_summary_out.tsv (only lowest common ancestor and taxonomic level)
- working_data/ASV_taxa_seq_counts.tsv (taxa, sequence and count data): sent to user's email address (-E)
- working_data/ps_taxamat.tsv (taxon path in phyloseq format)
- working_data/ps_countmat.tsv (count data in phyloseq format)
- working_data/ps_phylogeny.rds (sequence phylogeny in phyloseq format) \n
Here is an example of how you might run the script on Bessemer:\n
sbatch b2t_scripts/03_run_make_summary_files.sh -E user@university.ac.uk \n\n\n"

## parse arguments
while getopts E: flag; do
    case "${flag}" in
        E) email=${OPTARG};;
    esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if  [ -z "${email}" ]; then
   printf "${USAGE}" >&2; exit 1
fi

## Set directories
MAIN_DIR=$PWD
OUT_DIR="blast_out"

## Replace spaces in the taxonomizr_taxon_path.tsv file with underscores:
awk 'NR>1' ${MAIN_DIR}/${OUT_DIR}/taxonomizr_taxon_path.tsv | tr ' ' '_' > ${MAIN_DIR}/${OUT_DIR}/taxonomizr_taxon_path_us.tsv

## Generate the summary file (only contains the lca for each asv, not full taxon path):
awk -v var="NA" -F '\t' '
{ \
    if (NF == 8) { \
        if ($(NF) != var) {print $1 "\t" "1_species" "\t" $(NF)} \
        else if ($(NF-1) != var) {print $1 "\t" "2_genus" "\t" $(NF-1)} \
        else if ($(NF-2) != var) {print $1 "\t" "3_family" "\t" $(NF-2)} \
        else if ($(NF-3) != var) {print $1 "\t" "4_order" "\t" $(NF-3)} \
        else if ($(NF-4) != var) {print $1 "\t" "5_class" "\t" $(NF-4)} \
        else if ($(NF-5) != var) {print $1 "\t" "6_phylum" "\t" $(NF-5)} \
        else if ($(NF-6) != var) {print $1 "\t" "7_superkingdom" "\t" $(NF-6)} \
        else {print $1 "\t" "8_unknown" "\t" "NA"} \
    } \
    else {print $1 "\t" "8_unknown" "\t" "check_taxonomizr"} \
}' ${MAIN_DIR}/${OUT_DIR}/taxonomizr_taxon_path_us.tsv > ${MAIN_DIR}/${OUT_DIR}/taxonomizr_summary_out.tsv


## Generate a mini summary, taxonomizr_stats.txt, saved to blast_out
awk '{print $2}' ${MAIN_DIR}/${OUT_DIR}/taxonomizr_summary_out.tsv | sort | uniq -c > ${MAIN_DIR}/${OUT_DIR}/taxonomizr_stats.txt

## make arg string to pass to R script
ARGS=""
if [ "$email" ]; then ARGS="$ARGS -E $email"; fi

# Run the R script to combine MEGAN, sequence and counts into a single summary file
Rscript $PWD/b2t_scripts/03_make_summary_files.R
