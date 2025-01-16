#!/bin/bash

#SBATCH --job-name=03_make_summary
#SBATCH -o 03_make_filtpseq_o%j.out # Standard output
#SBATCH -e 03_make_filtpseq_e%j.err # Standard error
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
## Using the 03_run_make_filtered_phyloseq_files.sh script ##
###################################################
\n03b_run_make_filtered_phyloseq_files.sh -E [email address] \n
The script is an optional extra step for creating phyloseq files following additonal filtration steps of the ASVs. The user should apply their own filters to the main 
summary output file from script 03_run_make_summary_files.sh (a file called ASV_taxa_seq_counts.tsv). Filtration may be a necessity if you have many uninformative ASVs, 
and you wish to create a phylogeny for downstream analyses with the R package phyloseq. Filters could include removing any ASVs in which taxonomic assignment is for 
non-target organisms (e.g. removing prokaryotes if the targets were eukaryotes), or if taxonomic assignment for some ASVs is poor (for example, if assignment is only 
to a higher level, such as phylum). Once filtered, the file should be saved in the same .tsv format with the following file name:\n
ASV_taxa_seq_counts_filtered.tsv\n
These are the files created by the script (with brief description): \n
- working_data/ps_taxamat_filt.tsv (taxon path in phyloseq format)
- working_data/ps_countmat_filt.tsv (count data in phyloseq format)
- working_data/ps_phylogeny_filt.rds (sequence phylogeny in phyloseq format) \n
Here is an example of how you might run the script on Bessemer:\n
sbatch b2t_scripts/03b_run_make_filtered_phyloseq_files.sh -E user@university.ac.uk \n\n\n"

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

## Set directory
MAIN_DIR=$PWD

## make arg string to pass to R script
ARGS=""
if [ "$email" ]; then ARGS="$ARGS -E $email"; fi

# Run the R script to combine MEGAN, sequence and counts into a single summary file
Rscript $PWD/b2t_scripts/03b_make_filtered_phyloseq_files.R $ARGS
