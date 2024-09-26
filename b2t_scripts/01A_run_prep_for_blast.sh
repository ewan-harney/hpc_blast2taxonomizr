#!/bin/bash

#SBATCH --job-name=01A_blast_prep
#SBATCH -o 01A_blast_prep_o%j # Standard output
#SBATCH -e 01A_blast_prep_e%j # Standard error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem=4GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile
source ~/.bash_profile

USAGE="
################################################
## Using the 01A_run_prep_for_blast.sh script ##
################################################
\n01A_run_prep_for_blast.sh -F [relative path to fasta file] \n
The script assumes you have a relatively large number of ASVs to identify (more that 1000) and splits them into chunks of 100 ASV sequences.
It makes a directory called split_fasta and saves fasta chunks here. It also creates a list of the files (required for step 01B), as well as
symbolic links to ncbi taxa databases (required for script 01B), and a blast_logs directory (where output from script 01B will be saved).
Here is an example of how you might run the script on Bessemer:\n
sbatch b2t_scripts/01A_run_prep_for_blast.sh -F working_data/06_ASV_seqs.fasta\n\n\n"

## List arguments
while getopts F: flag; do
	case "${flag}" in
		F) FASTA_PATH=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${FASTA_PATH}" ]; then
   printf "${USAGE}" >&2; exit 1
fi

## PREP STEPS
## Define path variables
MAIN_DIR=$PWD
PREP_DIR="split_fasta"
## create the directory to put split files, and a new directory to store the logs produced by the 01B script
mkdir -p ${MAIN_DIR}/${PREP_DIR}
mkdir -p ${MAIN_DIR}/blast_logs
## make a symlink to ncbi taxdb files, which need to be in the current directory for the 01B script
ln -s /shared/genomicsdb2/shared/ncbi_nt/current/taxdb* .
## Change to the split_fasta directory
cd ${MAIN_DIR}/${PREP_DIR}

## split fasta file into chunks each with 100 sequences and generate list.txt file
awk 'BEGIN {n=0;} /^>/ {if(n%100==0){file=sprintf("chunk%d.fa",n);} print > file; n++; next;} { print >> file; }' < ${MAIN_DIR}/${FASTA_PATH}
NUM=$(ls chunk* | wc -l )
ls chunk* > ${MAIN_DIR}/${PREP_DIR}/${PREP_DIR}_list_of_${NUM}.txt
