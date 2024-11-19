#!/bin/bash

#SBATCH --job-name=01B_blastn_array
#SBATCH -o blast_logs/01B_blastn_array_%a_o%j # Standard output
#SBATCH -e blast_logs/01B_blastn_array_%a_e%j # Standard error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile
source ~/.bash_profile

USAGE="
##############################################
## Using the 01B_run_blastn_array.sh script ##
##############################################
\n01B_run_blastn_array.sh -B [absolute path to blast db] -N [number of files in the array] \n
This script should only be used following script 01A_run_prep_for_blast.sh. The job is submitted as a slurm array, speeding up the identification 
of sequences.It will make a directory called blast_out and save the results here. Before running the script you must determine the size of the 
array. This information is the number appearing in the name of a file created by 01A_run_prep_for_blast.sh, and can be seen using the 'ls' command:\n
ls split_fasta/split_fasta_list_of_*.txt \n
The number (where the * appears) must be entered TWO TIMES: 1) as the array limit AND 2) as a mandatory argument for -N.\n\n
###############
## IMPORTANT ##
###############\n
Array job submission differs to normal batch job submission. Here is an example of how you might run the script on Bessemer if you have an array 
size of 24:\n
sbatch --array=1-24 b2t_scripts/01B_run_blastn_array.sh -B /shared/genomicsdb2/shared/ncbi_nt/current/nt -N 24 \n\n\n"

## List arguments
while getopts B:N: flag; do
	case "${flag}" in
		B) DATABASE=${OPTARG};;
		N) NUM=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${DATABASE}" ] || [ -z "${NUM}" ]; then
   printf "${USAGE}" >&2; exit 1
fi

## PREP STEPS
## Define path variables
MAIN_DIR=$PWD
PREP_DIR="split_fasta"
OUT_DIR="blast_out"
mkdir -p ${MAIN_DIR}/${OUT_DIR}
DATA=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat ${PREP_DIR}/${PREP_DIR}_list_of_${NUM}.txt))
FASTA=$(echo "$DATA" | cut -f1 )

## run blast 
blastn -query ${MAIN_DIR}/${PREP_DIR}/${FASTA} -task blastn -db ${DATABASE} -out ${MAIN_DIR}/${OUT_DIR}/${FASTA}_blast.out.tab -num_threads 4 -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxid ssciname scomnames sblastname sskingdoms stitle"
