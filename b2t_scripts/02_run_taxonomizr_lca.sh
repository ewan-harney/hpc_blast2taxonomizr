#!/bin/bash

#SBATCH --job-name=02_taxonomizr_lca
#SBATCH -o 02_taxonomizr_lca_o%j # Standard output
#SBATCH -e 02_taxonomizr_lca_e%j # Standard error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem=4GB
#SBATCH --time=24:00:00

USAGE="
###############################################
## Using the 02_run_taxonomizr_lca.sh script ##
###############################################
\n02_run_taxonomizr_lca.sh -P [blast percent ID value 0-100] -T [top percent value 1-10] -B [full path to accessionTaxa.sql file] -E [email address] \n
The script takes the output from blast (script 01 or script 01A and 01B) which are assumed to be in the directory blast_out, and applies the 
following steps: \n
Step 1: Merge results if blast was run in array mode (automatically detected).
Step 2: Filter blast results by percentage ID (0-100). The user provides a minimum percentage ID (-P) for 
        filtering blast results. We recommend 85-95.
Step 3: The user can further adjust the sensitivity of the LCA algorithm by providing a top percent value
        (-T) between 1 and 10 (decimals permitted). 
Step 4: Run the taxonomizr lca (lowest common ancestor) algorithm; this requires the user to provide the 
        full path to the accessionTaxa.sql database (-B). The output of taxonomizr will be saved to 
        blast_out and also sent to the user's email address (-E).\n
Here is an example of how you might run the script on Bessemer:\n
sbatch b2t_scripts/02_run_taxonomizr_lca.sh -P 90 -T 2 -B /shared/genomicsdb2/shared/r_taxonomizr/current/accessionTaxa.sql -E user@university.ac.uk \n\n\n"

## load profile and environment
source ~/.bash_profile
mamba activate taxonomizr

## parse arguments
while getopts E:B:P:T: flag
do
	case "${flag}" in
		E) email=${OPTARG};;
		B) database=${OPTARG};;
		P) BPI=${OPTARG};;
		T) topperc=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if  [ -z "${email}" ] || [ -z "${database}" ] || [ -z "${BPI}" ] || [ -z "${topperc}" ]; then
   printf "${USAGE}" >&2; exit 1
fi

echo "blast percent identity value = " ${BPI}
echo "taxonomizr lca top percent value = "${topperc}

## Define path variables
MAIN_DIR=$PWD
OUT_DIR="blast_out"

## Step 1 : Check if run in array mode and, if so, merge the chunks to create all_blast.out.tab
if [ -f "${MAIN_DIR}/${OUT_DIR}/chunk0.fa_blast.out.tab" ]; then
  echo "Blast was run in array mode, merging chunks..."
  cat ${MAIN_DIR}/${OUT_DIR}/chunk* > ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab
fi

## Step 2: Remove additional taxa information and filter by user specified blast percentage identity (BPI)
cut -f1-12 ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab | awk -v var="${BPI}" '$3 >= var' > ${MAIN_DIR}/${OUT_DIR}/filtered_blast.out.tab

## Make a symlink for the most recent version of the accessionTaxa.sql database for taxonomizr to use:
ln -s ${database} accessionTaxa.sql

## build up arg string of email and topperc to pass to R script
ARGS=""
if [ "$email" ]; then ARGS="$ARGS -E $email"; fi
if [ "$topperc" ]; then ARGS="$ARGS -T $topperc"; fi

## Step 3: Run taxonomizr lca by calling Rscript
Rscript $PWD/b2t_scripts/02_taxonomizr_lca.R $ARGS
