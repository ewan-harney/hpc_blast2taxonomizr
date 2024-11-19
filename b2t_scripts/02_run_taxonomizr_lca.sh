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

## load profile and environment
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/taxonomizr

function usage {
    echo "
###############################################
## Using the 02_run_taxonomizr_lca.sh script ##
###############################################

Usage: $0 -P [blast percent ID] -T [top percent] -G ['search term' to exclude] -L [minimumn alignment length] -B [absolute path to accessionTaxa.sql] -O [input/out pirectory for blast results] -E [email address]
"
    echo "This script merges results if blast was run in array mode (automatically detected), then applies various filters (-P is mandatory; -G and -L are optional), before running the taxonomizr algorithm (a value of -T must be provided). The absolute path to the database (-B), and input/output directory (-O) can be provided as optional arguments (see below for defaults). Results are emailed to the user (-E).
    "
    echo "  -P Filter blast results by minimum percentage identity from 0-100. We recommend 97."
    echo "  -G (optional) Remove blast hits containing a search term e.g. 'uncultured', or 'uncultured\|environmental'."
    echo "  -L (optional) Filter results by a minimum alignment length (expressed as a percentage of ASV length) from 0-100. We recommend 75."
    echo "  -T Top percent value for the taxonomizr LCA algorithm (1-10). We recommend 2 (decimals permitted)."
    echo "  -B (optional) Absolute path to accessionTaxa.sql database. Defaults to '/shared/genomicsdb2/shared/r_taxonomizr/current/accessionTaxa.sql'."
    echo "  -O (optional) Relative path to blast results. Intermediate and final output will also be written here. Defaults to 'blast_out'."
    echo "  -E User email address.
    "
    echo "Here is an example of how you might run the script on Bessemer for a 200 bp expected amplicon, assuming default database and in/out paths:
    "
    echo "sbatch b2t_scripts/02_run_taxonomizr_lca.sh -P 97 -G 'uncultured' -L 150 -T 2 -E user@university.ac.uk
    "
    exit 1
}

## parse arguments
while getopts E:B:O:P:T:G:L: flag
do
    case "${flag}" in
        E) email=${OPTARG};;
        B) database=${OPTARG};;
        O) customdir=${OPTARG};;
        P) BPI=${OPTARG};;
        T) topperc=${OPTARG};;
        G) grep=${OPTARG};;
        L) minlen=${OPTARG};;
    esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${email}" ] || [ -z "${BPI}" ] || [ -z "${topperc}" ]; 
then
    echo "Some or all of the mandatory parameters are empty";
    usage
fi

## Echo filter parameters
echo "TAXONOMIZR PARAMETERS SET:
"
if [ "$grep" ];
then
    echo "  -G: excluding results containing the string: " ${grep}
fi

if [ "$minlen" ];
then
    echo "  -L: blast minimum alignment length percentage: " ${minlen}
fi

echo "  -P: blast percent identity value: " ${BPI}
echo "  -T: taxonomizr lca top percent value: "${topperc}

## Define path variables
## Location of project directory is MAIN_DIR
MAIN_DIR=$PWD
## Is blastpath custom or default?
if [ "$customdir" ];
then
    BLASTPATH=${customdir}
    echo "
Using custom input/output directory for blast results:" ${BLASTPATH}
else
    BLASTPATH="blast_out"
    echo "
Using default input/output directory for blast results:" ${BLASTPATH}
fi

## Check if symlink to accessionTaxa.sql exists. If not, make it, and  state whether path is default or custom
MYLINK=accessionTaxa.sql
if [ -L ${MYLINK} ] && [ -e ${MYLINK} ]; then
    echo "Symbolic link to accessionTaxa.sql already exists."
else
    echo "Symbolic link to accessionTaxa.sql not found, creating new link."
    if [ "$database" ];
    then
        DB=${database}
        echo "User has specified custom location of accessionTaxa.sql:" ${DB}
        ln -s ${DB} ${MYLINK}
    else
        DB="/shared/genomicsdb2/shared/r_taxonomizr/current/accessionTaxa.sql"
        echo "Using default location of accessionTaxa.sql:" ${DB}
        ln -s ${DB} ${MYLINK}
    fi
fi

## Exit if symlink is broken
if [ ! -e ${MYLINK} ]; then
   echo "Symbolic link not working, please check path" >&2; exit 1
else
   echo "Symbolic link is working"
fi

## Change to BLASTPATH directory
cd ${MAIN_DIR}/${BLASTPATH}

## Step 1: Check if run in array mode and, if so, merge the chunks to create all_blast.out.tab
if [ -f "chunk0.fa_blast.out.tab" ]; 
then
    echo "
Blast was run in array mode, merging chunks...
"
    cat chunk* > all_blast.out.tab
fi

## Apply optional filters (-G and -L) if specified. 
## Always filter by blast percent identity (-P), and use cut to remove additional taxa information before taxonomizr 
if [ "$grep" ] && [ ! "$minlen" ] ;
then
    grep -v "${grep}" all_blast.out.tab | cut -f1-12 | awk -v varPI="${BPI}" '$3 >= varPI' > filtered_blast.out.tab
elif [ ! "$grep" ] && [ "$minlen" ] ;
then
    awk -F "\t" '{print $0 "\t" $4/$13*100}' all_blast.out.tab | awk -F "\t" -v varML="${minlen}" '$20 >= varML' | cut -f1-12 | awk -v varPI="${BPI}" '$3 >= varPI' > filtered_blast.out.tab
elif [ "$grep" ] && [ "$minlen" ] ;
then
    grep -v "${grep}" all_blast.out.tab | awk -F "\t" '{print $0 "\t" $4/$13*100}' | awk -F "\t" -v varML="${minlen}" '$20 >= varML' | cut -f1-12 | awk -v varPI="${BPI}" '$3 >= varPI' > filtered_blast.out.tab
else
    cut -f1-12 all_blast.out.tab | awk -v varPI="${BPI}" '$3 >= varPI' > filtered_blast.out.tab
fi

## build up arg string of email and topperc to pass to R script
ARGS=""
if [ "$email" ]; then ARGS="$ARGS -E $email"; fi
if [ "$topperc" ]; then ARGS="$ARGS -T $topperc"; fi

echo "Filtering of blast results complete, proceeding to taxonomizr_lca.R script 
"

## Step 3: Run taxonomizr lca by calling Rscript
Rscript ${MAIN_DIR}/b2t_scripts/02_taxonomizr_lca.R $ARGS
