## Assigning taxonomy with BLASTn and taxonomizr on UoS BESSEMER.
This short HPC tutorial uses BLASTn to query the identity of nucleotide sequence data against an ncbi database, and applies the taxonomizr LCA (lowest common ancestor) algorithm to provide a probable taxonomic assignment.
<br></br>
<font size="4">
<details><summary><font size="6"><b>1) About, credits, and other information</b></font></summary>
<br></br>

The workflow was designed for eukaryotic taxonomic assignment. It is an alternative to the _11 Assign taxonomy to ASVs_ step in Katy Maher's [dada2 pipeline](https://github.com/khmaher/HPC_dada2), which uses curated databases for bacterial and microbial sequence data.
<br></br>

For metabarcoding projects featuring non-microbial eukaryotic data, using blastn to query sequences against ncbi databases provides a reliable means of taxonomic identification. However, the ncbi database will probably contain multiple sequences that match the query sequence, making it hard to tell which taxa the query sequence belongs to. THE LCA (lowest common ancestor) algorthm assigns taxonomy based on the lowest common taxonomic ancestor shared by the different taxonomic groups identified by BLAST. 
<br></br>

This pipeline was originally designed to work with the [MEGAN](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) lowest common ancestor (LCA) algorithm (pipeline still available at the [hpc_blast2megan](https://github.com/ewan-harney/hpc_blast2megan) repository). Unfortunately the taxonomy database for MEGAN hasn't been updated since February 2022 (this was true as of September 2024), which can lead to missing results. The R package [taxonomizr](https://github.com/sherrillmix/taxonomizr) provides a useful alternative because it incorporates a function to download the latest taxonomy data from ncbi. We run that function every month, maintaining an up-to-date version of the database on the University of Sheffield's [BESSEMER](https://docs.hpc.shef.ac.uk/en/latest/bessemer/index.html) HPC system. The blast2taxonomizr pipeline then uses the `condenseTaxa()` function of taxonomizr to provide LCA information.
<br></br>

Although designed to follow the dada2 pipeline, this workflow can be applied to query any nucleotide sequence data in fasta format. The code has been written for use with BESSEMER but should be applicable to any GNU/Linux based HPC system once the appropriate modifications are made (your mileage may vary).
<br></br>

Code which the user must run is highlighted in a code block like this:

```
I am code - you must run me
```

When we refer to programs, functions, filepaths, directory names, and file names within normal text we use inline code, like this:

* Running `blastn` will create the `all_blast.out.tab` file within the `blast_out` directory
<br></br>

When a specific button should be pressed this is shown in double quotes, like so:

* press "y" then "enter"
<br></br>

Contact: Ewan Harney //  e.harney@sheffield.ac.uk
<br>
</details>
<br>

<details><summary><font size="6"><b>2) Getting set up</b></font></summary>
<br></br>

<font size="4"><b>2.1) Access the HPC</b></font>
<br></br>
This workflow assumes you have already been using BESSEMER to run the dada2 pipeline. If that's not the case and you wish to get set up on this particular HPC, please refer to sections 2.1 - 2.4 of the [dada2 pipeline](https://github.com/khmaher/HPC_dada2)
<br></br>

<font size="4"><b>2.2) Navigate to your working directory</b></font>
<br></br>
Navigate to your project directory. If you have been running the dada2 analysis you likely have a `my_project` directory within the `/fastdata` directory on BESSEMER. Within `my_project` is the `working_data` directory, which contains the sequence data in a file called `06_ASV_seqs.fasta`. If you have not been running the dada2 pipeline, you can navigate to the directory containing your sequence data or a parent directory, whichever you prefer (the relative path to the fasta file will be specified when running the script).
<br></br>

<font size="4"><b>2.3) Copy blast2taxonomizr scripts</b></font>
<br></br>
Clone (download) this github repository, and copy the `b2t_scripts` directory contained within to your current location.
  
```
git clone "https://github.com/ewan-harney/hpc_blast2taxonomizr"
cp -r hpc_blast2taxonomizr/b2t_scripts .
```

Check the contents of the b2t_scripts directory. There should be 7 script files in the directory: 5 .sh files and 2 .R scripts:

```
ls b2t_scripts
```

<font size="4"><b>2.4) A note on editing scripts</b></font>
<br></br>
Unlike the scripts in the dada2 pipeline, or the 02 and 03 scripts in this pipeline, the 01 blast scripts (01, 01A, 01B) do not require an email address as a command line argument. By default the user will not receive email confirmation of job completion; however, this can be easily altered through a small change to the resource request section at the top of the .sh scripts. A script can be viewed and edited with the `nano` command and the relative or absolute path to the script, e.g. :

```
nano b2t_scripts/01_run_blastn_simple.sh
```

This will start the text editor `nano`. Notice that the first line of the script is `#!/bin/bash`, followed by an empty line, and then several lines commencing with `#SBATCH`. These `#SBATCH` arguments are used by slurm when the script is submitted (with `sbatch` or `qsub`) and allow the user to control certain parameters relating to job submission. Notice that the last `#SBATCH` line is `#SBATCH --mail-user=user@uni.ac.uk`. 

Using the arrow key, go to this line and change user@uni.ac.uk to your own email address. In my case I would edit the code so that it reads `#SBATCH --mail-user=e.harney@sheffield.ac.uk`.

Once you have made this change, you will need to save it. Notice at the bottom of the screen are lines of commands, such as "^G Get Help" and "^X Exit" etc. The "^" means holding down the control (windows) or command (mac) key. Pressing the "X" key whilst holding down control/command will allow you to save (if there have been changes) and exit. After pressing "^X" you will be prompted to save the changes (options are "y" for yes, "n" for no and "^C" for cancel). Press "y". You will then be given the chance to rename the file if you want. In our case we want to keep the old name, so simply press "enter" to save the file with the same name and exit `nano`.
<br>
</details>
<br>

<details><summary><font size="6"><b>3) Blast sequence data against an ncbi database</font></b></summary>
<br></br>

<font size="4"><b>3.1) Determine how many sequences are in your fasta file</b></font>
<br></br>

Depending on the size and number of sequences in our fasta file and the size of the database being used, this step can be quite slow. If your fasta file contains thousands of sequences, we can speed things up by slitting the fasta file into chunks and running `blastn` in parallel using the array functionality of slurm. This workflow therefore contains 2 different options for running blastn:

1. If you have < 1000 sequences we suggest running the single script: `01_run_blastn_simple.sh`
2. If you have > 1000 sequences we suggest splitting the file into chunks running it in array mode with the A and B scripts: `01A_run_prep_for_blast.sh` and `01B_run_blastn_array.sh`

<br>
To see how many sequences are in your fasta file, run the following:

```
grep -c '>' working_data/06_ASV_seqs.fasta
```

Although we suggest 1000 sequences as a threshold, you can run an array with less than 1000 sequences or the simple blast more than 1000 sequences (although this may take days to run). In section 3.2 we describe how to run `blastn` in simple mode, and in section 3.3 we describe how to run it in array mode.
<br></br>

<font size="4"><b>3.2) Running blastn in simple mode</b></font>
<br></br>

Running `blastn` in simple mode will create a new directory called `blast_out` in your current directory, as well as symbolic links to the ncbi taxadb files `taxdb.btd` and `taxdb.bti`. It will then run blastn and the output will be saved as `blast_out/all_blast.out.tab`. 
<br></br>

<b>To run `01_run_blastn_simple.sh` you need to provide:</b>
* the relative path to the fasta file containing the sequence data (-F)
* the location of an ncbi database on the HPC (-B)
<br></br>

It is most likely that you will use the nt database, which contains all nucleotide sequences available on GenBank. However, you can also supply a smaller or bespoke indexed database. Here we assume you are using nt. This is an example of how to submit the job on BESSEMER:
  
```
sbatch b2t_scripts/01_run_blastn_simple.sh -F working_data/06_ASV_seqs.fasta -B /shared/genomicsdb2/shared/ncbi_nt/current/nt
```

<font size="4"><b>3.3) Running blastn in array mode</b></font>
<br></br>

Slurm job arrays allow batch jobs to be broken down into parts and run in parallel, saving time for the user. Submitting these scripts is somewhat different to submitting normal sbatch jobs. For more information on arrays refer to the Sheffield HPC documentation on [advanced job submission](https://docs.hpc.shef.ac.uk/en/latest/hpc/scheduler/advanced_job_submission_and_control.html#gsc.tab=0). 
<br></br>

Running `blastn` in array mode requires running 2 scripts one after the other: `01A_run_split_fasta.sh` then `01B_run_blastn_array.sh`.
<br></br>

The `01A_run_prep_for_blast.sh` splits the input fasta file into chunks each containing 100 sequences which are written to a new directory called `split_fasta` (files will be named `chunk0.fa`,`chunk100.fa` etc). It also creates symbolic links to the ncbi taxadb files `taxdb.btd` and `taxdb.bti`, and makes a directory called `blast_logs`; all of these will be used in the following script. Finally it creates a text file called `split_fasta_list_of_X.txt` with the names of all the chunk.fa files to be used in script 01B. In your file the 'X' will be the total number of chunk.fa files and is an important parameter to set in the `01B_run_blastn_array.sh` script.
<br></br>

<b>To run `01A_run_prep_for_blast.sh` you need to provide:</b>
* the relative path to the fasta file containing the sequence data (-F)
<br></br>

An example command if you have run the dada2 pipeline might be:
  
```
sbatch b2t_scripts/01A_run_prep_for_blast.sh -F working_data/06_ASV_seqs.fasta
```
  
The `01B_run_blastn_array.sh` script will then use an array to simultaneously blast multiple chunk.fa files against an ncbi database. This script will create a new directory called `blast_out` in your current directory and write the output of blasting each chunk against the database to a seperate output in that folder (files will be named `chunk0.fa_blast.out.tab`, `chunk100.fa_blast.out.tab` etc).
<br></br>

<b>To run `01B_run_blastn_array.sh` you need to provide:</b>
* the location of an ncbi database on the HPC (-B)
* the number of input files to be run on the array (-N)
<br></br>

As stated in section 3.2, it is most likely that you will use the ncbi database nt. The number -N is contained in the file name of `split_fasta_list_of_X.txt` (in place of the 'X'). This can be viewed with the following command:
  
```
ls split_fasta/split_fasta*
```
  
<b>IMPORTANT!</b> Check carefully the number that appears in this file name! This is the number that <b>MUST</b> be used when running `01B_run_blastn_array.sh`. For example, if our original sequence.fasta file had contained 2350 sequences, it would have been split into 24 chunks, with the txt file named `split_fasta_list_of_24.txt`. The number (in this example <b>24</b>) should appear twice when we submit this job; as the array limit (following `sbatch --array=1-` ) and as the value of `-N`. For example:
  
```
sbatch --array=1-24 b2t_scripts/01B_run_blastn_array.sh -B /shared/genomicsdb2/shared/ncbi_nt/current/nt -N 24
```
  
<b>On the other hand</b> if our original sequence.fasta file had only contained 1780 sequences, it would have been split into 18 chunks, with the txt file named `split_fasta_list_of_18.txt`. In this case our `array` and `-N` value would be set to <b>18</b> and we would run the script like so:
  
```
sbatch --array=1-18 b2t_scripts/01B_run_blastn_array.sh -B /shared/genomicsdb2/shared/ncbi_nt/current/nt -N 18
```
  
Note that error and output log files for each subjob of the array will be written to the directory `blast_logs`.
<br></br>

<font size="4"><b>3.4) Monitoring and assessing the result of blastn</b></font>
<br></br>

Running `blastn` against the nt database can take a while. To follow the status of the job run the following command: 

```
squeue --me
```
  
For more information about the `squeue` output, refer to the Sheffield HPC documentation on [squeue](https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/Common-commands/squeue.html#gsc.tab=0). `squeue` will show the status of the job, and in the case of an array, how many of the subjobs have been submitted and how many are still queued.
<br></br>  
If `blastn` was run in simple mode, `blast_out` should contain a single file called `all_blast.out.tab`, and if it was run in array mode, it will contain several chunk.fa_blast.out.tab files. Look at the contents of (one of) these file(s) with:
  
```
head blast_out/all_blast.out.tab 
```
  
or 

```
head blast_out/chunk0.fa_blast.out.tab
```
  
Information about blast tabular output can be found at the [Metagenomics wiki](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6). The column headers in your file(s) correspond to:

qseqid / *saccver* / pident / length / mismatch / gapopen / qstart / qend / sstart / send / evalue / bitscore / *staxid* / *ssciname* / *scomnames* / *sblastname* / *sskingdoms* / *stitle*

(italics highlight differences to the default blast tabular output).
<br></br> 

All the rows displayed by head (the top 10) probably show results for the same sequence (ASV_1 if following the dada2 pipeline) because usually queries match many sequences in the nt database. Sometimes the alignment will be much better for one species than any other, but in other cases the alignments from many sequences will be comparable, and we may need to class the sequence at a higher taxonomic level (e.g. genus or family). To provide a likely taxonomic assignmnet for each ASV we will apply a lowest common ancestor (LCA) algorithm.
<br>
</details>
<br>

<details><summary><font size="6"><b>4) Run the taxonomizr LCA algorithm</font></b></summary>
<br></br>

For this step we will use [taxonomizr](https://github.com/sherrillmix/taxonomizr), an R package designed to parse NCBI taxonomy files and help with taxonomy assignment. The `condenceTaxa()` function calculates the lowest common ancestor or [LCA](https://en.wikipedia.org/wiki/Lowest_common_ancestor) using multiple blast results.
<br></br>

The `02_run_taxonomizr_lca.sh` script takes the output from `blast_out`, and does the following:
<br></br>

1. If blast was run in array mode, chunks are merged (this is automatically detected)
2. Blast results are filtered by pident (percentage of identical positions: column 3), which is provided by the user
3. The sensitivity of the LCA algorithm is further adjusted be varying the top percent arguement (see below).
4. The output is a taxon path file (all the taxonomic levels to the lowest common ancestor), saved to `blast_out` and emailed to the user.
<br></br>

This script carries out some preparation on the command line and then calls an R script called `02_taxonomizr_lca.R`. It produces an intermediate file `filtered_blast.out.tab` and a final output file called `taxonomizr_taxon_path.tsv`.
<br></br>

<b>To run the `02_run_taxonomizr_lca.sh` script you must provide:</b>
* Minimum percentage identity (0-100) for the blast results to be considered by taxonomizr (-P)
* The Top Percent parameter (1-10) for lca calculation (-T)
* Full path to the taxonomizr accessionTaxa.sql database (-B)
* An email address (-E)
<br></br>

<b>Additionally, there is an optional argument to exclude hits based on a search term:</b>
* Search term enclosed within single quotes, e.g. 'uncultured' (-G)
<br></br>

Taxonomic assignment may be improved by removing uncultured sequences (which rarely feature taxonomic information beyond domain). Note that this will only work if blast was run with the files 'taxdb.btd' and 'taxdb.bti' present, which include additional taxonomic information.
<br></br>

Blast will potentially output hundred of hits for each ASV. The minimum percentage identity threshold (-P) can be used to reduce the number of hits that are considered by taxonomizr. We suggest using a relatively high value to start with (90 or 95), which can be reduced if high numbers of NAs appear in the summary file.
<br></br>

We also provide a 'Top Percent' value (-T) to further adjust sensitivity. This filter removes any blast hits where the [bit score](https://www.metagenomics.wiki/tools/blast/evalue) is less than *100 minus T%* of the highest scoring blast hit. 
<br></br>

Top Percent can be set to anything from 0 to 100, but values of 1-10 are most appropriate. Setting Top Percent to a low value (1.5, 2) will retain fewer blast hits and probably result in better taxonomic assignment. However if the ASVs derive from organisms with poor representation in the reference database then it may be better to set Top Percent higher (5-10). The idea for the Top Percent filter is taken from MEGAN: for more information see the [MEGAN manual](https://software-ab.cs.uni-tuebingen.de/download/megan6/manual.pdf).
<br></br>

If running the analysis on BESSEMER, an up-to-date version of the taxonomizr database should be available at `/shared/genomicsdb2/shared/r_taxonomizr/current/accessionTaxa.sql`. You can also generate your own version using R (note that the file is rather large). For more details please refer to [taxonomizr](https://github.com/sherrillmix/taxonomizr). For deciding which values of -P and -T to use, we recommend initially using relatively strict (high) values for both, such as:
* -P 95
* -T 2
<br></br>

Thus you might run the job like so:

```
sbatch b2t_scripts/02_run_taxonomizr_lca.sh -P 95 -T 2 -B /shared/genomicsdb2/shared/r_taxonomizr/current/accessionTaxa.sql -E user@university.ac.uk
```
If you wish to exclude hits containing 'uncultured in their name you would include the optional -G argument:
```
sbatch b2t_scripts/02_run_taxonomizr_lca.sh -P 95 -T 2 -G 'uncultured' -B /shared/genomicsdb2/shared/r_taxonomizr/current/accessionTaxa.sql -E user@university.ac.uk
```
</details>
<br>

<details><summary><font size="6"><b>5) Check results and create summary files</font></b></summary>
<br><br>
  
<font size="4"><b>5.1) Run the `03_run_make_summary_files.sh` script </b></font>
<br><br>

The `taxonomizr_taxon_path.tsv` file produced in the previous step provides us with the taxon path to the LCA of each ASV. However, we may want to create addtional summary files, as well as files formatted for downstream analysis with the community analysis R package [phyloseq](https://joey711.github.io/phyloseq/). The final script of this pipeline, `03_run_make_summary_files.sh`, creates these files. As well as using the `taxonomizr_taxon_path.tsv` file, this step requires the `06_ASV_seqs.fasta` and `06_ASV_counts.tsv` files from the [dada2 pipeline](https://github.com/khmaher/HPC_dada2), which it expects to find in `working_data`.

The only arguement that needs to be provided to run the script is an email address, where the final summary file will be sent. The script can be run as follows:
  
```
sbatch b2t_scripts/03_run_make_summary_files.sh -E user@university.ac.uk
```
<br><br>

The script will generate some intermediate files which will be saved in `blast_out` (`taxonomizr_taxon_path_us.tsv`, `taxonomizr_summary_out.tsv`, and `taxonomizr_stats.txt`). It then calls an R script, `03_make_summary_files.R` which will combine taxonomizr results with dada2 results to produce a complete summary of the results from both pipelines (`ASV_taxa_seq_counts.tsv`) and files which can be used as input for downstream phyloseq analysis (`ps_taxamat.tsv`, `ps_countmat.tsv`, and `ps_phylogeny.rds`) all of which will be saved in `working_data`. When the script finishes running, the main summary (`ASV_taxa_seq_counts.tsv`) will be emailed to the user.
<br><br>

<font size="4"><b>5.2) Understanding the various output files </b></font>
<br><br>
Here we provide some information about the output generated by `03_run_make_summary_files.sh`. Initially three files are created in `blast_out`: 

* `taxonomizr_taxon_path_us.tsv` : similar to the original taxonomizr_taxon_path.tsv file, but with all spaces replaced by underscores (which makes it easier to parse).
* `taxonomizr_summary_out.tsv` : contains only the LCA and the LCA rank for each ASV (not the full taxon path)
* `taxonomizr_stats.txt` : provides an indication of taxonomizr performance.
<br></br>

These files are generated within a few seconds of the script beginning, so they can be inspected while you are waiting for the script to finish and send you the complete summary file. In particular, `taxonomizr_stats.txt`, is useful to view, as it shows how many ASV sequences have been assigned to each taxonomic rank. Ideally we would like most sequences to be are resolved to the species or genus level. You can view the file by running the following:
<br><br>
```
head blast_out/taxonomizr_stats.txt
```
It is common to have some ASVs with poor taxonomic assignment in these kinds of data sets, whether that is taxa only being assigned to higher levels (e.g. phylum or superkingdom), or not being assinged at all (here classed as NA). Deciding what is a 'good' result will depend on many factors including experiment type, sampling strategy, primers used etc. If you are not satisfied with your results, you can try tweaking the values of -P and -T in `02_run_taxonomizr_lca.sh`.
<br><br>

Reducing -P allows sequences with lower blast percentage identity to be considered by taxonomizr, and may reduce the number of sequences without assignment.  
Reducing -T will make the taxonomizr lca algorithm more stringent, potentially providing more specific taxonomic assignement.
<br><br>

The remaining files are created in working_data.
* `ASV_taxa_seq_counts.tsv` : complete summary of taxonomic results (lca taxon and taxon path to the lca taxon), sequence, and count results
* `ps_taxamat.tsv` : ASV taxon paths in phyloseq format
* `ps_countmat.tsv` : ASV count data in phyloseq format
* `ps_phylogeny.rds` : phylogenetic tree prepared according to protocol of [Callahan et al. 2016](https://f1000research.com/articles/5-1492/v1), see subsection _Construct the phylogenetic tree_,

Once the script finishes running the user will receive an email with the `ASV_taxa_seq_counts.tsv` file attached. Phyloseq phylogeny construction can take a while (depending on how many sequences are included in the analysis, up to 30 minutes), although the `ASV_taxa_seq_counts.tsv` file can usually be viewed within a couple of minutes of the script running.
</font>
<br>
</details>
<br>
