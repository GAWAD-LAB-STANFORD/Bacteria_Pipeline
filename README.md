# Guide to Bacteria_Pipeline

- [Purpose](#purpose)
- [How To Run](#how-to-run)
- [What It Does Exactly](#what-it-does-exactly)
- [Resources](#resources)
- [TODO and Notes](#todo-and-notes)

## Purpose
- This pipeline is built to identify bacterial species from pair-end fastq.gz files and remove human contamination

## How To Run
- All scripts are controlled by the master script, pipeline_control.sh
- To make job submission easier, use the submit_all.sh script to submit the master script
    - For help, use the *-h* or *--help* option like:

    ```bash
    sh submit_all.sh --help
    ```

- You can use ~ if your file or folder is in your home directory, and you can exclude a path if your file or folder is in the results directory, but otherwise use absolute paths
- You can have '/' or nothing at the end of a directory path, either is fine:
    - -r /home/groups/cgawad/results/
    - -r /home/groups/cgawad/results

### submit_all.sh
- **Run this script in order to run the entire pipeline**
- Required arguments: -p/--project and either -f/--fastq_dir or -r/--results_dir
    - Specify where your fastq.gz files are located (or will be put after you opt for the auto-demultiplexing) using *-f* or *--fastq_dir* and/or specify where to output the results using *-r* or *--results_dir*
        - If you do not specify a fastq directory, the program will assume the fastq.gz files are in the results directory you specified, and will end the program if no fastq.gz files are found
        - If you do not specify a results directory, the program will make a new folder with the current date in the name within the fastq directory 
    - Specify the project name for the final resulting VCF that will be made using *-p* or *--project*
- Optional arguments: -b/--run_dir, --sample_sheet, --R1_suffix, --R2_suffix, --err_out_dir, --contig_len_min, --slurm
    - You can have the script demultiplex your BCL files into fastq.gz files by specifying a run folder using *-b* or *--run_dir*. The program will look for a sample sheet called SampleSheet.csv in the first level within the run_dir or you can specify a different sample sheet with *--sample_sheet*. The program will make the fastq directory if it does not exist and tell you the sizes of undeteremined vs fully demultiplexed reads
    - If your read 1 and read 2 fastq.gz files differentiate themselves by some pattern other than _L001_R1_001.fastq.gz and _L001_R2_001.fastq.gz or _R1_001.fastq.gz and _R2_001.fastq.gz, use *--R1_suffix* and *--R2_suffix* options to let the pipeline know
    - You can specify a directory to output the standard error and out print statements of all jobs to using *--err_out_dir*
    - You can change the default minimum length for contigs that are kept from 5000 to a chosen number using *--contig_len_min*
    - Besides the already implemented job name and standard error and output print statements, you can specify additional slurm commands for the pipeline job following the use of the *--slurm* option. If you use this option, **make sure it is the last one you use**
        - A useful example would be setting a future time to run the job and asking for email notifications like so:

        ```
        ... --slurm --begin=now+12hours --mail-type=ALL 
        ```

- Get some example submissions by using *-h* or *--help* options like:

    ```bash
    sh submit_all.sh --help
    ```

## What It Does Exactly

### submit_all.sh
Will submit the pipeline_control.sh master script to run the entire pipeline. Why have an additional script instead of just submitting pipeline_control.sh directly? It makes it just a little easier for a novice Slurm user

### pipeline_control.sh
- **0_demultiplexer.sh** - Demultiplexing will be done if specified
- **1_process_sample.sh** - For each sample, a job will be run to do the following:
    - Trimmomatic trimming
    - Count the reads before trimming, the reads left after trimming, and the reads that were removed after trimming
    - Align all reads to human genome using BWA ALN
    - De novo assembly of contigs using SPAdes from non-human reads
    - Align non-human reads to contigs using BWA MEM
    - Collect alignment metrics
    - Count human, contig, and unaligned reads
    - Run kraken2 on non-human reads for microbes, plasmids, and phages
- Summarize metrics like human and contig alignment metrics, read targets, and kraken2 report output
    - Create json jtree files from kraken2 report outputs for later Rscript
- Only keep long contigs (default 5kb), remove low complexity regions of contigs, and consolidate them into files (320 contigs each)
- **2_blast_contigs.sh** - For each file of long contigs, a job will run to classify them using BLAST+ (blastn) against entire nucleotide database, plasmid database, and phage database
- Parse BLAST results
- Graph human contamination, kraken results, BLAST results, and other measures

## Resources
- How do I create a database for BLASTing?
    - Download NCBI database fasta files using a command like this:

    ```bash
    rsync --copy-links --recursive --times --verbose --progress rsync://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria*.genomic.fna.gz ./
    ```

    - Combine all those fasta files into an uncompressed single fasta file
    - Make a BLAST database from that combined single fasta file using a command like this:

    ```bash
    makeblastdb -in combined_bacteria_sequences.fasta -dbtype nucl -out bacteria
    ```

    - If you use NCBI's own BLAST database instead of creating your own (their databases are more condensed and take up less space), then download and update using a command like this:

    ```bash
    perl update_blastdb.pl --decompress nt
    ```

- How do I get annotations for my sequences?
    - Download NCBI database gbff annotation files using a command like this:

    ```bash
    rsync --copy-links --recursive --times --verbose --progress rsync://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria*.genomic.gbff.gz ./
    ```

    - Combine all those gbff files into an uncompressed single gbff file (or to do this faster with parallelization, do this step for each gbff file individually and combine the parsed annotations at the very end), then run a script like the one I made (parse_genbank_annotations.py) with your modifications

    ```bash
    python3 parse_genbank_annotations.py uncompressed_bacteria.gbff parsed_bacteria_annotations.tsv
    ```

## TODO and Notes