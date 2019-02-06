# minion-assembler
Bash script to automate assembly of bacterial genomes

Script is run in bash with ./Minion_assembly.sh

All changes are made to variables within the script in order to select outputs, variables listed below.

experiment_name=Experiment_name_here
directory=~/Documents/Minion/$experiment_name
create_folder_structure=false
sort_barcodes_into_folders=false
demultiplexing=false
read_trimming=false
read_trim_length=2000
unicyler_hybrid_assembly=true
miniasm_assembly=true
nanopolish_miniasm=false
indel_checking=false
annotate=false			#annotations at the moment are only based on hybrid assemblies (trimmed or untrimmed) due to the below.
pangenome=false 		#pangenomes should be produced with hybrid assemblies as truncated proteins caused by indels in non-hybrid assemblies affect this
auto_shutdown=false
trimmed=true
threads=48
raw_minion_fastq=$directory/raw_files/fastq
short_reads_gzipped="short reads to be placed in $directory/sample/name_of_sample/reads/ with naming format ending in R1_001.fastq.gz and R2_001.fastq.gz"
