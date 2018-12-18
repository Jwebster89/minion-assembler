#!/bin/bash

############################################################################################################################################################################################################
###   	 ~BEFORE RUNNING~                                                                                                                                                                                ###
### Before running this script, please create a folder in $directory with $experiment_name with subdirectories '/raw_files/fastq' E.g. '~/Documents/Minion/My_Experiment/raw_files/fastq'. Then place    ###
### all generated fastq files in this folder. The same can be done for a fast5 folder in '/My_Experiment/raw_files/fast5'                                                                                ###
###                                                                                                                                                                                                      ###
############################################################################################################################################################################################################                                                                                                                                                                                                      

########### 	~Experiment details~     ########

experiment_name=281118_Pantoea_Astars
directory=~/Documents/Minion/$experiment_name

#################################################


### To create a folder structure for the current experiment with sample names. Provide a file in $directory with samples in barcoded order on separate lines called 'sample_names.txt'                   
### Sorting barcodes into folders will move all barcodes/trimmed barcodes into the created folder structure below. This should be done initially after demultiplexing and trimming, but should be turned 
### off when re-assembling.                                                                                                                                                                              

create_folder_structure=false
sort_barcodes_into_folders=false

#### To perform any of the following programs/assemblies, set variable to TRUE. Otherwise, set to FALSE.

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

### if reads have been previously trimmed and you are re-running an assembly on trimmed reads, set the following to 'true', otherwise leave 'false'.

trimmed=true


threads=48

#### Input files required
raw_minion_fastq=$directory/raw_files/fastq
short_reads_gzipped="short reads to be placed in $directory/sample/name_of_sample/reads/ with naming format ending in R1_001.fastq.gz and R2_001.fastq.gz"






###################################################################################################################
###### The following should not be changed unless you know what you're doing or you want to break this script #####
###################################################################################################################
if $unicyler_hybrid_assembly ; then
	echo "Have short reads been placed in the $directory/_sample_/reads/ for each sample? (y/n): "
	echo "reads should have a naming format ending in R1_001.fastq.gz and R2_001.fastq.gz"
	read var
fi

if $create_folder_structure; then
	mkdir $directory/samples
	for sample in `cat $directory/sample_names.txt`; do
		mkdir $directory/samples/$sample/
		mkdir $directory/samples/$sample/reads
	done
fi

if $demultiplexing; then
	cat $raw_minion_fastq/*fastq > $raw_minion_fastq/$experiment_name.cat.fastq
	porechop -i $raw_minion_fastq/$experiment_name.cat.fastq --discard_middle -t $threads -b $directory/barcodes
	rm $raw_minion_fastq/$experiment_name.cat.fastq
fi

if $read_trimming ; then
for sample in $directory/barcodes/BC*; do
	f=$(basename $sample)
	~/bin/Filtlong/bin/filtlong --min_length $read_trim_length --keep_percent 90 $directory/barcodes/$f > $directory/barcodes/trimmed.$f
	trimmed=true
done
fi

if $sort_barcodes_into_folders; then
	if $trimmed; then
		COUNTER=1
		for item in $directory/barcodes/trimmed.BC*; do
			f=$(basename $item)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			cp $directory/barcodes/$f $directory/samples/$location/reads/
			COUNTER=$[$COUNTER +1]
		done
	else
		COUNTER=1
		for item in $directory/barcodes/BC*; do
			f=$(basename $item)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			cp $directory/barcodes/$f $directory/samples/$location/reads/
			COUNTER=$[$COUNTER +1]
		done
	fi
fi

if $unicyler_hybrid_assembly ; then
	if [ $var == "y" ]; then
		if $trimmed; then
			COUNTER=01
			for trimmed_sample in $directory/barcodes/trimmed.BC*; do
				# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
				f=$(basename $trimmed_sample)
				location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
				unicycler -1 $directory/samples/$location/reads/*R1_001.fastq.gz -2 $directory/samples/$location/reads/*R2_001.fastq.gz -t $threads  \
				-l $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq -o $directory/samples/$location/trimmed_unicycler_hybrid
				cp $directory/samples/$location/trimmed_unicycler_hybrid/assembly.fasta $directory/samples/$location/trimmed_unicycler_hybrid/${location}_hybrid_trimmed.fasta
				assembly-stats $directory/samples/$location/trimmed_unicycler_hybrid/${location}_hybrid_trimmed.fasta > $directory/samples/$location/trimmed_unicycler_hybrid/assembly_stats_${location}.txt
				
				COUNTER=$[$COUNTER +1]
			done
		else
			COUNTER=01
			for sample in $directory/barcodes/BC*; do
				# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
				f=$(basename $sample)
				location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
				unicycler -1 $directory/samples/$location/reads/*R1_001.fastq.gz -2 $directory/samples/$location/reads/*R2_001.fastq.gz -t $threads  \
				-l $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq -o $directory/samples/$location/unicycler_hybrid
				cp $directory/samples/$location/trimmed_unicycler_hybrid/assembly.fasta $directory/samples/$location/trimmed_unicycler_hybrid/${location}_hybrid.fasta
				assembly-stats $directory/samples/$location/trimmed_unicycler_hybrid/${location}_hybrid_trimmed.fasta > $directory/samples/$location/trimmed_unicycler_hybrid/assembly_stats_${location}.txt

				COUNTER=$[$COUNTER +1]
			done
		fi
	fi
fi

if $miniasm_assembly; then
	if $trimmed; then
		COUNTER=01
		for trimmed_sample in $directory/barcodes/trimmed.BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
			f=$(basename $trimmed_sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			mkdir $directory/samples/$location/miniasm_longread/

			~/bin/minimap2/minimap2 -t $threads -x ava-ont $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq | gzip -1 > $directory/samples/$location/miniasm_longread/minimap.gz

			~/bin/miniasm/miniasm -f $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.gz > $directory/samples/$location/miniasm_longread/miniasm.gfa

			awk '/^S/{print ">"$2"\n"$3}' $directory/samples/$location/miniasm_longread/miniasm.gfa > $directory/samples/$location/miniasm_longread/miniasm.fasta
			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/miniasm.fasta $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap.racon.paf

			racon -t $threads $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.racon.paf $directory/samples/$location/miniasm_longread/miniasm.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly.fasta

			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/consensus_assembly.fasta $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap_2.racon.paf

			racon -t $threads $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap_2.racon.paf $directory/samples/$location/miniasm_longread/consensus_assembly.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta

			assembly-stats $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta > assembly_stats_${location}.txt
			cp $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta $directory/samples/$location/miniasm_longread/${location}_miniasm_trimmed.fasta
		
			COUNTER=$[$COUNTER +1]
		done
	else
		COUNTER=01
		for sample in $directory/barcodes/BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
			f=$(basename $sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			mkdir $directory/samples/$location/miniasm_longread/

			~/bin/minimap2/minimap2 -t $threads -x ava-ont $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq | gzip -1 > $directory/samples/$location/miniasm_longread/minimap.gz

			~/bin/miniasm/miniasm -f $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.gz > $directory/samples/$location/miniasm_longread/miniasm.gfa

			awk '/^S/{print ">"$2"\n"$3}' $directory/samples/$location/miniasm_longread/miniasm.gfa > $directory/samples/$location/miniasm_longread/miniasm.fasta
			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/miniasm.fasta $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap.racon.paf

			racon -t $threads $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.racon.paf $directory/samples/$location/miniasm_longread/miniasm.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly.fasta

			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/consensus_assembly.fasta $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap_2.racon.paf

			racon -t $threads $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap_2.racon.paf $directory/samples/$location/miniasm_longread/consensus_assembly.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta
		
			assembly-stats $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta > assembly_stats_${location}.txt
			cp $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta $directory/samples/$location/miniasm_longread/${location}_miniasm.fasta

			COUNTER=$[$COUNTER +1]
		done
	fi
fi
############
if $nanopolish_minasm; then
	#cat $raw_minion_fastq/*fastq > $raw_minion_fastq/$experiment_name.cat.fastq
	if $trimmed; then
		COUNTER=1
		for trimmed_sample in $directory/barcodes/trimmed.BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
			f=$(basename $trimmed_sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			#nanopolish index -d $directory/raw_files/fast5/pass $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq			
			#bwa index $directory/samples/$location/miniasm_longread/${location}_miniasm_trimmed.fasta
			mkdir $directory/samples/$location/nanopolish/
			#bwa mem -x ont2d -t $threads $directory/samples/$location/miniasm_longread/${location}_miniasm_trimmed.fasta $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq | samtools sort -o $directory/samples/$location/nanopolish/reads.sorted.${location}.BAM
			#samtools index $directory/samples/$location/nanopolish/reads.sorted.${location}.BAM
			nanopolish variants --consensus $directory/samples/$location/nanopolish/nano.${location}.fasta -r $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq -b $directory/samples/$location/nanopolish/reads.sorted.${location}.BAM -g $directory/samples/$location/miniasm_longread/${location}_miniasm_trimmed.fasta -q dcm,dam --min-candidate-frequency 0.1
			COUNTER=$[$COUNTER +1]
		done
	fi
fi

if $indel_checking; then
	for item in $directory/barcodes/BC*; do
		location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
		mv -f ~/bin/ideel/genomes/* ~/bin/ideel/completed_genomes/
		cp $directory/samples/$location/miniasm_longread/${location}_miniasm.fasta ~/bin/ideel/genomes/${location}_miniasm.fa
		cp $directory/samples/$location/trimmed_unicycler_hybrid/${location}_hybrid_trimmed.fasta ~/bin/ideel/genomes/${location}_hybrid_trimmed.fa
		cp $directory/samples/$location/trimmed_unicycler_hybrid/${location}_hybrid.fasta ~/bin/ideel/genomes/${location}_hybrid.fa
		COUNTER=$[$COUNTER +1]
	done
	~/bin/ideel/snakemake
fi

if $annotate; then
	mkdir $directory/analyses/annotations/gffs/
	if $trimmed; then
		COUNTER=1
			for trimmed_sample in $directory/barcodes/trimmed.BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many annotations to produce.
			f=$(basename $trimmed_sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			prokka --outdir $directory/analyses/annotations/$location --force --prefix trimmed.$location --cpus $threads $directory/samples/$location/trimmed_unicycler_hybrid/assembly.fasta
			cp $directory/analyses/annotations/$location/*gff $directory/analyses/annotations/gffs/
			COUNTER=$[$COUNTER +1]
		done
	else
		COUNTER=1
		for trimmed_sample in $directory/barcodes/BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many annotations to produce.
			f=$(basename $trimmed_sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			prokka --outdir $directory/analyses/annotations/$location --force --prefix $location --cpus $threads $directory/samples/$location/unicycler_hybrid/assembly.fasta
			cp $directory/analyses/annotations/$location/*gff $directory/analyses/annotations/gffs/
			COUNTER=$[$COUNTER +1]
		done
	fi
	
fi

if $pangenome; then
	roary -e --maft -p $threads -f $directory/analyses/annotations/gffs/roary $directory/analyses/annotations/gffs/*.gff
fi

if $auto_shutdown; then
	shutdown -h now


fi
