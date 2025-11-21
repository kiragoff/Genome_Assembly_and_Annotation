#!/bin/bash
set -e
##not actually an executable script, but my qiime2 documentation
##DOCUMENTATION for assembly of NA isolates sequenced by Quebec collaborators
##ASSEMBLY OF NANOPORE/MINION/LONG READ/HIGH MOLECULAR WEIGHT SEQUENCING

##GENERAL NOTES
##I am working from a folder called NAIsolates. I have my reads in a folder inside of it named 'BaseCalled'
##The name of the initial folder isn't important, just know that this is your base folder you'll be moving in and out of
##If you name your folder with the basecalled reads something else, just make sure you replace that everywhere else here I'm pulling from the BaseCalled folder
##I know from having done this assembly before that these are over-sequenced: sequencing depth ranges from 300-800x. I've gotten clean assemblies with 20x from different sequencing runs.
##Basecalling was done in Guppy in Super Accurate Mode (SUP)
##Basecalling with Dorado is more accurate - if I had the raw data, I would re-basecall with Dorado

#Re-basecall, if needed! Ideally things are called in Dorado in super accurate mode
#This requires the pod5 files
#Here I'm feeding it an array of barcodes to use
mkdir basecalled
cd basecalled
for f in {11..24}
do
	/path/to/dorado basecaller --emit-fastq sup ../Pod5/barcode$f/ | gzip > barcode$f.fastq.gz
done

#If you didn't re-basecall and gzip above, gzip any fastq files to save dramatic amounts of space
cd basecalled
gzip *.fastq
#OR if you have pigz installed (multi-threads gzipping)
#pigz *.fastq

##RENAME
#Rename files to something short but relevant 
#If I had isolate numbers associated with the barcodes, I would use the isolate codes instead of bcx
#You can skip this step but I find it makes things easier to deal with
#This is mostly here in case you start with a file named something like GROW_STRAINS_UofC_20240119_MP_barcode21.fastq.gz
for f in {11..24}
do
	mv barcode$f.fastq.gz bc$f.fastq.gz
done

## QC TRIMMING
#We are using Porechop-abi https://github.com/bonsai-team/Porechop_ABI
#This trims adapters and barcodes from the sequences
#The abi flag (-abi) means porechopabi is guessing the adapters from the reads themselves, rather than us supplying them
#The discard middle flag (--discard_middle) finds sequences with adapters in the middle and chops them into separate reads/treats them as chimeric sequences
#If you run into memory issues, add a threads flag (ie -t 9) to limit number of threads
#Use individual calls if you're adding just one thing to the folder, etc

mkdir ../qctrim
cd ../qctrim/
conda activate porechopabi
for f in {11..24}
do
	porechop_abi -abi --discard_middle -i ../basecalled/bc$f.fastq.gz -o bc$f.qctrim.fastq
done
conda deactivate
pigz *.fastq


##QC DATA
#Basic quality check with FastQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#If you get an "out of memory error" use '-t x' 

mkdir ../fastqc
conda activate fastqc
fastqc -t 9 -o ../fastqc/ *.fastq.gz
conda deactivate

#look at your output

##REMOVE POOR QUALITY READS, SHORT READS, AND DNA CONTROL SEQUENCES	
#get rid of low quality (q<=10) reads, dna control sequences (dcs) and length <=1000 with Chopper
#https://github.com/wdecoster/chopper
#Increase your length filter above 1000 if your fragment distribution allows
#make sure DCS.fasta is in base folder (if using relative path) or databases (if using absolute path)
#Add -t flag if memory errors: time chopper -t 9 -c /mnt/sdb/databases/DCS.fasta -q 10 -l 1000 -i $f > $base.clean.fastq
#This limits us to using 9 threads
	
mkdir ../clean
cd ../clean

conda activate chopper
for f in {11..24}
do
	chopper -t 9 -c /mnt/sdb/databases/DCS.fasta -q 10 -l 1000 -i ../qctrim/bc$f.qctrim.fastq.gz > bc$f.clean.fastq
done
conda deactivate
pigz *.fastq



CONTINUE UPDATING SCRIPT FROM HERE, NOTE CHANGES IN FOR LOOPS ARRAYS, FOLDER NAMES, ADD PLASMID MODE TO FLYE

##ASSEMBLY
##Assemble with Flye https://github.com/fenderglass/Flye
##First pass is meta flag on everything
mkdir ../Assembly
cd ../Assembly
conda activate flye
for f in ../CleanData/*.clean.fastq
do
	base="$(basename $f .clean.fastq)";
	echo "\n\n$base\n\n"
	time flye --meta --threads 9 --nano-hq $f -o ./$base.meta.9thread
done

#Flye seems to have trouble with extremely high coverage when run in single-genome/non-meta mode.
#Previously only BC21 assembled with normal script (328x vs 866x, 814x, and 622x) in single genome mode
#ERROR: No disjointigs were assembled - please check if the read type and genome size parameters are correct
#There is a step that filters out repetitive contigs - with the higher coverage 
#--asm-coverage 50 tells flye to keep only the longest contigs with at least 50x coverage; requires an estimate of genome size
#If you have a normal level of coverage, just rerun the for loop after removing the meta flag (--meta) 
#and changing the output name so it doesn't overwrite your previous assembly (./$base.meta.9thread -> ./$base.single.9thread or something like that)
##If you prefer, you can also use BBNorm from BBTools to decrease the numbers of reads prior to assembly


flye --threads 9 --asm-coverage 50 --genome-size 5300000 --nano-hq ../CleanData/bc21.clean.fastq.gz -o bc21.single.9thread
flye --threads 9 --asm-coverage 50 --genome-size 3100000 --nano-hq ../CleanData/bc22.clean.fastq.gz -o bc22.single.9thread
flye --threads 9 --asm-coverage 50 --genome-size 5800000 --nano-hq ../CleanData/bc23.clean.fastq.gz -o bc23.single.9thread
flye --threads 9 --asm-coverage 50 --genome-size 5400000 --nano-hq ../CleanData/bc24.clean.fastq.gz -o bc24.single.9thread

##Single genome mode assumes sequencing depth for all contigs, meta mode does not. 
##Running both is a good way to check contamination and extrachromosomal elements.
##I recommend running in meta mode, then in single genome mode when you've confirmed that you are looking at a single isolate
##Contamination can look like: multiple contigs (especially circularized) > 1mbp or multiple contigs in the 100s of kbp but nothing over 1mbp
##Differences in GC ratios can also be a flag
##In meta mode, you might see a large number of non-circularized, shorter contigs - the higher the genome coverage the more you're like to see
##These are usually mis-assemblies - you can BLAST them or map them back against the main chromosome to check, LCB alignments are also a good check
##Or you can see if they disappear upon single-genome assembly
##Running with more than one thread causes some heuristic effects, especially for complex metagenomes. 
#For single genomes, after the Medaka polishing step, OrthoANIu comparison of assemblies done with 1 thread or 9 show ANI values >= 99.90%, indicating the genomes are functionally identical
##There can be a few more differences in the extra chromosomal elements, however
##Even Kit 14 v10 chemistry basecalled in Doradao in super accurate mode has a 0.05% error rate - expect it to be higher with different chemistry, different accuracy of basecalling, or Guppy instead of Dorada



mkdir dir1 dir2 etc

##Collect and rename assembly files
mkdir ../AllAssemblies
cd ../AllAssemblies

for f in ../Assembly/*
do
	base=$(basename $f)
	cp $f/assembly.fasta $base.assembly.fasta
done


## Copying and renaming assembly graphs into one folder
mkdir ../AssemblyGraphs
cd ../AssemblyGraphs
for f in ../AllAssemblies/*.fasta
do
	base=$(basename $f .assembly.fasta)
	cp ../Assembly/$base/assembly_graph.gfa $base.assembly_graph.gfa
done

##Bandage is good tool for viewing the assembly graphs
##After opening, go File -> Load Graph
##Click 'Draw graph' on the graphic interface
##You can click it again and it will visualize differently
##Lets you look at what's circularized, anything that's got loops, optional overlays for depth, etc
#Paper: https://academic.oup.com/bioinformatics/article/31/20/3350/196114





##GENOME POLISHING
##Genome Polishing with Medaka
##This requires two loops because it's using the same set of initial reads (in the BaseCalled folder) to polish the meta and non-meta assemblies
##Three loops with the manual ECE removal I did here to check our results against the Quebec pipeline
##Polishing looks for local misassemblies and other errors and corrects them
##With long read only sequencing like this, we're using the original reads to polish the assemblies
##The gold standard has long been using short read sequencing to polish long read assemblies
##Short read sequencing has a much lower error rate
##But current chemistry and basecallers have dramatically lowered the error rate in nanopore long-read sequencing
##This polished file should be used for all downstream analysis


mkdir ../Mapping
mkdir ../medaka
cd ../medaka
export TF_FORCE_GPU_ALLOW_GROWTH=true
conda activate medaka
for f in ../BaseCalled/*.fastq.gz
do
	base=$(basename $f .fastq.gz)
	time mini_align -i ../BaseCalled/$base.fastq.gz -r ../Assembly/$base.meta.9thread/assembly.fasta -P -m -p ../Mapping/$base.meta.9thread -t 8
	time samtools index ../Mapping/$base.meta.9thread.bam
	time medaka consensus ../Mapping/$base.meta.9thread.bam $base.meta.9thread.hdf --threads 8 --model r1041_e82_400bps_sup_v4.3.0 --batch_size 10
	time medaka stitch $base.meta.9thread.hdf ../Assembly/$base.meta.9thread/assembly.fasta $base.meta.9thread.polished.assembly.fasta
done

for f in ../BaseCalled/*.fastq.gz
do
	base=$(basename $f .fastq.gz)
	time mini_align -i ../BaseCalled/$base.fastq.gz -r ../Assembly/$base.single.9thread/assembly.fasta -P -m -p ../Mapping/$base.single.9thread -t 8
	time samtools index ../Mapping/$base.single.9thread.bam
	time medaka consensus ../Mapping/$base.single.9thread.bam $base.single.9thread.hdf --threads 8 --model r1041_e82_400bps_sup_v4.3.0 --batch_size 10
	time medaka stitch $base.single.9thread.hdf ../Assembly/$base.single.9thread/assembly.fasta $base.single.9thread.polished.assembly.fasta
done

conda deactivate

##In this case, because I am looking at the effects of the Quebec pipeline, I am also doing manual removal of the extraneous extrachromosomal elements from the polished metagenomic assembly and running them through the rest of the pipeline
##To do manual removal of ECE from meta assemblies without overloading text editor (will only happen for some files)
##Open fasta file of meta assembly with sublime, View -> uncheck Word Wrap -> Delete the ECE you want to get rid of (use All_assembly_info.txt to figure out what to pull and what to leave)
##Basically, keep anything circularized or with high coverage


##Gather all polished assemblies

mkdir AllPolished
cd medaka
cp /medaka/*.fasta /AllPolished






##GENOME COMPLETENESS AND CONTAMINATION
##CheckM2 for genome completeness and contamination
mkdir ../CheckM2
cd ../CheckM2
conda activate checkm2
time checkm2 predict --input ../medaka --output-directory . -t 6 --allmodels -x fasta 

for f in ../medaka/*.polished.assembly.fasta
do
	base=$(basename $f .polished.assembly.fasta)
	checkm2 predict --input $f --output -t 6 --allmodels -x fasta 
done

##What do these measures mean?
##It looks at what % of marker genes are present
##95% completeness means that 95% of the marker genes CheckM2 is looking for are present
##30% contamination means that 30% of the marker genes have more than one copy, 100% means there's multiples of everything
##I'm usually happen as long as contamination <10%
##Original CheckM2 paper uses three categories for MAGs: 
	##high quality: >90% complete, <5% contaminated
	##medium quality: 50%-90% complete, <10% contaminated
	##low quality: <50% complete, <10% contaminated)
##We expect better quality measures from single-genome sequencing
##Generally, once a genome is circularized, it's complete, no matter what % of marker genes have been found


##Copy the CheckM2 document into the main folder so that the assembly information, sourmash, and quality info are in the same place

cp CheckM2/quality_report.tsv .





##TAXONOMY
##Sourmash for taxonomy 
##https://sourmash.readthedocs.io/en/latest/
##The representative sequence database is also in the database folder


mkdir ../SourmashRep
cd ../SourmashRep
conda activate sourmash
for f in ../medaka/*.polished.assembly.fasta
do
	base=$(basename $f .polished.assembly.fasta)
	echo 
	echo "########################$base"
	echo
	time sourmash sketch dna -p k=51 $f -o $base.sig.gz --name $base
	time sourmash gather $base.sig.gz /mnt/sdb/databases/gtdb-rs214-k51.zip --save-matches $base.matches.zip -o $base.csv
	time sourmash tax annotate --gather-csv $base.gather.csv --taxonomy-csv /mnt/sdb/databases/gtdb-rs214.lineages.csv
done
conda deactivate

##Showing here with representative database and full database, just for reference

mkdir ../Sourmash
cd ../Sourmash
conda activate sourmash
for f in ../medaka/*.polished.assembly.fasta
do
	base=$(basename $f .polished.assembly.fasta)
	echo 
	echo "########################$base"
	echo
	time sourmash sketch dna -p k=51 $f -o $base.sig.gz --name $base
	time sourmash gather $base.sig.gz /mnt/sdb/databases/gtdb-rs214-reps.k51.zip --save-matches $base.matches.zip -o $base.csv
	time sourmash tax annotate --gather-csv $base.gather.csv --taxonomy-csv /mnt/sdb/databases/gtdb-rs214.lineages.csv
done
conda deactivate

##Gather all SourmashRep taxonomy into one file

cd ..
cat SourmashRep/*csv > All_SourmashRep.csv

##If you are opening in LibreOffice, make sure that upon import you check "comma" as a separator option or it will import in one big chunk (csv means comma separated values, tsv is tab separated)

cd ..
cat SourMash/*csv > All_Sourmash.csv



##ALTERNATE TAXONOMY APPROACHES
##SourMash is a time consuming step. Alternatively - or additionally, because they do taxonomy in different ways - you can upload the polished assembly to TYGS
## https://tygs.dsmz.de/user_requests/new
##This will give you whole-genome based taxonomy, with pairwise comparisons to reference genomes - I find this to be useful
##You can also request a full 16S report, but it doesn't do it by default
##You don't need to specify reference genomes
##TYGS will email you when it's done, and you can take a look at phylogenetic trees and pairwise comparisons to other genomes. 
##The number we care about here is the dDDH d4 value: if this is <70%, it's a potential new species
##GTDB-tk is a standard, but it's too memory intensive to do on this computer. As of 2024, it is very possible to publish genomes with the TYGS speciation.
##GTDB-Tk can be run on the Alliance Canada servers, or on ARC
##If you are concerned about contamination based on the meta mode assembly results, try running your meta assembly through Lemur and see what signatures it finds for each contig





##ANNOTATION
##Annotation with Bakta
##When you upload genomes to NCBI for publication, the genomes will be re-annotated with PGAP.
##PGAP also requires more processing power than we have. Bakta is a good compromise.
##Prokka is also a (more common) option, but Prokka is no longer maintained, and Bakta provides slightly better annotation with lower processing requirements
##Bakta is also available as an online tool, and you can upload polished assemblies to https://bakta.computational.bio/
##However, you have to download and rename each individual file for each individual assembly, as well as copying the statistics data into a new file

mkdir ../Bakta
cd ../Bakta
conda activate bakta

for f in ../medaka/*polished.assembly.fasta
do
	base=$(basename $f .polished.assembly.fasta)
	bakta --db /mnt/sdb/databases/db --output ./$base $f
done
conda deactivate


##The path variable doesn't want to stay exported, so we're using the database flag
##If there are memory issues or panic, add a thread flag

##For each genome, Bakta outputs 
	## A .txt file of summary statistics
	## png and svg of the annotate assembly
	## fna, faa, ffn, gbff, gff3, embl, and tsv files
	## separate faa and tsvs of all hypothetical proteins
	## several detailed logs

##To gather your selected format together, use a copy loop 
##Here I'm pulling all of the statistics together by copying all the text files (.txt) into a new folder
##This is optional but I find it more convenient
##If you have files in the base Bakta folder, or non-Bakta folders, the command line will throw and error
##Everything will copy just fine despite that


mkdir ../BaktaStats
cd ../BaktaStats

for f in ../Bakta/*
do
	base=$(basename $f)
	cp $f/$base.polished.assembly.txt $base.bakta.stats.txt
done


##Repeat as desired to pull each type of file you want 
##$base.bakta.stats.txt changes the name of the copied file to [foldername].bakta.stats.txt
##When you're changing the file type you're moving, change this as well.
##For example, relocating and renaming all the gff3 files:

mkdir ../BaktaGFF3
cd ../BaktaGFF3

for f in ../Bakta/*
do
	base=$(basename $f)
	cp $f/$base.polished.assembly.gff3 $base.polished.bakta.gff3
done

##And for GBFFs, FAAs, FFNs

mkdir ../BaktaGBFF
cd ../BaktaGBFF
for f in ../Bakta/*
do
	base=$(basename $f)
	cp $f/$base.polished.assembly.gbff $base.polished.bakta.gbff
done

mkdir ../BaktaFAA
cd ../BaktaFAA

for f in ../Bakta/*
do
	base=$(basename $f)
	cp $f/$base.polished.assembly.faa $base.polished.bakta.faa
done

mkdir ../BaktaFFN
cd ../BaktaFFN

for f in ../Bakta/*
do
	base=$(basename $f)
	cp $f/$base.polished.assembly.ffn $base.polished.bakta.ffn
done


##Move pngs and svgs into the same folder

mkdir ../BaktaVis
cd ../BaktaVis
for f in ../Bakta/*
do
	base=$(basename $f)
	cp $f/$base.polished.assembly.svg $base.polished.bakta.svg
	cp $f/$base.polished.assembly.png $base.polished.bakta.png
done


##Bakta isn't great for finding CRISPRs
##You can use this online tool instead: https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Index






##ADDITIONAL CONSIDERATIONS
##If you are trying to figure out if you have sequenced the same organism multiple times, or if they are two closely related strains, do an ANI comparison of their genomes
## https://www.ezbiocloud.net/tools/ani
## Can only do pairwise comparisons, but uses the OrthoANI algorithm. 
##>=95% often considered the same species. 
##>=99.9% is an arbitrary cutoff I've been using as a strain differentiator
##Even with current chemistry and basecalling, there's still a 0.05% error rate, and there are some heuristic impacts of multithreading the assembly process
##Whether they are actually different strains will depend on where the differences are - if they're impacting genes, etc





Overview by Kira Goff, July 2024
