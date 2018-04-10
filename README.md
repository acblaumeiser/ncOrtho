## ncOrtho  
NcOrtho is a software to predict orthologous microRNAs within a 
species of interest for a reference microRNA precursor sequence. Based on
the reference sequence and a user defined set of species a core orthologous 
set of microRNA precursor sequences is constructed based on gene order. 
This set is used to construct a covariance model which is in turn used
to search your species of interest.
For workflow details see the user manual (UserManual_ncOrtho.pdf).

## Installation
The tool comes in two Perl scripts: the ncOrtho1.0.0_pre.pl script to 
precompute/hash your input data and the ncOrtho1.0.0_main.pl script, which
implements the prediction algorithm. Both scripts rely on secondary programs
which must be present on your system. Also you have to specify the respective
path in the head of both scripts. 

List of programs needed:
- Infernal-package (cmbuild, cmcalibrate, cmsearch)
- blastn and makeblastdb
- t-coffee

## Usage
Prerequisite for the ncOrtho-1.0.0_pre.pl script is that all required input data
is present in the following directory structure:

data/root/genome
data/root/gtf
data/core/genome
data/core/gtf
data/oma
data/miRNAs
data/interest

The precomputing script can now be called in the following way:

/home/homer/ncOrtho/ncOrtho-1.0.0_pre.pl -root_genome /home/homer/ncOrtho/example/data/root/genome/Homo_sapiens.GRCh38.dna.toplevel.fa.reduced_chr9 -root_gtf /home/homer/ncOrtho/example/data/root/gtf/Homo_sapiens.GRCh38.76.gtf -core_genome_folder /home/homer/ncOrtho/example/data/core/genome/ -core_gtf_folder /home/homer/ncOrtho/example/data/core/gtf/ -oma_ortho_folder /home/homer/ncOrtho/example/data/oma/

Next the main script can be called by typing:

/home/homer/ncOrtho/ncOrtho-1.0.0_main.pl -root_genome /home/homer/ncOrtho/example/data/root/genome/Homo_sapiens.GRCh38.dna.toplevel.fa.reduced_chr9 -root_gtf_hash_file /home/homer/ncOrtho/example/data/root/gtf/Homo_sapiens.GRCh38.76.gtf.hash -core_genome_hash_file /home/homer/ncOrtho/example/data/core/genome/core_genome.hash -core_gtf_hash_file /home/homer/ncOrtho/example/data/core/gtf/core_gtf.hash -oma_hash_file /home/homer/ncOrtho/example/data/oma/oma_ortho.hash -nc_rna /home/homer/ncOrtho/example/data/miRNAs/hsa-let-7a-1.fa -interest_genome /home/homer/ncOrtho/example/data/interest/Mus_musculus.GRCm38.dna.toplevel.fa_chr13 -out /home/homer/ncOrtho/example/search/test

## History
first release in version 1.0.0

## Credits
Dept. for Applied Bioinformatics
Institute for Cell Biology and Neurosciences
Goethe University Frankfurt am Main
