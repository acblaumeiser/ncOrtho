#! /usr/bin/perl
## Author : Daniel Amsel - daniel.amsel@gmx.de
## Modified - Ingo Ebersberger
## Modified by Andreas Blaumeiser
## Last modified: 04/10/2018

use strict;
use warnings;
use Getopt::Long;
use Storable;
#use Bio::DB::Fasta;

##################################################################################
############################# PATHVARIABLES TO EDIT ##############################
##################################################################################

# Infernal
my $cmbuild     = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmbuild';			# path or command to cmbuild
my $cmcalibrate = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmcalibrate';		# path or command to cmcalibrate
my $cmsearch    = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch';		# path or command to cmsearch

# BLAST+
my $blastn      = 'blastn';		# path or command to blastn
my $formatdb	= 'makeblastdb';	# path or command to makeblastdb

# T-COFFEE
my $tcoffee     = 't_coffee';		# path or command to t_coffee

##################################################################################
######################## DO NOT TOUCH THE CODE BELOW THIS LINE ###################
##################################################################################

###############################
####### INPUT VARIABLES #######
###############################
my $root_genome;             ##
my $root_gtf_hash_file;      ##
###############################
my $core_genome_hash_file;   ##
my $core_gtf_hash_file;      ##
###############################
my $oma_hash_file;	     ##
###############################
my $ncRNA;                   ##
my $rna_start;               ##
my $rna_stop;                ##
my $rna_chr;                 ##
###############################
my $ukn_genome;              ##
###############################
my $outpath;                 ##
###############################
my $createCM = 1;            ##
###############################
my $max_intra_prot = 0;	     ##
my $min_seq_len = 0.9;	     ##
###############################
my $help;                    ##
if (@ARGV==0) {              ##
        $help = 1;           ##
}                            ##
###############################
my $outfile = 1;	     ##
my $cpu = 8;	             ##
###############################

my $helpmessage = 	"
===============================================================================
| 			NON CODING RNA SEARCH - De'wI'	 		      |
|=============================================================================|
|									      |	
|	-root_genome <> 		=	path/to/root/genome.fa	      |		
|	-root_gtf_hash_file <>		=	path/to/root/gtf_hash_file    |
|									      |	
|	-core_genome_hash_file <>	=	path/to/core/genome_hash_file |
|	-core_gtf_hash_file <>		=	path/to/core/gtf_hash_file    |
|									      |	
|	-oma_hash_file <>		=	path/to/oma/db_hash_file      |
|									      |	
|	-nc_rna <>			=	path/to/query/nc_rna.fa       |
|									      |	
|	-interest_genome <>		=	path/to/interest/genome.fa    |
|									      |
|	-out|outpath <>			=	path/to/output                |
|									      |	
|-----------------------------------------------------------------------------|
|-----------------------------------------------------------------------------|
|									      |
| OPTIONAL:								      |
|									      |	
|	-rna_start <>			=	start position of ncRNA	      |
|	-rna_stop <>			=	stop  position of ncRNA       |
|	-rna_chr <>			=	chromosome     of ncRNA       |
|	-cpu <>				=	nr of cores		      |	
|	-mip <>				=	max. inserted proteins        |
|		maximum number of allowed proteins between up- and down-      |
|		stream protein sequences - e.g. <0> or <2> or ...             |
|	-msl <>				=	min. sequence length          |
|		minimum sequence length of query RNA that must be found       |
|		by BLASTN - e.g. <0.123> or  <0.9> or ...                     |
| 									      |
|	-h|help				=	prints this help              |
|									      |
|     !!!!! Please note to provide the full paths to the algorithm !!!!!      |
|=============================================================================|
===============================================================================\n";

#################################################################################################
GetOptions      (       'root_genome=s'                 =>      \$root_genome                   ,
                        'root_gtf_hash_file=s'          =>      \$root_gtf_hash_file            ,

                        'core_genome_hash_file=s'       =>       \$core_genome_hash_file	,
                        'core_gtf_hash_file=s'          =>      \$core_gtf_hash_file            ,

			'oma_hash_file=s'		=>	\$oma_hash_file			,		
	
                        'nc_rna=s'                      =>      \$ncRNA                         ,

                        'interest_genome=s'             =>      \$ukn_genome                    ,

                        'out|outpath=s'                 =>      \$outpath                       ,

                        'rna_start=s'                   =>      \$rna_start                     ,
                        'rna_stop=s'                    =>      \$rna_stop                      ,
                        'rna_chr=s'                     =>      \$rna_chr                       ,
		
			'mip=i'				=>	\$max_intra_prot		,
			'msl=f'				=>	\$min_seq_len			,
			'cpu=i'				=>	\$cpu				,
			
			'outfile=i'			=>	\$outfile			,
			
                        'h|help'                        =>      \$help                          )

                        || die "Error in command line arguments!\n$!\n"                         ;
#################################################################################################

###########################
if ($help) {
        print $helpmessage;
        exit;
}

###############################
## checking for interest genome
if (!-e "$ukn_genome"){
	die "Error while checking for query genome.\nCould not find: $ukn_genome\nExiting...\n";
}

###########################
## checking for outdir
if (!-e "$outpath"){
	print "$outpath does not exist... Trying to generate it.\n";
	!`mkdir -p $outpath` or die "Could not generate outdir $outpath.\n"; 
}

chdir $outpath; # change to specified output folder

###########################
###########################
## specifying the output file names
my $alignment 		= $outpath."/seq.aln";
my $stockholm_aln 	= $outpath."/rna_aln.sto";
my $covariance_model 	= $outpath."/rna.cm";
my $cmsearch_out = $outpath."/cmsearch.out";
###########################
################################################
## checking whether a covariance model already exists for this miRNA
if (-e "$outpath/rna.cm" and -e "$outpath/seq.aln" and -e "$outpath/rna_aln.sto"){ 
	print "A covariance model already exists in $outpath.\nSkipping the covariance model generation.\n";
	$createCM = 0;

}
###########################
###################################################################
my %root_gtf_hash               = %{retrieve($root_gtf_hash_file)};
my %core_gtf_hash               = %{retrieve($core_gtf_hash_file)};
my %core_genome_hash		= %{retrieve($core_genome_hash_file)};
my %oma_hash			= %{retrieve($oma_hash_file)};
###################################################################        

#################################################################################################
################### create output folder for nc_RNA ukn genome combination ######################
#################################################################################################
my @query_ncrna_name_split      = split('/',$ncRNA);
my $query_ncrna_name            = $query_ncrna_name_split[-1];
##################################################################################
my @ukn_genome_name_split       = split('/',$ukn_genome);
my $ukn_genome_name             = $ukn_genome_name_split[-1];
##################################################################################
       
my $ncRNA_seq;
open(NCRNA,"<",$ncRNA);
while(<NCRNA>){
        next if (/^>/);
        $ncRNA_seq.=$_;
}
close(NCRNA);
my $ncRNA_seq_len = length($ncRNA_seq);

#################################################################################################
################# BLASTN search for ncRNA if no start, stop, chr defined by user ################
#################################################################################################
my $rna_name;
my $rna_strand;
my $tmp_start;
my $tmp_stop;
###################################################################################
if ((not defined $rna_start) or (not defined $rna_stop) or (not defined $rna_chr)){
        my $ncPosBLASTout = "$outpath/ncRNA_rootGenome.blast.out";
        if ((-e $ncRNA) && (-e $root_genome)){
		#system("$blastn -p blastn -d $root_genome -i $ncRNA -o $ncPosBLASTout -m 8");
                system("$blastn -task blastn -db $root_genome -query $ncRNA -out $ncPosBLASTout -outfmt 6 -num_threads $cpu");
        }
       ###################################################################################
       if (-z $ncPosBLASTout){
                die print STDERR "No sequence found in BLASTn search - Position finding not possible for $ncRNA\n";
        }
###################################################################################
	my $rna_score;
        ($rna_start,$rna_stop,$rna_chr,$rna_name,$rna_strand,$rna_score)	= &blast_parser($ncPosBLASTout);
}
else{
	my @tmp_rna_name = split('/',$ncRNA);
	$rna_name = $tmp_rna_name[-1];
}
print ">$rna_name :\n";
print "Found RNA at $rna_start, $rna_stop, $rna_chr\n";
print STDERR ">$rna_name :\n";

if ($createCM) { ### a covariance model does not yet exist for the miRNA

###################################################################################
#################### SEARCH UP AND DOWN STREAM SEQUENCE OF NCRNA ##################
###################################################################################

# |---DS---|=============|---ncRNA----|=========|---US---| #


	print "Filtering chromosomes and sorting sequences by start position in ROOT GTF file...";
	my @root_keys;
	my %chr_hash;
# 
	foreach(keys %root_gtf_hash){
		my @tmp_array = @{$root_gtf_hash{$_}};  # (47365496 47384273 10 ENST00000585316 ENSG00000265763)
		if 	($tmp_array[2] eq $rna_chr){
			if (not exists $chr_hash{$tmp_array[4]}) {
				$chr_hash{$tmp_array[4]}=\@tmp_array;
			}	
		}
		@root_keys = sort { ${$chr_hash{$a}}[0] <=> ${$chr_hash{$b}}[0]} keys(%chr_hash);
	}
	print "done\n";
	my @US=(); # array of all Upstream genes of ncRNA
	my @DS=(); # array of all Downstream genes of ncRNA
	my @MS=(); # ncRNA matches within Protein ( e.g. intron )
	
	foreach(@root_keys){
		my $ensg = $_;
		my @tmp_array = @{$chr_hash{$ensg}};

		if 	($tmp_array[0] > $rna_stop){
			#print "Stop position of intergenic region $tmp_array[0]\n";
			push (@US,\@tmp_array);
		}
		elsif 	($tmp_array[1] < $rna_start){
			#print "Start position of intergenic region $tmp_array[1]\n";
			push (@DS,\@tmp_array);
		}
		elsif 	(($tmp_array[0] <= $rna_start) && ($tmp_array[1] >= $rna_stop)){
			#print "+++++ INTRON HIT $rna_name : @tmp_array\n";  # debug
			push (@MS,\@tmp_array);
		}
	}
	my $rna_file = $outpath."/core_orthos.fa";
	my $rna_file_seq = "";
	# a test that prints the number of upstream annotations
	#my $test = scalar(@US);
	#print "$test\n";
	if ((scalar(@US) == 0) && (scalar(@MS)== 0)){ 
		print "ERROR : Did not find Upstream Sequence! Terminating\n"; 
		# in the root genome, no gene was found upstream of the miRNA position, therefore a shared syntenic region cannot be constructed
		exit;
	}	
	elsif ((scalar(@DS) == 0) && (scalar(@MS)== 0)){
		print "ERROR : Did not find Downstream Sequence! Terminating\n"; 
		# in the root genome, no gene was found downstream of the miRNA position, therefore a shared syntenic region cannot be constructed
		exit;
	}	

	# routine for the case that the miRNA of the root genome is located within a gene
	elsif (scalar(@MS)> 0){
		my %MS_hash;
		print "ncRNA is located within PROTEIN\n";
		foreach(keys %oma_hash){
			my $species = $_;
			my %tmp_hash = %{$oma_hash{$species}};
			my $ms_index=0;		# if more than one gene is in @MS

			if (exists $tmp_hash{${$MS[$ms_index]}[4]}){
				print "$species has ortholog of ${$MS[$ms_index]}[4]\n";
				$MS_hash{$species} = $tmp_hash{${$MS[$ms_index]}[4]};
				print "==== @{$tmp_hash{${$MS[$ms_index]}[4]}}\n";
			}
			else{
				$ms_index+=1;
				my $found_ms_seq = 0;
				my $ms_size = @MS-1;
				while (($ms_index <= $ms_size) && ($found_ms_seq == 0)){
					if (exists $tmp_hash{${$MS[$ms_index]}[4]}){
						$MS_hash{$species} = $tmp_hash{${$MS[$ms_index]}[4]};
						$found_ms_seq = 1;
					}
					$ms_index+=1;
				}
				if (not exists $MS_hash{$species}){
					print "ncRNA within Protein: No ortholog found for $species\n";
					delete $core_gtf_hash{$species};
                   	delete $core_genome_hash{$species};
				}
			}
		}
		# for each potential coortholog, a blast run is done, then the one with the highest score is selected
		foreach (keys %MS_hash){
			my $species = $_;
			my @tmp_array = @{$MS_hash{$species}};
			print "co ortholog array :@tmp_array\n";
			my %tmp_hash = %{$core_gtf_hash{$species}};
			my $core_genome_file = $core_genome_hash{$species};
			print "$core_genome_file\n";
			my $core_rna_start;
			my $core_rna_stop;
			my $core_rna_chr;
			my $core_rna_strand;
			my $core_rna_score = 0;
			my $highscore_outfile = "";
    	   
			foreach(@tmp_array){
				my $ortholog = $_;
				my @pos = @{$tmp_hash{$ortholog}};
				print "positions: @pos\n"; # 43620095 43652148 1 ENSMMUT00000014574 ENSMMUG00000014970
               	# parse intergenic sequence between up and downstream   
                my $inter_seq   = &genome_parser($core_genome_file,$pos[0],$pos[1],$pos[2],2);
				my $len_seq = length($inter_seq);
               	print "INTERSEQ-LEN: $len_seq\n";
	
	            my $outfile     = $outpath."/".$species."_".$pos[4].".interseq"; # G.gorilla_ENSGGOG00000034972.interseq
	            open(OUTPUT,">",$outfile);
	            print OUTPUT ">$pos[2]\n$inter_seq";
	            close(OUTPUT);
			
				# search for the initial ncRNA within the identified ortholog PCT sequences of each core species
	            system("$formatdb -dbtype nucl -in $outfile"); 
        	   	my $interseq_blast_out = $outfile.".blastout";
				#system("$blastn -p blastn -d $outfile -i $ncRNA -o $interseq_blast_out -m 8"); #perform the blastn search
#                system("$blastn -db $root_genome -query $ncRNA -out $ncPosBLASTout -outfmt 6");
				system("$blastn -task blastn -db $outfile -query $ncRNA -out $interseq_blast_out -outfmt 6 -num_threads $cpu"); #perform the blastn search

	            if (-z $interseq_blast_out){
	            	print "No BLAST HIT found. Skipping\n";
					next;
	            }
	            my ($pred_rna_start,$pred_rna_stop,$pred_rna_chr,$pred_rna_name,$pred_rna_strand,$pred_rna_score) = &blast_parser($interseq_blast_out);

				# incooporate length criterion on the matched region by the blast search
				if (($pred_rna_stop-$pred_rna_start)<($ncRNA_seq_len*$min_seq_len)){
               		print STDERR "§§§§ The potential ortholog of $species was too short.\n";	###!!!### dies here
               	    next;
               	}
				if ($core_rna_score==0){
					$core_rna_start		= $pred_rna_start;
					$core_rna_stop		= $pred_rna_stop;
					$core_rna_chr		= $pred_rna_chr;
					$core_rna_strand	= $pred_rna_strand;
					$core_rna_score		= $pred_rna_score;
					$highscore_outfile 	= $outfile;
				}
				else{
					if($pred_rna_score > $core_rna_score){
				        $core_rna_start         = $pred_rna_start;
                        $core_rna_stop          = $pred_rna_stop;
                        $core_rna_chr           = $pred_rna_chr;
                	    $core_rna_strand        = $pred_rna_strand;
		                $core_rna_score         = $pred_rna_score;
						$highscore_outfile	= $outfile;
					}
				}
			}
			if(not defined $core_rna_score){
		       	print STDERR "No BLASTN HIT found , skipping $species!\n";
			}
			if ($highscore_outfile eq ""){
				next;
			}
            my $pred_rna = &genome_parser($highscore_outfile,$core_rna_start,$core_rna_stop,$core_rna_chr,$core_rna_strand);
			$rna_file_seq.=">$species\n$pred_rna\n";
		}
	}

	# the shared syntenic region does exist in the root genome, now we check if it is present in each of the core species
	# the value in MIP is the number of insertions which are allowed within the syntenic region of the core genome
	else{	
		my %US_hash;
		my %DS_hash;
		print "ncRNA in intergenic region, construct shared syntenic region: US--------miRNA-----------DS\n";
		#my $free_insertions = $max_intra_prot; # allowed insertions of proteins between up and downstream - default 0
		#print "INITIAL NUMBER OF FREE INSERTIONS: $free_insertions\n";
		foreach(keys %oma_hash){ # wird für jede species ausgeführt
			my $us_index = 0;
			my $ds_index = -1;
			my $species = $_;
			print "###################################################################################\n";
			print "######### Construction of shared syntenic region for species:\t$species \n";
			print "###################################################################################\n";
			my %tmp_hash = %{$oma_hash{$species}}; # contains oma orthologs for that species
			my $free_insertions = $max_intra_prot; # allowed insertions of proteins between up and downstream - default 0 ---------------- Mirko edit
			print "INITIAL NUMBER OF FREE INSERTIONS FOR SPECIES $species:\t $free_insertions\n"; #----------------- Mirko edit
		
			print "initial US gene of root species: ${$US[$us_index]}[4] \n";

##################################################################################################
# check the upstream region of the miRNA 
##################################################################################################

			# if the hit is found right away
			if (exists $tmp_hash{${$US[$us_index]}[4]}){
				print "$species: ortholog exists in the first place, no insertions needed\n ";
				$US_hash{$species} = $tmp_hash{${$US[$us_index]}[4]};
			
				print " +++++ print current US hash INITIAL +++++\n";
				my $key = $species;
				print "KEY: $key ## ";
				my @values = @{$US_hash{$species}};
				foreach (@values) {
					print "$_\t";
				}
				print"\n";

			}
		# or else we need to allow insertions
			else{
				$us_index+=1;
				my $found_us_seq = 0;
				my $us_size = @US-1;
				while (($free_insertions > 0) && ($us_index <= $us_size) && ($found_us_seq==0)){
					$free_insertions-=1;
					print"UPDATE ON FREE INSERTIONS WHILE CHECKING US:  $free_insertions\n"; ## testing
					if (exists $tmp_hash{${$US[$us_index]}[4]}){
						$US_hash{$species} = $tmp_hash{${$US[$us_index]}[4]};
						$found_us_seq = 1;
					
						print " +++++ print current US hash UPDATE +++++\n";
						my $key = $species;
						print "KEY: $key ## ";
						my @values = @{$US_hash{$species}};
						foreach (@values) {
						print "$_\t";
						}
						print"\n";

					}
					$us_index+=1;
				}
				# if no insertions are left anymore, cancel routine and delete species from potential core ortholog set, cause shared synteny can not be established.
				if (not exists $US_hash{$species}){
					if ($free_insertions == 0){
						print "UPSTREAM : No free insertions remaining. Did not find orthologous sequence.\n";
					}	
					elsif ( $us_index == scalar(@US)){
						print "UPSTREAM : No potential orthologs left on upstream chromosome\n";
					}
					print "UPSTREAM : No ortholog found for $species! Skipping species!\n";
					delete $core_gtf_hash{$species};
					delete $core_genome_hash{$species};
				}
			}
##################################################################################################
# check the downstream region of the miRNA 
##################################################################################################
		
			print "initial DS gene of root species: ${$DS[$ds_index]}[4] \n";

			# if the hit is found right away
			if (exists $tmp_hash{${$DS[$ds_index]}[4]}){
				$DS_hash{$species} = $tmp_hash{${$DS[$ds_index]}[4]};
				print " +++++ print current DS hash INITIAL +++++\n";
				my $key = $species;
				print "KEY: $key ## ";
				my @values = @{$DS_hash{$species}};
				foreach (@values) {
					print "$_\t";
				}
				print"\n";
			}

			# or else we need to allow insertions		
			else{
				$ds_index-=1;
				my $found_ds_seq = 0;
				my $ds_size = @DS-1;
				while (($free_insertions > 0) && ($ds_index <= $ds_size) && ($found_ds_seq==0)){
					$free_insertions-=1;
					print"UPDATE ON FREE INSERTIONS WHILE CHECKING DS  $free_insertions\n";
					if (exists $tmp_hash{${$DS[$ds_index]}[4]}){
						$DS_hash{$species} = $tmp_hash{${$DS[$ds_index]}[4]};
						$found_ds_seq = 1;
						print " +++++ print current US hash UPDATE +++++\n";
						my $key = $species;
						print "KEY: $key ## ";
						my @values = @{$DS_hash{$species}};
						foreach (@values) {
							print "$_\t";
						}
						print"\n";
					}
					$ds_index+=1;
				}

				# if no insertions are left anymore, cancel routine and delete species from potential core ortholog set, cause shared synteny can not be established.
				if (not exists $DS_hash{$species}){
					if ($free_insertions == 0){
						print "DOWNSTREAM : No free insertions remaining. Did not find orthologous sequence.\n"; 
					}
					elsif ( $ds_index == scalar(@DS)){
						print "DOWNSTREAM : No potential orthologs left on downstream chromosome\n";
					}
					print "DOWNSTREAM : No ortholog found for $species! Skipping species!\n";
					delete $core_gtf_hash{$species};
					delete $core_genome_hash{$species};
				}
			}
		}

		if (!keys %US_hash){
			print STDERR "Absolutely no ortholog UPSTREAM sequences found! Exiting...\n";
			die;
		}
		elsif (!keys %DS_hash){
			print STDERR "Absolutely no ortholog DOWNSTREAM sequences found! Exiting...\n";
			die;
		}
	
	# extract positions (start,stop,chr) from core gtf file based on ENS ids from the constructed shared syntenic regions	

	########################################################################
	#### TO DO : HANDLE CO ORTHOLOGS ( array size > 1 ) ####################
	
		foreach(keys %core_gtf_hash){
			my $species 	= $_;
			my @US_array 	= @{$US_hash{$species}};
			my @DS_array 	= @{$DS_hash{$species}};
		
			my %tmp_hash = %{$core_gtf_hash{$species}};
			if    ((scalar(@US_array) == 1) && (scalar(@DS_array) == 1)){
				$US_hash{$species} = $tmp_hash{$US_array[0]};
				$DS_hash{$species} = $tmp_hash{$DS_array[0]};
			}
			elsif ((scalar(@US_array) == 1) && (scalar(@DS_array)  > 1)){
	        	$US_hash{$species} = $tmp_hash{$US_array[0]};
				my $US_start= ${$tmp_hash{$US_array[0]}}[0];
				my $highscore_stop=0;
				my $highscore_ensg="";
				foreach(@DS_array){
					my $tmp_ensg = $_;
					my $tmp_stop = ${$tmp_hash{$tmp_ensg}}[1];
					if (($highscore_ensg eq "") && ($tmp_stop < $US_start)){
						$highscore_ensg = $tmp_ensg;
						$highscore_stop = $tmp_stop;
					}
					else{
						if (($tmp_stop > $highscore_stop) && ($tmp_stop < $US_start)){
							$highscore_stop = $tmp_stop;
							$highscore_ensg = $tmp_ensg;
						}	
					}
				}
				$DS_hash{$species} = $tmp_hash{$highscore_ensg};
			}
			elsif ((scalar(@US_array)  > 1) && (scalar(@DS_array) == 1)){
	        	$DS_hash{$species} = $tmp_hash{$DS_array[0]};
	            my $DS_stop= ${$tmp_hash{$DS_array[0]}}[0];
	            my $highscore_start=0;
	            my $highscore_ensg="";
	            foreach(@US_array){
	            	my $tmp_ensg = $_;
	            	my $tmp_start = ${$tmp_hash{$tmp_ensg}}[1];
	            	if (($highscore_ensg eq "") && ($tmp_start > $DS_stop)){
	                 	$highscore_ensg = $tmp_ensg;
	                	$highscore_start = $tmp_start;
	            	}
	            	else{
	            		if (($tmp_start > $highscore_start) && ($tmp_start > $DS_stop)){
	            			$highscore_start = $tmp_start;
	            			$highscore_ensg = $tmp_ensg;
	            		}
	            	}
	            }
	            $US_hash{$species} = $tmp_hash{$highscore_ensg};
			}
			elsif ((scalar(@US_array)  > 1) && (scalar(@DS_array) > 1)){
				my %tmp_us_hash; # {ensg}=$start
				my %tmp_ds_hash; # {ensg}=$stop
				foreach (@US_array){
					$tmp_us_hash{${$tmp_hash{$_}}[4]}=${$tmp_hash{$_}}[0];
				}
				foreach (@DS_array){
					$tmp_ds_hash{${$tmp_hash{$_}}[4]}=${$tmp_hash{$_}}[1];
				}
			}
		}

		### check how many proteins are located between up and downstream
		my %count_hash;
		foreach(keys %core_gtf_hash){
			my $species  = $_;
			my $us_start = ${$US_hash{$_}}[0];
			my $ds_stop  = ${$DS_hash{$_}}[1];
			my $ds_chr   = ${$DS_hash{$_}}[2];
			my $counter = 0;
		
			# check if all values needed to construct the blast library are indeed present.
			# a shared syntenic region might be constructed based on the oma (orthology) data, even though a gene used 
			# in that region might not be annotated in the gtf file. In that case the region has to be skipped because of missing data.
	
			if(not defined $us_start){
				print "Error: For species $species Upstream Start Pos. doesn't exist; skipping species cause of missing GTF annotations!\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
	            delete $core_genome_hash{$species};
				next;
			}
			elsif (not defined $ds_stop){
				print "Error: For species $species Downstream Stop Pos. doesn't exist; skipping species cause of missing GTF annotations!\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
				delete $core_genome_hash{$species};
				next;
			}
			elsif (not defined $ds_chr){
				print "Error: For species $species Chromosome doesn't exist; skipping species cause of missing GTF annotations!\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
				delete $core_genome_hash{$species};
				next;
			}
			else { 
				# also added by mirko
				print "species $species\t start $us_start\t stop $ds_stop \t chr $ds_chr\n";
				my %tmp_hash = %{$core_gtf_hash{$species}};
				foreach(keys %tmp_hash){
					my @tmp_array = @{$tmp_hash{$_}};
					if (($tmp_array[0] > $ds_stop) && ($tmp_array[1] < $us_start)&& ($ds_chr eq $tmp_array[2])){
						$counter+=1;
					}
				}
				$count_hash{$species}=$counter;
			}	
		}
		
		# if there is one core species in the core set, that has more inter proteins than the threshold allows (default 0)
		# delete the species and only take the species that fit the criteria
		foreach(keys %count_hash){
			my $species = $_;
			if ($count_hash{$species}>$max_intra_prot){
				print "$species deleted because of intergenic protein coding sequence\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
				delete $core_genome_hash{$species};
			}
		}
		# count the species that left after deleting the ones over the threshold
		my @carray=keys(%count_hash);
		my $csize=scalar(@carray); # how many species are left	
		my $take='no';
		if (keys %count_hash){
			foreach(keys %count_hash){
				if ($count_hash{$_} <= $max_intra_prot  ){
					$take='yes';
				}
				else{
					$take='no';
				}
			}
		}
		else{
			print STDERR "All core species over threshold! ";
			die "Terminating!\n";
		}
		### parse the intergenic region in the genome of each core species ###
		foreach(keys %core_genome_hash){	
			my $species	= $_;
			my $genome_file = $core_genome_hash{$species};
			my $start	= ${$DS_hash{$_}}[0];
			my $stop	= ${$US_hash{$_}}[1];
			my $chr		= ${$US_hash{$_}}[2];
			my $strand = 2; # set strand to 2, because for intergen region, it is unimportant
			# parse intergenic sequence between up and downstream
			#print "---> start genome parser with: genome:$genome_file\tstart:$start\tstop:$stop\tchr:$chr\tstrand:$strand\n";

			my $inter_seq	= &genome_parser($genome_file,$start,$stop,$chr,$strand);
			#print "$inter_seq\n";		
			my $len_seq = length($inter_seq);
			print "INTERSEQ-LEN: $len_seq\n";
			my $outfile	= $outpath."/".$species.".interseq";	
			open(OUTPUT,">",$outfile);
			print OUTPUT ">$chr\n$inter_seq";
			close(OUTPUT);
			#system("$formatdb -p F -i $outfile");
			system("$formatdb -dbtype nucl -in $outfile");
			#blastn search in intergenic sequence for ncRNA
			my $interseq_blast_out = $outfile.".blastout";
			system("$blastn -task blastn -db $outfile -query $ncRNA -out $interseq_blast_out -outfmt 6 -num_threads $cpu");

			if (-z $interseq_blast_out){
				print STDERR "No BLASTN HIT found , skipping $species!\n";
				next;
			}
			my ($pred_rna_start,$pred_rna_stop,$pred_rna_chr,$pred_rna_name,$pred_rna_strand,$pred_rna_score) = &blast_parser($interseq_blast_out);
			# if the found sequence has less than a $min_seq_len of the the query RNA, it is discarded
			if (($pred_rna_stop-$pred_rna_start)<($ncRNA_seq_len*$min_seq_len)){
				delete $core_genome_hash{$species};
				print STDERR "§§§§ The potential ortholog of $species was too short\n";
				next;
			}
			my $pred_rna = &genome_parser($outfile,$pred_rna_start,$pred_rna_stop,$pred_rna_chr,$pred_rna_strand);
			$rna_file_seq.=">$species\n$pred_rna\n";
		}

	}	
	if ($rna_file_seq eq ""){
		# if this is inactive, t_coffee might run into a core dump, besides... alignment with  1 seq does not make sense...
		print STDERR "§§ The sequence file for the alignment only contains 1 sequence - Exiting $ncRNA!\n";
		die;
	}

	# generates core orthologs outfile 
	$rna_file_seq.=">query\n$ncRNA_seq";  
	open(RNAOUT,">",$rna_file);
	print RNAOUT $rna_file_seq;
	close(RNAOUT);

	print "STEP 01\n";
	system("$tcoffee -cpu $cpu -special_mode=rcoffee -in $rna_file -output=clustalw_aln > $alignment"); 

	print "STEP 02\n";
	system("$tcoffee -other_pg seq_reformat -in $alignment -action +add_alifold -output stockholm_aln > $stockholm_aln");

	print "STEP 03\n";
	system("$cmbuild -F $covariance_model $stockholm_aln");
#	system("$cmbuild $covariance_model $stockholm_aln");
	print "STEP 04\n";
	system("$cmcalibrate --cpu $cpu $covariance_model");
}

else {
	print "Skipping STEPs 01 - 04\n"; 
}

print "STEP 05\n";

system("$cmsearch --cpu $cpu --noali --tblout $cmsearch_out $covariance_model $ukn_genome");

################ parse CMSEARCH output #################
my $cm_chr;
my $cm_start;
my $cm_stop;
my $cm_strand;
my $cm_score=undef;
my $cm_relevant;
my %cm_hash;
my $cm_index = 1; # key for cmsearch output hit 

# hard cut off for cm_search.out file 						## not in use ##
my $cmsearch_out_small = $outpath."/cmsearch_small.out";
system("head -n 52 $cmsearch_out > $cmsearch_out_small");
open(CMOUT,"<",$cmsearch_out_small);

#open(CMOUT,"<",$cmsearch_out);
while(<CMOUT>){
        next if (/^#/);
        my @split_line = split;
        $cm_relevant = $split_line[16];
        if (not defined $cm_score){
                $cm_score = $split_line[14];
                $cm_chr = $split_line[0];
                my $cm_tmp_start = $split_line[7];
                my $cm_tmp_stop  = $split_line[8];
                $cm_strand = $split_line[9];
                if ($cm_strand eq '+'){
                        $cm_strand = 0;
                }
                else{
                        $cm_strand = 1;
                }
                if ($cm_tmp_start < $cm_tmp_stop){
                        $cm_start = $cm_tmp_start;
                        $cm_stop  = $cm_tmp_stop;
                }
                else{
                        $cm_start = $cm_tmp_stop;
                        $cm_stop  = $cm_tmp_start;
                }
                my @tmp_array=($cm_start,$cm_stop,$cm_chr,$cm_strand,$cm_score);
                $cm_hash{$cm_index}=\@tmp_array;
                $cm_index+=1;
        }
        elsif($cm_relevant eq '!'){
                $cm_score = $split_line[14];
                $cm_chr = $split_line[0];
                my $cm_tmp_start = $split_line[7];
                my $cm_tmp_stop  = $split_line[8];
                $cm_strand = $split_line[9];
                if ($cm_strand eq '+'){
                        $cm_strand = 0;
                }
                else{
                        $cm_strand = 1;
                }
                if($cm_tmp_start < $cm_tmp_stop){
                        $cm_start = $cm_tmp_start;
                        $cm_stop  = $cm_tmp_stop;
                }
                else{
                        $cm_start = $cm_tmp_stop;
                        $cm_stop  = $cm_tmp_start;
                }
                my @tmp_array=($cm_start,$cm_stop,$cm_chr,$cm_strand,$cm_score);
                $cm_hash{$cm_index}=\@tmp_array;
                $cm_index+=1;
        }
	#last if $. == $max_line;
}
close(CMOUT);

my %result_hash;	# {index} = @(candidate_name, start, stop, chr, strand, score, RBBH_result);

	############# RECIPROCAL BLAST SEARCH #############

foreach(keys %cm_hash){
	my $cm_index  = $_;
	my @tmp_array 	= @{$cm_hash{$cm_index}};
	my $cmstart 	= $tmp_array[0];
	my $cmstop  	= $tmp_array[1];
	my $cmchr   	= $tmp_array[2];
	my $cmstrand	= $tmp_array[3];
	my $cmscore	= $tmp_array[4];
########
	my $uknRNA = &genome_parser($ukn_genome,$cmstart,$cmstop,$cmchr,$cmstrand);
	my $uknRNA_file = $outpath."/ukn_rna$cm_index.fa";
	open(UKNRNA,">",$uknRNA_file);
	print UKNRNA ">ukn_rna$cm_index\n$uknRNA";
	close(UKNRNA);

	my $reciproc_blast_out = $outpath."/reciproc_blast$cm_index.out";

	system("$blastn -task blastn -db $root_genome -query $uknRNA_file -out $reciproc_blast_out -outfmt 6 -num_threads $cpu");
	
	if (-z $reciproc_blast_out){ 
		print STDERR "No ortholog sequence found in unknown genome ukn_rna$cm_index.fa trying next ukn_rna\n";
		next;
	}
	my ($predict_rna_start,$predict_rna_stop,$predict_rna_chr,$predict_rna_name,$predict_rna_strand,$predict_rna_score) = &blast_parser($reciproc_blast_out); 
	
	############# RECIPROCAL BEST HIT CHECK ##############
	### CHANGED eq, lt, gt to ==, <, > in this version by Andreas ############
	### simplified code ###

	my $rbbh = 0; # equals 1 if rbbh criterion fulfilled, 0 otherwise

	if ($rna_chr eq $predict_rna_chr){ # candidate bbh and mirna gene are on the same contig, now they have to overlap
		if(($predict_rna_start <= $rna_start and $rna_start <= $predict_rna_stop) or ($predict_rna_start <= $rna_stop and $rna_stop <= $predict_rna_stop)){
			$rbbh = 1;
		}
		elsif(($rna_start <= $predict_rna_start and $predict_rna_start <= $rna_stop) or ($rna_start <= $predict_rna_stop and $predict_rna_stop <= $rna_stop)){
			$rbbh = 1;
		}
	}
	if ($rbbh){
		#print "@ Reciprocal Best BLAST Hit was successful! - Orthologous sequence for $query_ncrna_name  at: $cmstart - $cmstop - $cmchr\n";
#		my @result_array = ("ukn_rna$cm_index", $cmstart, $cmstop, $cmchr, $cmstrand, $cmscore, 'SUCCESS');
		my @result_array = ('SUCCESS', "ukn_rna$cm_index", $cmstart, $cmstop, $cmchr, $cmstrand, $cmscore);
		$result_hash{$cm_index}=\@result_array;
	}
	else{ # candidate bbh and mirna coordinates do not match, false candidate
		#print STDERR "RBBH failed for $query_ncrna_name  at: $cmstart - $cmstop - $cmchr!\n";
#		my @result_array = ("ukn_rna$cm_index", $cmstart, $cmstop, $cmchr, $cmstrand, $cmscore, 'FAIL');
		my @result_array = ('FAIL');
		$result_hash{$cm_index}=\@result_array;
	}
}
my $no_hit_found = 0;
foreach (keys %result_hash){
	my $key = $_;
	my @tmp_array = @{$result_hash{$key}};
	if($tmp_array[0] eq 'SUCCESS'){
		print "#### Potential (CO-)Ortholog found: @tmp_array\n";
		$no_hit_found = 1;
		#### mirko edit #####
		# generate an additional outfile if specified by user
		if ($outfile){
			open(R_OUT, '>>',$outpath."/results.out");
			print R_OUT "@tmp_array\n";
			close(R_OUT);
		}
	}
}

if ($no_hit_found == 0){
	print "#### No potential ortholog satisfied the reciprocal best hit criterion.\n";
}
############################################################################################################
######################################### S U B R O U T I N E S ############################################
############################################################################################################
sub genome_parser{
        my $file        		= $_[0];
        my $intergenic_region_start 	= $_[1];
	my $intergenic_region_stop  	= $_[2];
	my $chr         		= $_[3];
	my $strand			= $_[4];
        my $take_seq    = 0;
        my $tmp_chromosome;
        my $tmp_seq;

	#print "#------------# DEBUG GENOME PARSER #------------#\n";
	#print "file: $file\t start:$intergenic_region_start\t stop:$intergenic_region_stop\t chr:$chr\t strand:$strand\n";	

	## construct single genomic sequence from file
	open(GENOME,"<",$file);
        while(<GENOME>){
                if ($_=~/^>/){
                        if ($_=~/^>$chr/){
                                $take_seq = 1;
                        }
                        else{
                                $take_seq = 0;
                        }
                }
                else{
                        if($take_seq == 1){
                                chomp;
                                $tmp_chromosome.=$_;
                        }
                }
        }
        close(GENOME);
	## new parsing routine
#	my $db = Bio::DB::Fasta->new("$file");
#	$tmp_seq = $db->seq($tmp_chromosome,$intergenic_region_start => $intergenic_region_stop);

	## old parsing routine
	if ($intergenic_region_start < $intergenic_region_stop){
		$intergenic_region_start = $intergenic_region_start-1;
        	my $length = $intergenic_region_stop-$intergenic_region_start;	
		$tmp_seq = substr($tmp_chromosome,$intergenic_region_start,$length);
	}
	elsif($intergenic_region_start > $intergenic_region_stop){
		$intergenic_region_stop = $intergenic_region_stop-1;
        	my $length = $intergenic_region_start-$intergenic_region_stop;	
		$tmp_seq = substr($tmp_chromosome,$intergenic_region_stop,$length);
	}
	else{
		print "ERROR in start and stop position!\n";	
	}


	if ($strand == 0){	# only switch Nucleotides
		$tmp_seq =~ tr/tT/uU/;
	}	
	elsif( $strand == 1){	# switch Nucleotides and reverse string
		$tmp_seq = reverse $tmp_seq;
		$tmp_seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;	
		$tmp_seq =~ tr/tT/uU/;
	}
	
        return($tmp_seq);
}

sub blast_parser {
	my $ncPosBLASTout = $_[0];

	my $sub_rna_start;
	my $sub_rna_stop;
	my $sub_rna_chr;
	my $sub_rna_name;

	my $sub_reverse = 0;
        open(BLASTOUT,"<$ncPosBLASTout");
        my $sub_score = undef;


	my $sub_tmp_start;
	my $sub_tmp_stop;

        while(<BLASTOUT>){
                my @sub_blast_line     = split;
                next unless (not defined $sub_score);      
                $sub_rna_name     	= $sub_blast_line[0];
                $sub_tmp_start     	= $sub_blast_line[8];
                $sub_tmp_stop      	= $sub_blast_line[9];
                $sub_rna_chr        	= $sub_blast_line[1];
		$sub_score		= $sub_blast_line[11];
        }
        close(BLASTOUT);

        if($sub_tmp_start < $sub_tmp_stop){
                $sub_rna_start = $sub_tmp_start;
                $sub_rna_stop  = $sub_tmp_stop;
        }
        else{
                $sub_rna_start = $sub_tmp_stop;
                $sub_rna_stop  = $sub_tmp_start;
		$sub_reverse = 1;		# == - strand
        }
	return ($sub_rna_start,$sub_rna_stop,$sub_rna_chr,$sub_rna_name,$sub_reverse,$sub_score);
}
