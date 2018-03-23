#! /usr/bin/perl
# Author : Daniel Amsel - daniel.amsel@gmx.de
use strict;
use warnings;
use Getopt::Long;
use Storable;

##################################################################################
############################# PATHVARIABLES TO EDIT###############################
##################################################################################

my $formatdb = 'formatdb';		# path or command to formatdb	

##################################################################################
######################## DO NOT TOUCH THE CODE BELOW THIS LINE####################
##################################################################################


#######################
##### VARIABLES #######
#######################
my $root_genome;     ##
my $root_gtf;        ##
#######################
my $core_genome_path;##
my $core_gtf_path;   ##
#######################
my $oma_ortho_path;  ##
#######################
my $help;            ##
if (@ARGV==0) {      ##
        $help = 1;   ##
}                    ##
#######################

my $helpmessage = 	"
===============================================================================
| 			PRECOMPUTING SCRIPT FOR NCRNA_FINDER	 	      |
|=============================================================================|
|									      |	
|	-root_genome <> 		=	path/to/root/genome/.fa	      |		
|	-root_gtf <>			=	path/to/root/gtf/.gtf  	      |
|									      |	
|	-core_genome_folder <>		=	path/to/core/genome	      |
|	-core_gtf_folder <>		=	path/to/core/gtf    	      |
|									      |	
|	-oma_ortho_folder <>		=	path/to/oma/		      |
|	-h|help				=	prints this help              |
|									      |
| !!!!!	Please note to provide the full paths to the algorithm !!!!!	      |
|=============================================================================|
===============================================================================\n";


##################################################################################
GetOptions      (       'root_genome=s'         =>      \$root_genome           ,
                        'root_gtf=s'            =>      \$root_gtf              ,
                        'core_genome_folder=s'  =>     	\$core_genome_path      ,
                        'core_gtf_folder=s'     =>      \$core_gtf_path         ,
			'oma_ortho_folder=s'	=>	\$oma_ortho_path	,
			'h|help'		=>	\$help			);
##################################################################################


###########################
if ($help) {
        print $helpmessage;
        exit;
}
###########################


#######################################################################
my $root_gtf_hash_file          = $root_gtf.".hash";
my $core_gtf_hash_file          = $core_gtf_path."core_gtf.hash";
my $core_genome_hash_file	= $core_genome_path."core_genome.hash";
my $oma_ortho_hash_file		= $oma_ortho_path."oma_ortho.hash";
#######################################################################
my %root_gtf_hash;
my %core_gtf_hash;
my %core_genome_hash;
my %oma_ortho_hash;
##################

#####################################################
############### make BLAST DBs ######################
#####################################################
system("$formatdb -p F -i $root_genome");
print "finished root genome blast db\n";				
#####################################################
############## SAVE ROOT GTF TO HASH ################
 #%root_gtf_hash{ENSG} = @(start,stop,chr,enst,ensg)#
#####################################################
if (-e $root_gtf_hash_file){						
        unlink $root_gtf_hash_file;
}

%root_gtf_hash = %{&gtf_parser($root_gtf)};
store \%root_gtf_hash,$root_gtf_hash_file;
print "finished root gtf hashing\n";					
################################################
########## SAVE CORE GTF TO HASH ###############				
################################################
if (-e $core_gtf_hash_file){
        unlink $core_gtf_hash_file;
}
opendir DIR, $core_gtf_path;
my @core_gtf_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
closedir DIR;

foreach (@core_gtf_files){
	my $core_gtf_file 	= $_;
     	my $species_name    	= &directory_parser($core_gtf_file);
        $core_gtf_hash{$species_name}       	= &gtf_parser($core_gtf_path.$core_gtf_file); 
}
store \%core_gtf_hash,$core_gtf_hash_file;
print "finished core gtf hashing\n";				# testprint
################################################
######### SAVE CORE GENOME TO HASH #############				
################################################
if (-e $core_genome_hash_file){
        unlink $core_genome_hash_file;
}
opendir DIR, $core_genome_path;
my @core_genome_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
closedir DIR;
foreach(@core_genome_files){
	my $species_name = &directory_parser($_);
	$core_genome_hash{$species_name}	= $core_genome_path.$_;			
}

store \%core_genome_hash,$core_genome_hash_file;

print "finished core genome hashing\n";				# testprint
################################################
########## SAVE OMA ORTHOLOGS TO HASH ##########					
################################################
if (-e $oma_ortho_hash_file){
        unlink $oma_ortho_hash_file;
}

opendir DIR, $oma_ortho_path;
my @oma_ortho_files = grep { $_ ne '.' && $_ ne '..' && $_ ne '*~'} readdir DIR;
closedir DIR;


foreach(@oma_ortho_files){
	my $oma_ortho_file	= $_; #N.leucogenys
	print "$oma_ortho_file\n";
	my $oma_file_path = $oma_ortho_path."/".$oma_ortho_file;
	my %tmp_hash = %{&oma_parser($oma_file_path)};

#	my %tmp_hash = %{&oma_parser($oma_ortho_path.$oma_ortho_file)};
	$oma_ortho_hash{$oma_ortho_file} = \%tmp_hash;
}

store \%oma_ortho_hash,$oma_ortho_hash_file;
print "finished oma hashing\n";					# testprint



print "Pre-computation finished!\n";

#########################################################################################
############################### S U B R O U T I N E S ###################################
#########################################################################################
sub directory_parser{  	
	my $core_gtf_file       = $_[0];
        my @split_file_name     = split('\.',$core_gtf_file);  #Nomascus_leucogenys.Nleu1.0.76.gtf
        my $species_name        = $split_file_name[0];         #Nomascus_leucogenys
        print "Initial species name: $species_name\n";
        $species_name           =~ s/_/ /g;                    #Nomascus leucogenys
        my @species_name_split  = split(' ',$species_name);
        $species_name    = substr($species_name, 0, 1).".".$species_name_split[1]; #N.leucogenys
        print "Parsed species name: $species_name\n";
	return $species_name;
}


###############
sub oma_parser{
	my %tmp_hash;
	my $oma_file	= $_[0];
	open(OMA,"<",$oma_file);
	while(<OMA>){
		my @split_line = split;
		if (not exists $tmp_hash{$split_line[0]}){			
			my @tmp_array = ($split_line[1]);
			$tmp_hash{$split_line[0]}=\@tmp_array;
		}
		else{
			my @tmp_array = @{$tmp_hash{$split_line[0]}};
			push (@tmp_array,$split_line[1]);
			$tmp_hash{$split_line[0]}=\@tmp_array;
		}
	}
	close(OMA);
	return (\%tmp_hash);
}
###############  changed by mirko
sub gtf_parser{
        my %tmp_ensg_hash;
        my $gtf_file = $_[0];
        my %gtf_hash;
        open(GTF,"<",$gtf_file) || die "can not open file : $!\n";
        while(<GTF>){                
		my $line=$_;
		if ($line =~ /^#/){
			next;
		}
                my @split_line = split ("\t", $line);				
		
                next unless ($split_line[2] =~ /^transcript/);	
		
		my @gtf_attributes = split ("; ", $split_line[8]);		
		my %tag_value_pairs;
		foreach (@gtf_attributes){
			my @tag_value = split (" ",$_);
			my $tag = $tag_value[0];
			my $value = $tag_value[1];
			$value =~ s/"//g;
			$tag_value_pairs{$tag} = $value;				
		}
		next unless $tag_value_pairs{'gene_biotype'} eq "protein_coding"; # geht nicht bei core files
		my $ENSG = $tag_value_pairs{'gene_id'};
		my $ENST = $tag_value_pairs{'transcript_id'};

                my $chr   = $split_line[0];
                my $start = $split_line[3];
                my $stop  = $split_line[4];

                if (not exists $tmp_ensg_hash{$ENSG}){
                        my @tmp_array = ($start,$stop,$chr,$ENST,$ENSG);
                        $tmp_ensg_hash{$ENSG}=\@tmp_array;
                }
                else{
                        my $saved_start = ${$tmp_ensg_hash{$ENSG}}[0];
                        my $saved_stop  = ${$tmp_ensg_hash{$ENSG}}[1];
                        my $saved_len   = $saved_stop - $saved_start + 1;		# take only the longest available transcript into account
                        my $new_len     = $stop - $start + 1;
                        if ($new_len > $saved_len){
                                my @rewrite_array       = ($start,$stop,$chr,$ENST,$ENSG);
                                $tmp_ensg_hash{$ENSG}   = \@rewrite_array;
                        }
                }
        }
	return \%tmp_ensg_hash;
}




