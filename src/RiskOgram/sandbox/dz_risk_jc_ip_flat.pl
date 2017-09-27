#!/usr/bin/perl

# This is a wrapper program that is used to predict disease risk from personal 
# genome sequence
# 
# The original code was built by Rong Chen, Buttelab, Stanford (09/05/2012).
# This is modified by Jieming Chen, Buttelab, UCSF (03/15/17).
# The new modifications allow the use of flat files to replace mySQL tables.
# TODO: currently there is no partial replacement of tables.

#############################################################################
###############################    Header     ###############################
#############################################################################
use warnings; # JC
use strict; # JC
use Getopt::Std;
use DBI;
use English;

#############################################################################
###########################    Main Loop     ###############################
#############################################################################
my $start = localtime(time);
my $current_dir = $ENV{"PWD"};

## Get arguments
my $options = ();
my %options;
getopts("h:d:u:p:i:j:k:o:q:r:g:s:t:", \%options);
unless ((defined $options{h})
	&& (defined $options{d})
	&& (defined $options{u}) 
	&& (defined $options{p})
	&& (defined $options{i})
	&& (defined $options{j})
	&& (defined $options{k})
	&& (defined $options{o})
	&& (defined $options{q})
	&& (defined $options{r})){
    &printMessage;
    exit;
}

my $host = $options{h};
my $db2 = $options{d};
my $user = $options{u};
my $passwd = $options{p};
my $in = $options{i};
my $LR_t = $options{j};
my $R2_t = $options{k}; # table to retrieve pairwise R2 values
my $out_snps = $options{o};
my $out_LR = $options{q};
my $out_snps_all = $options{r};
my $gwas_study_cnt_cutoff;
my $gwas_sample_size_cutoff;

if (defined $options{g}){
    $gwas_study_cnt_cutoff = $options{g};
}
if (defined $options{s}){
    $gwas_sample_size_cutoff = $options{s};
}

## Define variables
my $LR_patient = $out_snps_all . '_tmp';
my $query = '';

## Formulate different LR_patient.pl queries based on optional args -g -s -t
if (defined $options{t}){
		    my $Pmin = $options{t};
		    if ((defined $options{g}) && (defined $options{s})){
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff -s $gwas_sample_size_cutoff -t $Pmin};
		
						print "Ran LR_patient.pl script...1\n"; ##debug
				}
				elsif (defined $options{g}){
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff and -t $Pmin};
		    
		    		print "Ran LR_patient.pl script...2\n"; ##debug
				}
				elsif (defined $options{s}){
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -s $gwas_sample_size_cutoff and -t $Pmin};
		    
		    		print "Ran LR_patient.pl script...3\n"; ##debug
				}
				else{
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -t $Pmin};
		     	
						print "Ran LR_patient.pl script...4\n"; ##debug
				}
}
else {
		    if ((defined $options{g}) && (defined $options{s})){
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff -s $gwas_sample_size_cutoff};
		    }
		    elsif (defined $options{g}){
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff};
		    }
		    elsif (defined $options{s}){
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -s $gwas_sample_size_cutoff};
		    }
		    else{
						$query = qq{LR_patient_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient};
		    }
}

## Run queries
# Run LR_patient.pl to get LR - outputs -o patient_LR
my $time1 = localtime(time);
print "$time1\nRunning LR_patient.pl script...\n"; ##debug
system "$query";
my $time2 = localtime(time);
print "$time2\n";

# Run LR_meta.pl to get disease info (Varimed) - outputs -q dz_risk
#print "Running LR_meta.pl script...\n"; ##debug
#system "LR_meta_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $LR_patient -o $out_snps_all";
#
## Run LR_R2.pl to get and account for LD - outputs -r full_snp_list used after removing those SNPs in LD
#print "Running LR_R2.pl script...\n"; ##debug
#system "LR_R2_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $out_snps_all -j $R2_t -o $out_snps";
#
## Run combine_LR_single.pl to compile the above
#print "Running combine_LR_single.pl script...\n"; ##debug
#system "combine_LR_single_jc_ip.pl -h $host -d $db2 -u $user -p $passwd -i $out_snps -o $out_LR";
#
### Connect to the database
#my $db = DBI->connect("DBI:mysql:database=$db2;host=$host",
#                        $user, $passwd, {'RaiseError' => 1});
#
### Drop tmp table
#$db->do("drop table if exists $LR_patient");
#
### Print to log
#my $end = localtime(time);
#print "$out_LR and $out_snps have been successfully created by $0 from $start to $end.\n";
#exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This wrapper program (sandbox version) is to predict disease risk for any personal genome sequence
    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name for input genome sequence
    -j table name for gender and ethnicity specific LR
    -k table name for ethnicity-specific Linkage Disequilibrium R2
    -o output disease SNP LR table name
    -q output disease risk table name
    -r output table with full SNP list
    -g cutoff for the number of GWAS studies (>=, optional)
    -s cutoff for the total sample size of GWAS studies (>=, optional)
    -t cutoff for Pmin (<=, optional)
    
  Usage: $0 -h buttelab-db1 -d proj_dz_risk -u rchen1 -p * -i p_genotype -j proj_patient_risk.LR_final -k var_ld_data.chr_CEU_all -o p_snps -q p_LR -r p_snps_all -g 2 -s 2000 -t 0.000001
END
}
