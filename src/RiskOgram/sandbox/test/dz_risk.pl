#!/usr/bin/perl

# Rong Chen, Buttelab, Stanford
# This program is to predict disease risk from personal genome sequence
# 09/05/2012

#############################################################################
###############################    Header     ###############################
#############################################################################
use Getopt::Std;
use DBI;
use English;

#############################################################################
###########################    Main Loop     ###############################
#############################################################################
$start = localtime(time);
$current_dir = $ENV{"PWD"};

# get arguments
$options = ();
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
	&& (defined $options{r})
	){
    &printMessage;
    exit;
}
$host = $options{h};
$db2 = $options{d};
$user = $options{u};
$passwd = $options{p};
$in = $options{i};
$LR_t = $options{j};
$R2_t = $options{k}; # table to retrieve pairwise R2 values
$out_snps = $options{o};
$out_LR = $options{q};
$out_snps_all = $options{r};

if (defined $options{g}){
    $gwas_study_cnt_cutoff = $options{g};
}
if (defined $options{s}){
    $gwas_sample_size_cutoff = $options{s};
}

$LR_patient = $out_snps_all."_tmp";

if (defined $options{t}){
    $Pmin = $options{t};
    if ((defined $options{g}) && (defined $options{s})){
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff -s $gwas_sample_size_cutoff -t $Pmin};

	print "Ran LR_patient.pl script...1\n"; ##debug
    }
    elsif (defined $options{g}){
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff and -t $Pmin};
    	print "Ran LR_patient.pl script...2\n"; ##debug
    }
    elsif (defined $options{s}){
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -s $gwas_sample_size_cutoff and -t $Pmin};
    	print "Ran LR_patient.pl script...3\n"; ##debug
    }
    else{
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -t $Pmin};
     	
	print "Ran LR_patient.pl script...4\n"; ##debug
        }
}
else {
    if ((defined $options{g}) && (defined $options{s})){
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff -s $gwas_sample_size_cutoff};
    }
    elsif (defined $options{g}){
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -g $gwas_study_cnt_cutoff};
    }
    elsif (defined $options{s}){
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient -s $gwas_sample_size_cutoff};
    }
    else{
	$query = qq{LR_patient.pl -h $host -d $db2 -u $user -p $passwd -i $LR_t -j $in -o $LR_patient};
    }
}

system "$query";
system "LR_meta.pl -h $host -d $db2 -u $user -p $passwd -i $LR_patient -o $out_snps_all";
system "LR_R2.pl -h $host -d $db2 -u $user -p $passwd -i $out_snps_all -j $R2_t -o $out_snps";
system "combine_LR_single.pl -h $host -d $db2 -u $user -p $passwd -i $out_snps -o $out_LR";

# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

$db->do("drop table if exists $LR_patient");

$end = localtime(time);
print "$out_LR and $out_snps have been successfully created by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to predict disease risk for any personal genome sequence
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
