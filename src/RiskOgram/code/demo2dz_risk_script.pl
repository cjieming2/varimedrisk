#!/usr/bin/perl

# Rong Chen, SMI
# This program is to create a script to run dz_risk.pl for each patient in the patient_demo table
# 01/25/2011

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
getopts("h:d:u:p:i:j:k:o:q:r:g:s:", \%options);
unless ((defined $options{h})
	&& (defined $options{d})
	&& (defined $options{u}) 
	&& (defined $options{p})
	&& (defined $options{i})
	&& (defined $options{j})
	&& (defined $options{k})
	&& (defined $options{o})
	){
    &printMessage;
    exit;
}
$host = $options{h};
$db2 = $options{d};
$user = $options{u};
$passwd = $options{p};
$demo_t = $options{i};
$LR_t = $options{j};
$ld_t = $options{k};
$out = $options{o};


$tail ="";
if (defined $options{g}){
    $gwas_study_cnt_cutoff = $options{g};
    $tail = $tail." -g $gwas_study_cnt_cutoff";
}
if (defined $options{s}){
    $gwas_sample_size_cutoff = $options{s};
    $tail = $tail." -s $gwas_sample_size_cutoff";
}

if ($demo_t =~ /^(.+)\.(.+)$/){
    $db_demo = $1;
}
else {
    $db_demo = "";
}

# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

if ($ld_t =~ /^(.+chr_)(.+)$/){
    $ld_t_head = $1;
    $ld_t_tail = $2;
}
else {
    print "$LD_t\n";
    print "Error: non-standardized format for LD table\n";
    exit;
}

$query = qq{select distinct p_id, pop_hapmap3 from $demo_t};
$get_patient_id = $db->prepare($query);
$get_patient_id->execute();

open(OUT, ">$out") || die "Error $! when write $out";
print OUT "\#!/bin/tcsh\n\n";

$head = "bsub ~/bin/dz_risk.pl -h $host -d $db2 -u $user -p $passwd ";
while (($p_id, $pop_hapmap3)=$get_patient_id->fetchrow_array()){
    $p_genotype = $db_demo.".".$p_id."_genotype";
    $p_snps = $p_id."_snps";
    $p_LR = $p_id."_LR";
    $p_snps_all = $p_id."_snps_all";
    $p_ld_t = $ld_t_head.$pop_hapmap3."_".$ld_t_tail;
    print OUT "$head -i $p_genotype -j $LR_t -k $p_ld_t -o $p_snps -q $p_LR -r $p_snps_all $tail\n";
}
close(OUT);

$db->disconnect();
$end = localtime(time);
print "$out is successfully updated by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to create a script to run dz_risk.pl for each patient in the patient_demo table.
You will need a p_id and pop_hapmap3 columns in your patient demo table.
The genotype tables should be recored as p_id_genotype table in the same database as your patient demo table.
    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name for patient demographics
    -j table name to retreive LR data
    -k table name to retreive LD R2 data
    -o output script file
    -g cutoff for the number of GWAS studies (>=, optional)
    -s cutoff for the total sample size of GWAS studies (>=, optional)
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_CG_Diversity -u rchen1 -p ** -i var_CG_Diversity.patient_demo -j proj_patient_risk.LR_caucasian -k var_ld_data_LR.chr_all_LR -o run_dz_risk.csh -g 2 -s 2000
END
}
