#!/usr/bin/perl

# Rong Chen, SMI
# This program is to run dz_risk.pl on a list of patients
# 12/17/2010

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
$genotype_t = $options{i};
$LR_t = $options{j};
$ld_t = $options{k};
$out_snps = $options{o};
$out_LR = $options{q};
$out_snps_all = $options{r};

if (defined $options{g}){
    $gwas_study_cnt_cutoff = $options{g};
}
if (defined $options{s}){
    $gwas_sample_size_cutoff = $options{s};
}

# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

$query = qq{select distinct patient_id from $genotype_t};
$get_patient_id = $db->prepare($query);
$get_patient_id->execute();

# get list of dbSNPs in LR table
$db->do(qq{create temporary table LR(dbSNP bigint unsigned primary key)});
$db->do(qq{insert into LR select distinct dbSNP from $LR_t});
$db->do(qq{optimize table LR});

# get a smaller genotype data table $genotype_t2 with LR data
$genotype_t2 = $genotype_t ."_LR_all";
$db->do(qq{drop table if exists $genotype_t2});
$db->do(qq{create table $genotype_t2 like $genotype_t});
$db->do(qq{alter table $genotype_t2 disable keys});
$db->do(qq{insert into $genotype_t2 select a.* from $genotype_t a join LR b using(dbSNP)});
$db->do(qq{alter table $genotype_t2 enable keys});
$db->do(qq{optimize table $genotype_t2});

$p_cnt=0;
while (($patient_id) = $get_patient_id->fetchrow_array()){
    $p_genotype = $patient_id."_genotype";
    $db->do(qq{drop table if exists $p_genotype});
    $db->do(qq{create table $p_genotype like $genotype_t2});
    $db->do(qq{alter table $p_genotype disable keys});
    $db->do(qq{insert into $p_genotype select * from $genotype_t2 where patient_id=$patient_id});
    $db->do(qq{alter table $p_genotype enable keys});
    $db->do(qq{optimize table $p_genotype});

    $p_snps = $patient_id."_snps";
    $p_LR = $patient_id."_LR";
    $p_snps_all = $patient_id."_snps_all";

    if ((defined $options{g}) && (defined $options{s})){
	$query = qq{dz_risk.pl -h $host -d $db2 -u $user -p $passwd -i $p_genotype -j $LR_t -k $ld_t -o $p_snps -q $p_LR -r $p_snps_all -g $gwas_study_cnt_cutoff -s $gwas_sample_size_cutoff};
    }
    elsif (defined $options{g}){
	$query = qq{dz_risk.pl -h $host -d $db2 -u $user -p $passwd -i $p_genotype -j $LR_t -k $ld_t -o $p_snps -q $p_LR -r $p_snps_all -g $gwas_study_cnt_cutoff};
    }
    elsif (defined $options{s}){
	$query = qq{dz_risk.pl -h $host -d $db2 -u $user -p $passwd -i $p_genotype -j $LR_t -k $ld_t -o $p_snps -q $p_LR -r $p_snps_all -s $gwas_sample_size_cutoff};
    }
    else{
	$query = qq{dz_risk.pl -h $host -d $db2 -u $user -p $passwd -i $p_genotype -j $LR_t -k $ld_t -o $p_snps -q $p_LR -r $p_snps_all};
    }
    system "$query";
    if ($p_cnt == 0){
	$db->do(qq{drop table if exists $out_snps});
	$db->do(qq{create table $out_snps like $p_snps});
	$db->do(qq{alter table $out_snps add column patient_id int unsigned, add index(patient_id)});
	$db->do(qq{alter table $out_snps disable keys});
	
	$db->do(qq{drop table if exists $out_LR});
	$db->do(qq{create table $out_LR like $p_LR});
	$db->do(qq{alter table $out_LR add column patient_id int unsigned, add index(patient_id)});
	$db->do(qq{alter table $out_LR disable keys});

	$db->do(qq{drop table if exists $out_snps_all});
	$db->do(qq{create table $out_snps_all like $p_snps_all});
	$db->do(qq{alter table $out_snps_all add column patient_id int unsigned, add index(patient_id)});
	$db->do(qq{alter table $out_snps_all disable keys});
    }
    $db->do(qq{insert into $out_snps select *, $patient_id from $p_snps});
    $db->do(qq{insert into $out_LR select *, $patient_id from $p_LR});
    $db->do(qq{insert into $out_snps_all select *, $patient_id from $p_snps_all});

    $db->do(qq{drop table if exists $p_genotype});
    $db->do(qq{drop table if exists $p_snps});
    $db->do(qq{drop table if exists $p_LR});
    $db->do(qq{drop table if exists $p_snps_all});
    
    $p_cnt++;
}
$get_patient_id->finish();

$db->do(qq{drop table if exists $genotype_t2});
$db->do(qq{alter table $out_snps enable keys});
$db->do(qq{optimize table $out_snps});
$db->do(qq{alter table $out_LR enable keys});
$db->do(qq{optimize table $out_LR});
$db->do(qq{alter table $out_snps_all enable keys});
$db->do(qq{optimize table $out_snps_all});

$db->disconnect();
$end = localtime(time);
print "$out_LR and $out_snps are successfully updated by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to run dz_risk.pl and analyze the genetic risk for a patient genotype table.
    Available arguments:
    -h mysql host name
    -d database name for output tables
    -u user name
    -p password
    -i input table name for patient genotypes
    -j table name to retreive LR data
    -k table name to retreive LD R2 data
    -o output table for snps
    -q output table for LR
    -r output table for the full list of snps
    -g cutoff for the number of GWAS studies (>=, optional)
    -s cutoff for the total sample size of GWAS studies (>=, optional)
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_hapmap3 -u rchen1 -p ** -i var_hapmap3.CEU_west_genotype -j proj_patient_risk.LR_caucasian -k var_ld_data_LR.chr_CEU_all_LR -o CEU_west_snps -q CEU_west_LR -r CEU_west_snps_all -g 2 -s 2000
END
}
