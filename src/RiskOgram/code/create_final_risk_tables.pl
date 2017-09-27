#!/usr/bin/perl

# Rong Chen, SMI
# This program is to load LD data from LD data into mysql
# 08/19/2009

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
getopts("h:d:u:p:i:o:q:", \%options);
unless ((defined $options{h})
	&& (defined $options{d})
	&& (defined $options{u}) 
	&& (defined $options{p})
	&& (defined $options{i})
	&& (defined $options{o})
        && (defined $options{q})
	){
    &printMessage;
    exit;
}
$host = $options{h};
$db2 = $options{d};
$user = $options{u};
$passwd = $options{p};
$patient_t = $options{i};
$out_risk_t = $options{o};
$out_LR_t = $options{q};
# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

$query = qq{drop table if exists $out_risk_t};
$db->do($query);
$query = qq{create table $out_risk_t(patient_id int unsigned, index(patient_id),
				 broad_phenotype varchar(128), index(broad_phenotype),
				 risk_post float(10,7),
				 risk_pre float(10,7),
				 LR float(7,4), index(LR),
				 SNP_risk int unsigned, index(SNP_risk),
				 SNP_protective int unsigned, index(SNP_Protective))};
$db->do($query);
$query = qq{alter table $out_risk_t disable keys};
$db->do($query);
$query = qq{drop table if exists $out_LR_t};
$db->do($query);
$query = qq{create table $out_LR_t(patient_id int unsigned, index(patient_id),
				 broad_phenotype varchar(128), index(broad_phenotype),
				 dbSNP bigint(20), index(dbSNP),
				 genotype varchar(64), index(genotype),
				 LR float, index(LR),
				 study_cnt int(11), index(study_cnt),
				 sample_size int(11), index(sample_size),
				 chromosome varchar(8), index(chromosome),
				 end int(10) unsigned, index(end),
				 fxn_type varchar(64), index(fxn_type),
				 hapBlock varchar(64), index(hapBlock),
				 symbol varchar(64), index(symbol),
				 study_cnt_total int(10), index(study_cnt_total),
				 gwas varchar(8), index(gwas))};
$db->do($query);
$query = qq{alter table $out_LR_t disable keys};
$db->do($query);

$query = qq{select distinct id from $patient_t where pilot1='Y' or pilot2='Y' order by id};
$get_patient_id = $db->prepare($query);
$get_patient_id->execute();
while ($patient_id = $get_patient_id->fetchrow_array()){
    $in_risk = $patient_id."_risk";
    $in_LR = $patient_id."_LR";
    $query = qq{insert into $out_risk_t select '$patient_id', broad_phenotype, risk_post, risk_pre, LR, SNP_risk, SNP_protective from $in_risk};
    $db->do($query);
    $query = qq{insert into $out_LR_t select '$patient_id', broad_phenotype, dbSNP, genotype, LR,
		study_cnt, sample_size, chromosome, end, fxn_type, hapBlock,
		symbol, study_cnt_total, gwas from $in_LR};
    $db->do($query);
}
$get_patient_id->finish();
$query = qq{alter table $out_risk_t enable keys};
$db->do($query);
$query = qq{optimize table $out_risk_t};
$db->do($query);

$query = qq{alter table $out_LR_t enable keys};
$db->do($query);
$query = qq{optimize table $out_LR_t};
$db->do($query);

$end = localtime(time);
print "$out_risk_t and $out_LR_t are successfully updated by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to combine personal risk and LR tables into the final risk and dbSNP_LR table.
    Available arguments:
    -h mysql host name
    -d database name for output tables
    -u user name
    -p password
    -i table name for patient information
    -o output table name for disease risk
    -q output table name for dbSNP_LR
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_1000genome -u rchen1 -p ** -i var_1000genome.patient_demo -o patient_risk -q patient_dbSNP_LR
END
}
