#!/usr/bin/perl

# Rong Chen, Dept. of Pediatrics, Stanford University
# This program is to combine LR, snps, snps_all tables from individuals into LR, snps, snps_all tables
# 01/26/2011

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
getopts("h:d:u:p:i:o:q:r:", \%options);
unless ((defined $options{h})
	&& (defined $options{d})
	&& (defined $options{u}) 
	&& (defined $options{p})
	&& (defined $options{i})
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
$in_demo = $options{i};
$out_LR = $options{o};
$out_snps = $options{q};
$out_snps_all = $options{r};

# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

$query = qq{select distinct p_id, pop from $in_demo order by p_id};
$get_p_id = $db->prepare($query);
$get_p_id->execute();

$i=0;
while (($p_id, $pop) = $get_p_id->fetchrow_array()){
    $p_LR = $p_id."_LR";
    $p_snps = $p_id."_snps";
    $p_snps_all = $p_id."_snps_all";

    if ($i == 0){
	$db->do(qq{drop table if exists $out_LR});
	$db->do(qq{drop table if exists $out_snps});
	$db->do(qq{drop table if exists $out_snps_all});
	$db->do(qq{create table $out_LR like $p_LR});
	$db->do(qq{create table $out_snps like $p_snps});
	$db->do(qq{create table $out_snps_all like $p_snps_all});
	$db->do(qq{alter table $out_LR add column p_id varchar(64), add column pop varchar(64), add index(p_id), add index(pop)});
	$db->do(qq{alter table $out_snps add column p_id varchar(64), add column pop varchar(64), add index(p_id), add index(pop)});
	$db->do(qq{alter table $out_snps_all add column p_id varchar(64), add column pop varchar(64), add index(p_id), add index(pop)});
    }
    $db->do(qq{insert into $out_LR select *, \"$p_id\", \"$pop\" from $p_LR});
    $db->do(qq{insert into $out_snps select *, \"$p_id\", \"$pop\" from $p_snps});
    $db->do(qq{insert into $out_snps_all select *, \"$p_id\", \"$pop\" from $p_snps_all});

    # remove individual tables
    $db->do(qq{drop table if exists $p_LR});
    $db->do(qq{drop table if exists $p_snps});
    $db->do(qq{drop table if exists $p_snps_all});

    $i++;
}
$get_p_id->finish();

$db->do(qq{optimize table $out_LR});
$db->do(qq{optimize table $out_snps});
$db->do(qq{optimize table $out_snps_all});

$end = localtime(time);
print "$out_LR, $out_snps, and $out_snps_all are successfully updated by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to combine LR, snps, snps_all tables from individuals into LR, snps, snps_all tables
    Available arguments:
    -h mysql host name
    -d database name for output tables
    -u user name
    -p password
    -i input table name for patient demo
    -o output table name for LR
    -q output table name for snps
    -r output table name for snp_all
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_CG_Diversity -u rchen1 -p ** -i var_CG_Diversity.patient_demo -o patient_LR -q patient_snps -r patient_snps_all
END
}
