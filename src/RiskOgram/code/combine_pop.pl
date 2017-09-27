#!/usr/bin/perl

# Rong Chen, SMI
# This program is to combine the LR, snps, snps_all tables on different populations into a single table
# 10/21/2010

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
getopts("h:d:u:p:i:o:", \%options);
unless ((defined $options{h})
	&& (defined $options{d})
	&& (defined $options{u}) 
	&& (defined $options{p})
	&& (defined $options{i})
	&& (defined $options{o})
	){
    &printMessage;
    exit;
}
$host = $options{h};
$db2 = $options{d};
$user = $options{u};
$passwd = $options{p};
$in = $options{i};
$out = $options{o};
# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

@pop = qw(ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI);
$db->do(qq{drop table if exists $out});

for($i=0; $i<@pop; $i++){
    $in_t = $pop[$i].$in;
    if ($i == 0){
	$db->do(qq{create table $out like $in_t});
	$db->do(qq{alter table $out add column pop varchar(16), add index(pop)});
    }
    $db->do(qq{insert into $out select *, "$pop[$i]" from $in_t});
    $db->do(qq{drop table if exists $in_t});
}
$db->do(qq{optimize table $out});	
$db->disconnect();

$end = localtime(time);
print "$out has been successfully created by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to combine the LR, snps, snps_all tables on different populations into a single table
    Available arguments:
    -h mysql host name
    -d database name for output tables
    -u user name
    -p password
    -i input table name on individial population (skipping pop name)
    -o output table name
    
  Usage 1: $0 -h bmir-db1 -d proj_dz_risk_hapmap3_R3 -u rchen1 -p ** -i _AFR_only_T2D_LR -o patient_AFR_only_T2D_LR
  Usage 2: $0 -h bmir-db1 -d proj_dz_risk_hapmap3_R3 -u rchen1 -p ** -i _AFR_only_T2D_snps -o patient_AFR_only_T2D_snps
  Usage 3: $0 -h bmir-db1 -d proj_dz_risk_hapmap3_R3 -u rchen1 -p ** -i _AFR_only_T2D_snps_all -o patient_AFR_only_T2D_snps_all
END
}
