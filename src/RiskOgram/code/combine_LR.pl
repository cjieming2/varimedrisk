#!/usr/bin/perl

# Rong Chen, SMI
# This program is to report final risk for a patient by accumulating all individual SNPs and times it with the pre-test probability
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

$out_file = $current_dir."/".$out.".txt";
# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

$query = qq{drop table if exists $out};
$db->do($query);
$query = qq{create table $out(patient_id int unsigned, index(patient_id),
			      broad_phenotype varchar(128), index(broad_phenotype),
			      LR float(7,4), index(LR),
			      LR_max float(7,4), index(LR_max),
                              SNP_risk int unsigned, index(SNP_risk),
			      SNP_protective int unsigned, index(SNP_protective)
                              )
	    };
$db->do($query);

$query = qq{select distinct patient_id, broad_phenotype from $in};
$get_phenotype = $db->prepare($query);

$query = qq{select count(distinct dbSNP) from $in where patient_id=? and broad_phenotype=? and LR>=1};
$get_SNP_risk = $db->prepare($query);

$query = qq{select count(distinct dbSNP) from $in where patient_id=? and broad_phenotype=? and LR<1};
$get_SNP_protective = $db->prepare($query);

$query = qq{select exp(sum(log(LR))) from $in where patient_id=? and broad_phenotype=?};
$get_LR = $db->prepare($query);

$query = qq{select LR from $in where patient_id=? and broad_phenotype=? order by abs(log(LR)) desc limit 1};
$get_LR_max = $db->prepare($query);

open(OUT, ">$out_file") || die "Error: $! when write $out_file";
$get_phenotype->execute();
while (($patient_id, $phenotype) = $get_phenotype->fetchrow_array()){
    $get_SNP_risk->execute($patient_id, $phenotype);
    $SNP_risk = $get_SNP_risk->fetchrow_array();
    $get_SNP_risk->finish();
    
    $get_SNP_protective->execute($patient_id, $phenotype);
    $SNP_protective = $get_SNP_protective->fetchrow_array();
    $get_SNP_protective->finish();

    $get_LR->execute($patient_id, $phenotype);
    $LR = $get_LR->fetchrow_array();
    $get_LR->finish();

    $get_LR_max->execute($patient_id, $phenotype);
    $LR_max = $get_LR_max->fetchrow_array();
    $get_LR_max->finish();
    
    printf OUT "%d\t%s\t%.6f\t%.6f\t%d\t%d\n", $patient_id, $phenotype, $LR, $LR_max, $SNP_risk, $SNP_protective;
}
$get_phenotype->finish();
close(OUT);
$query = qq{load data local infile \'$out_file\' into table $out};
$db->do($query);
unlink($out_file);
$query = qq{alter table $out enable keys};
$db->do($query);
$query = qq{optimize table $out};
$db->do($query);

$end = localtime(time);
print "$out is successfully updated by $0 from $start to $end.\n";
exit(0);

#############################################################################
# Flip genotype for backward strand
#############################################################################
sub flip()
{
    my ($data)=@_;
    my @token2 = ();
    my @token = split(//, $data);
    for (my $i=0; $i<@token; $i++){
	if ($token[$i] eq 'A'){
	    $token2[$i] = 'T';
	}
	elsif ($token[$i] eq 'C'){
            $token2[$i] = 'G';
        }
	elsif ($token[$i] eq 'G'){
            $token2[$i] = 'C';
        }
	elsif ($token[$i] eq 'T'){
            $token2[$i] = 'A';
        }
	else {
	    print "Warning: can not flip strand for $data\n";
	}
    }
    return join ("", @token2);
}
#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to report final LR for every patient by accumulating all individual SNPs

    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name
    -o output table name
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_1000genome -u rchen1 -p ** -i patient_dbSNP_LR -o patient_LR
END
}
