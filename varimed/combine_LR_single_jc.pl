#!/usr/bin/perl

# This program is to report final risk for a patient by accumulating all individual SNPs 
#
# The original code was built by Rong Chen, Buttelab, Stanford (08/19/2009).
# This is modified by Jieming Chen, Buttelab, UCSF (03/15/17).
# TODO: it doesnt seem to deal with the pre-test probability
# 1) using only 1 LR if there are multiple same dbSNP IDs due to overlapping with multiple genes -- to prevent redundancy
# 2) this is assuming that distinct(dbSNP, LR) is valid, that is, this is only removing exact redundacies
# --> to-do: choose highest LR if same rsID

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

# get arguments
my $options = ();
my %options;
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
my $host = $options{h};
my $db2 = $options{d};
my $user = $options{u};
my $passwd = $options{p};
my $in = $options{i};
my $out = $options{o};

my $out_file = $current_dir."/".$out.".txt";
# Connect to the database.
my $db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

my $query = qq{drop table if exists $out};
$db->do($query);
$query = qq{create table $out(
			      broad_phenotype varchar(128), index(broad_phenotype),
			      LR double, index(LR),
			      LR_max float, index(LR_max),
                              SNP_risk int unsigned, index(SNP_risk),
			      SNP_protective int unsigned, index(SNP_protective)
                              )
	    };
$db->do($query);

$query = qq{select distinct broad_phenotype from $in};
my $get_phenotype = $db->prepare($query);

$query = qq{select count(distinct dbSNP) from $in where broad_phenotype=? and LR>=1};
my $get_SNP_risk = $db->prepare($query);

$query = qq{select count(distinct dbSNP) from $in where broad_phenotype=? and LR<1};
my $get_SNP_protective = $db->prepare($query);

## produces _dz_risk table with input from LR_R2.pl
## LR is exp(sum of all logLRs) of all SNPs
#$query = qq{select exp(sum(log(LR))) from $in where broad_phenotype=?}; # JC
$query = qq{select exp(sum(log(LR))) from (select distinct dbSNP, LR from $in where broad_phenotype=?) AS new}; # JC; derived table requires alias, so 'new' is a dummy alias
my $get_LR = $db->prepare($query);

## LR_max picks the top SNP with max absolute logLR
$query = qq{select LR from $in where broad_phenotype=? order by abs(log(LR)) desc limit 1};
my $get_LR_max = $db->prepare($query);

open(OUT, ">$out_file") || die "Error: $! when write $out_file";
$get_phenotype->execute();

my $SNP_risk;
my $SNP_protective;
my $LR;
my $LR_max;
while ((my $phenotype) = $get_phenotype->fetchrow_array()){
    $get_SNP_risk->execute($phenotype);
    $SNP_risk = $get_SNP_risk->fetchrow_array();
    $get_SNP_risk->finish();
    
    $get_SNP_protective->execute($phenotype);
    $SNP_protective = $get_SNP_protective->fetchrow_array();
    $get_SNP_protective->finish();

    $get_LR->execute($phenotype);
    $LR = $get_LR->fetchrow_array();
    $get_LR->finish();

    $get_LR_max->execute($phenotype);
    $LR_max = $get_LR_max->fetchrow_array();
    $get_LR_max->finish();
    
    printf OUT "%s\t%.12f\t%.6f\t%d\t%d\n", $phenotype, $LR, $LR_max, $SNP_risk, $SNP_protective;
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

my $end = localtime(time);
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
This program is to report final LR for a single patient by accumulating all individual SNPs.
There is no patient_id column.

    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name
    -o output table name
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_ref -u rchen1 -p ** -i ref_LR_uniq -o ref_LR_final
END
}
