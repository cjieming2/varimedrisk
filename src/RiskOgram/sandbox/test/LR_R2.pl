#!/usr/bin/perl

# Rong Chen, Stanford
# This program is to select a single SNP with the highestabs(log(LR)) for any pair of SNPs with R_square>=0.8
# 09/02/2009

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
getopts("h:d:u:p:i:j:o:", \%options);
unless ((defined $options{h})
	&& (defined $options{d})
	&& (defined $options{u}) 
	&& (defined $options{p})
	&& (defined $options{i})
	&& (defined $options{j})
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
$t_ld_key = $options{j};
$out = $options{o};

$LD_cutoff = 0.3;
$dist_cutoff = 37000;

if ($t_ld_key =~ /^(.+chr)(_.+)$/){
    $t_ld_head = $1;
    $t_ld_tail = $2;
}
else {
    print "Error: incorrect table name format for LD R_suqare data\n";
    exit(0);
}

%having_ld = ();
@chr = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M);
foreach $chr(@chr){
    $having_ld{$chr}=1;
}

# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});
$db->do("drop table if exists $out");
$db->do("create table $out like $in");
$get_dz_chr = $db->prepare("select distinct broad_phenotype, chromosome from $in order by broad_phenotype, chromosome");
$get_SNP = $db->prepare("select distinct dbSNP, end from $in where broad_phenotype=? and chromosome=? order by gwas_study_cnt desc, study_cnt_total desc, gwas_sample_size desc, study_cnt desc, sample_size desc;");

$get_dz_chr->execute();
while (($dz, $chr)=$get_dz_chr->fetchrow_array()){
    $t_ld = $t_ld_head.$chr.$t_ld_tail;
    $get_LD_SNP = $db->prepare("select distinct dbSNP2 from $t_ld a join $in b ON(a.dbSNP2=b.dbSNP) where a.dbSNP1=? and b.broad_phenotype=? and b.chromosome=? and a.R_square>=$LD_cutoff");
    $get_SNP_dist = $db->prepare("select distinct dbSNP from $in where broad_phenotype=? and chromosome=? and end<? and end>? and dbSNP != ? and dbSNP not in (select dbSNP2 from $t_ld where dbSNP1=?)");
    @dbSNP = ();
    %end = ();
    $get_SNP->execute($dz, $chr);
    while(($dbSNP, $end)=$get_SNP->fetchrow_array()){
	push(@dbSNP, $dbSNP);
	$end{$dbSNP} = $end;
    }
    $get_SNP->finish();
    
    if (exists $having_ld{$chr}){
	for($i=0; $i<@dbSNP; $i++){
	    %dbSNP_LD = ();
	    $get_LD_SNP->execute($dbSNP[$i], $dz, $chr);
	    while ($dbSNP1=$get_LD_SNP->fetchrow_array()){
		$dbSNP_LD{$dbSNP1}=1; 
	    }
	    $get_LD_SNP->finish();
	    $get_SNP_dist->execute($dz, $chr, $end{$dbSNP[$i]}+$dist_cutoff, $end{$dbSNP[$i]}-$dist_cutoff, $dbSNP[$i], $dbSNP[$i]);
	    while ($dbSNP1=$get_SNP_dist->fetchrow_array()){
                $dbSNP_LD{$dbSNP1}=2;
            }
            $get_SNP_dist->finish();
	    for($j=$i+1; $j<@dbSNP; $j++){
		if (exists $dbSNP_LD{$dbSNP[$j]}){
		    print "Remove rs$dbSNP[$j] on chromome $chr on $dz due to LD with rs$dbSNP[$i] type: $dbSNP_LD{$dbSNP[$j]}\n";
		    splice(@dbSNP, $j, 1);
		    $j--;
		}
	    }
	}
    }
    #print "@dbSNP\n";
    $SNP_list = "(".join(", ", @dbSNP).")";
    #print "$SNP_list\n";
    $insert_data = $db->prepare("insert into $out select * from $in where broad_phenotype=? and chromosome=? and dbSNP in $SNP_list");
    $insert_data->execute($dz, $chr);
    $insert_data->finish();
    #print "insert data broad_phenotype=$dz and chromosome=$chr and dbSNP in ($SNP_list)\n";
}
$get_dz_chr->finish();

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
This program is to keep a single SNP with the highest abs(log(LR)) from any SNP pair with R_square>=0.3 on any disease/chromosome pair.

    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name
    -j table to retrieve the R suqare value of any SNP pair
    -o output table name
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_westome -u rchen1 -p ** -i 20_LR -j var_ld_data.chr_CEU_all -o 20_LR_uniq
END
}
