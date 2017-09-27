#!/usr/bin/perl

# This program is to select a single SNP 
# with the largest (gwas_study_size, total_cnt_gwas, gwas_sample_size, study_cnt, sample_size) 
# for any pair of SNPs with R_square>=0.5, or distance of 40kb
#
# The original code was built by Rong Chen, Buttelab, Stanford (08/19/2009).
# This is modified by Jieming Chen, Buttelab, UCSF (03/15/17), by
# 1) including comments to make code more readable
# 2) adding "use warnings" and "use strict"
# 3) declaring variables, and adding try and if statements to catch exceptions


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
my $timestart = localtime(time);
my $current_dir = $ENV{"PWD"};

## Get arguments
my %options = ();
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
my $host = $options{h};
my $db2 = $options{d};
my $user = $options{u};
my $passwd = $options{p};
my $in = $options{i};
my $t_ld_key = $options{j};
my $out = $options{o};

## Variables
my $LD_cutoff = 0.5; # JC changed from 0.3 to 0.5 (Weronika used this)
my $dist_cutoff = 40000; # JC changed from 37000 to 40000 (arbitrary - whole num)
my $t_ld_head;
my $t_ld_tail;
my $insert_data;
my $SNP_list;

my @chr;

# Checking for correct LD database.tablename chr format in option j : e.g. var_ld_data.chr_CEU_all
if ($t_ld_key =~ /^(.+chr)(_.+)$/){
    $t_ld_head = $1;
    $t_ld_tail = $2;
}
else {
    print "Error: incorrect table name format for LD R_square data\n";
    exit(0);
}

# Initialize LD hash
my %having_ld = ();
@chr = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M);
foreach my $chr(@chr){
    $having_ld{$chr}=1;
}

## Connect to the database.
my $db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});
$db->do("drop table if exists $out");

## Create output table
$db->do("create table $out like $in");
my $get_dz_chr = $db->prepare("select distinct broad_phenotype, chromosome from $in order by broad_phenotype, chromosome");

## JC couldnt find study_cnt and sample_size
## not really required right now so gloss over
my $get_SNP = $db->prepare("select distinct dbSNP, end from $in where broad_phenotype=? and chromosome=? order by gwas_study_cnt desc, study_cnt_total desc, gwas_sample_size desc, study_cnt desc, sample_size desc;");

$get_dz_chr->execute();

## find LD table for each disease and chromosome
## population is given as an input
## find LD SNP pair
while (my ($dz, $chr)=$get_dz_chr->fetchrow_array()){
    
    ## declare variables and prepare queries
    my $t_ld = $t_ld_head.$chr.$t_ld_tail;
    my $get_LD_SNP = $db->prepare("select distinct dbSNP2 from $t_ld a join $in b ON(a.dbSNP2=b.dbSNP) where a.dbSNP1=? and b.broad_phenotype=? and b.chromosome=? and a.R_square>=$LD_cutoff");
    my $get_SNP_dist = $db->prepare("select distinct dbSNP from $in where broad_phenotype=? and chromosome=? and end<? and end>? and dbSNP != ? and dbSNP not in (select dbSNP2 from $t_ld where dbSNP1=?)");
    my @dbSNP = ();
    my %end = ();
    
    ## get all SNP info from dz LR table (VariMed)
    $get_SNP->execute($dz, $chr);
    while(my ($dbSNP, $end)=$get_SNP->fetchrow_array()){
				push(@dbSNP, $dbSNP);
				$end{$dbSNP} = $end;
    }
    
    $get_SNP->finish();
    
    ## for all the dbSNPs in VariMed
    if (exists $having_ld{$chr}){
				for(my $i=0; $i<@dbSNP; $i++){
	    			my %dbSNP_LD = ();
	    			
	    			## find if they match dbSNP1 in LD table, get dbSNP2
	    			$get_LD_SNP->execute($dbSNP[$i], $dz, $chr);
	    			
	    			## for all the SNPs in VariMed that has an LD with dbSNP1
	    			while (my $dbSNP1=$get_LD_SNP->fetchrow_array()){
								$dbSNP_LD{$dbSNP1}=1; 
	    			}
	    			
	    			$get_LD_SNP->finish();
	    			
	    			## for all the SNPs in VariMed that has an LD by distance
	    			$get_SNP_dist->execute($dz, $chr, $end{$dbSNP[$i]}+$dist_cutoff, $end{$dbSNP[$i]}-$dist_cutoff, $dbSNP[$i], $dbSNP[$i]);
	    			while (my $dbSNP1=$get_SNP_dist->fetchrow_array()){
                $dbSNP_LD{$dbSNP1}=2;
            }
            $get_SNP_dist->finish();
            
            ## if SNP exists to have LD either by dbSNP2 r2 (type 1) or distance (type 2)
            ## SNPs do not seem to be removed by comparing LR... 
	    			for(my $j=$i+1; $j<@dbSNP; $j++){
								if (exists $dbSNP_LD{$dbSNP[$j]}){
		    						print "Remove rs$dbSNP[$j] on chromosome $chr on $dz due to LD with rs$dbSNP[$i] type: $dbSNP_LD{$dbSNP[$j]}\n";
		    						splice(@dbSNP, $j, 1);
		    						$j--;
								}
	    			}
			}
    }
    
    #print "@dbSNP\n"; ## debug
    $SNP_list = "(".join(", ", @dbSNP).")";
    #print "$SNP_list\n"; ## debug
    $insert_data = $db->prepare("insert into $out select * from $in where broad_phenotype=? and chromosome=? and dbSNP in $SNP_list");
    $insert_data->execute($dz, $chr);
    $insert_data->finish();
    #print "insert data broad_phenotype=$dz and chromosome=$chr and dbSNP in ($SNP_list)\n";
}
$get_dz_chr->finish();

## Disconnect from database
$db->disconnect();

## Output time
my $timeend = localtime(time); ## change end variable diff from loop above
print "$out is successfully updated by $0 from $timestart to $timeend.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to keep a single SNP with the largest (gwas_study_size, total_cnt_gwas, gwas_sample_size, study_cnt, sample_size) for any pair of SNPs with R_square>=0.5, or distance of 40kb from any SNP pair with R_square>=0.5 on any disease/chromosome pair.

    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name
    -j table to retrieve the R square value of any SNP pair 
    	 (note naming convention in example; tables to be split by chr)
    	 (note also that only dbSNP1 is used for joining, so LD table has to be commutative, i.e. dbSNP1 has to be found in dbSNP2 as well)
    	 (Note also that R2 is fixed at 0.5 right now)
    -o output table name
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_westome -u rchen1 -p ** -i 20_LR -j var_ld_data.chr_CEU -o 20_LR_uniq
END
}
