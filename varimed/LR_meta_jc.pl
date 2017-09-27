#!/usr/bin/perl

# This program is to create an all_snps info merged table
# For each broad_phenotype, for each genotype, and each SNP, 
# one integrated LR is calculated for each disease/SNP pair
# by combining log(LR) across studies weighted by size of study

# The original code was built by Rong Chen, Buttelab, Stanford (08/19/2009).
# This is modified by Jieming Chen, Buttelab, UCSF (03/15/17)
# 1) adding comments to make code more readable
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
my $start = localtime(time);
my $current_dir = $ENV{"PWD"};

## Get arguments
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

## Connect to the database.
my $db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

## Create metadata data table from VariMed based on patient genotype
## indexes are created
my $query = qq{drop table if exists $out};
$db->do($query);
$query = qq{create table $out(broad_phenotype varchar(128), index(broad_phenotype),
			      dbSNP bigint, index(dbSNP),
			      genotype char(2), index(genotype), # JC varchar2char
			      LR float, index(LR),
			      study_cnt int,
			      sample_size int,
			      chromosome varchar(2) DEFAULT NULL,
			      end bigint(20) unsigned DEFAULT NULL,
			      fxn_type varchar(128) DEFAULT NULL,
			      symbol varchar(128) DEFAULT NULL,
			      study_cnt_total int(10) unsigned DEFAULT NULL,
			      gwas_study_cnt int(10) unsigned DEFAULT NULL,
			      gwas_sample_size int(10) unsigned DEFAULT NULL
			      )
			  };
$db->do($query);

## get_dz_snp query
$query = qq{select distinct broad_phenotype, dbSNP, genotype, chromosome, end, fxn_type, symbol, study_cnt_total, gwas_study_cnt, gwas_sample_size from $in};
my $get_dz_snp = $db->prepare($query) or die "Couldn't prepare statement: + $DBI::errstr; stopped";

## get_study_id query
$query = qq{select distinct study_id from $in where broad_phenotype=? and dbSNP=? and genotype=?};
my $get_study_id = $db->prepare($query) or die "Couldn't prepare statement: + $DBI::errstr; stopped";

## get_study_cnt query
$query = qq{select count(distinct study_id) from $in where broad_phenotype=? and dbSNP=? and genotype=?};
my $get_study_cnt = $db->prepare($query) or die "Couldn't prepare statement: + $DBI::errstr; stopped"; 

## get_LR_count query
$query = qq{select count(*) from $in where broad_phenotype=? and dbSNP=? and study_id=? and genotype=?};
my $get_LR_count = $db->prepare($query) or die "Couldn't prepare statement: + $DBI::errstr; stopped";

## get_LR query
$query = qq{select LR, total from $in where broad_phenotype=? and dbSNP=? and study_id=? and genotype=?};
my $get_LR = $db->prepare($query) or die "Couldn't prepare statement: + $DBI::errstr; stopped";

## get_LR_single_study
$query = qq{select LR, total from $in where broad_phenotype=? and dbSNP=? and genotype=?};
my $get_LR_single_study = $db->prepare($query) or die "Couldn't prepare statement: + $DBI::errstr; stopped";

## Define variables
my $broad_phenotype;
my $dbSNP;
my $genotype;
my $LR_final;
my $study_cnt;
my $chr;
my $end;
my $fxn_type;
my $symbol;
my $study_cnt_total;
my $gwas_study_cnt;
my $gwas_sample_size;
my $size_sum;
my $size_final;
my $size1;
my $size1_1;
my $size1_max;
my $size1_sum;
my $LR_cnt;
my $LR_sum;
my $LR1_sum;
my $LR1;
my $LR1_1;
my $study_id;

## tmp file for upload to database later
open(OUT, ">$out_file") || die "Error: $! when write $out_file";

## get_dz_snp; get disease and study SNPs metadata for each dbSNP rsID found in the patient genotype table
$get_dz_snp->execute() or die "Couldn't execute statement: $DBI::errstr; stopped";

while (my ($broad_phenotype, $dbSNP, $genotype, $chr, $end, $fxn_type, $symbol, $study_cnt_total, $gwas_study_cnt, $gwas_sample_size)=$get_dz_snp->fetchrow_array()){
		
		## count number of studies for each broad_phenotype
    $get_study_cnt->execute($broad_phenotype, $dbSNP, $genotype) or die "Couldn't execute statement: $DBI::errstr; stopped";
    $study_cnt = $get_study_cnt->fetchrow_array();
    $get_study_cnt->finish();
    
    ## if just 1 study, print to output file
    if ($study_cnt==1){
				$get_LR_single_study->execute($broad_phenotype, $dbSNP, $genotype) or die "Couldn't execute statement: $DBI::errstr; stopped";
				($LR_final, $size_final)=$get_LR_single_study->fetchrow_array();
				$get_LR_single_study->finish();
				print OUT "$broad_phenotype\t$dbSNP\t$genotype\t$LR_final\t$study_cnt\t$size_final\t$chr\t$end\t$fxn_type\t$symbol\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
				next;
    }
    else { ## otherwise if many studies, print to stdout (which should be piped to a log file) 
    			 ## on how LR is combined
    			 ## cumulative for each broad phenotype in stdout/logfile
						$get_study_id->execute($broad_phenotype, $dbSNP, $genotype) or die "Couldn't execute statement: $DBI::errstr; stopped";
						$LR_sum = 0;
						$size_sum = 0;
						$size_final = 0;
						$LR_final = 0;
						while ($study_id=$get_study_id->fetchrow_array()){
						    print "\n$broad_phenotype, dbSNP=$dbSNP, $genotype, study_id=$study_id\n";
						    $get_LR_count->execute($broad_phenotype, $dbSNP, $study_id, $genotype) or die "Couldn't execute statement: $DBI::errstr; stopped";
						    $LR_cnt = $get_LR_count->fetchrow_array();
						    $get_LR_count->finish();
						    $LR1 = 0;
						    $size1 = 0;
						    # get $LR1 and $size1 for every dz/SNP/study_id
						    $get_LR->execute($broad_phenotype, $dbSNP, $study_id, $genotype) or die "Couldn't execute statement: $DBI::errstr; stopped";
						    
						    ## if only 1 LR value
						    if ($LR_cnt==1){
												($LR1, $size1)=$get_LR->fetchrow_array();
												$get_LR->finish();
			    			}
					    	else { ## if multiple LR values, combine all LR values to get a single LR
												$LR1_1 = 0;
												$size1_1 = 0;
												$LR1_sum = 0;
												$size1_sum = 0;
												$size1_max = 0;
												while (($LR1_1, $size1_1)=$get_LR->fetchrow_array()){
														
														## weight log LR with indiv study size then sum it
												    $LR1_sum += log($LR1_1)*sqrt($size1_1+1);
												    $size1_sum += sqrt($size1_1+1);
												    if ($size1_max < $size1_1){
																$size1_max=$size1_1;
				    								}
												}
												## then normalize LR sum by all studies size from that broad phenotype
												$get_LR->finish();
												$LR1 = exp($LR1_sum/$size1_sum);
												$size1 = $size1_max;
			    			}
						    print "size_max: $size1, LR_sum_prev: $LR_sum\n";
						    $LR_sum += log($LR1)*sqrt($size1); ## there is no +1 here
						    $size_sum += sqrt($size1);
						    $size_final += $size1;
						    print "size_max: $size1, size_sum_sqrt: $size_sum, LR1_curr=ln(LR)*sqrt(size): $LR1, LR_sum_curr: $LR_sum\n";
						 }
    }
    $get_study_id->finish();
    
    ## once LR combined, if study has samples present, print to tmp file
    ## else print as error to STDOUT 
    if ($size_sum>0){
					$LR_final = exp($LR_sum/$size_sum);
					print OUT "$broad_phenotype\t$dbSNP\t$genotype\t$LR_final\t$study_cnt\t$size_final\t$chr\t$end\t$fxn_type\t$symbol\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
    }
    else { 
					print "size_sum= $size_sum (<=0) in $broad_phenotype, dbSNP=rs$dbSNP, $genotype, LR_average_combined=$LR_final, $study_cnt$study_cnt, size_final_combined=$size_final (<0)\n";
    }
}
$get_dz_snp->finish();
close(OUT);

$query = qq{load data local infile \'$out_file\' into table $out};
$db->do($query);
unlink($out_file);

$query = qq{optimize table $out};
$db->do($query);

$end = localtime(time);
print "$out is successfully updated by $0 from $start to $end.\n";
exit(0);

#############################################################################
# print usage message
#############################################################################
sub printMessage()
{
    print<<END;
This program is to compute the final LR for each disease/SNP pair
    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input table name
    -o output table name
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_westome -u rchen1 -p ** -i 20_snps_1 -o 20_snps_2
END
}
