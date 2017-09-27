#!/usr/bin/perl

# Rong Chen, SMI
# This program is to create a LR table for case vs. control studies
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
$query = qq{create table $out(broad_phenotype varchar(128), index(broad_phenotype),
			      dbSNP bigint, index(dbSNP),
			      genotype varchar(64), index(genotype),
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

$query = qq{select distinct broad_phenotype, dbSNP, genotype, chromosome, end, fxn_type, symbol, study_cnt_total, gwas_study_cnt, gwas_sample_size from $in};
$get_dz_snp = $db->prepare($query);

$query = qq{select distinct study_id from $in where broad_phenotype=? and dbSNP=? and genotype=?};
$get_study_id = $db->prepare($query);

$query = qq{select count(distinct study_id) from $in where broad_phenotype=? and dbSNP=? and genotype=?};
$get_study_cnt = $db->prepare($query);

$query = qq{select count(*) from $in where broad_phenotype=? and dbSNP=? and study_id=? and genotype=?};
$get_LR_count = $db->prepare($query);

$query = qq{select LR, total from $in where broad_phenotype=? and dbSNP=? and study_id=? and genotype=?};
$get_LR = $db->prepare($query);

$query = qq{select LR, total from $in where broad_phenotype=? and dbSNP=? and genotype=?};
$get_LR_single_study = $db->prepare($query);

open(OUT, ">$out_file") || die "Error: $! when write $out_file";
$get_dz_snp->execute();
while (($broad_phenotype, $dbSNP, $genotype, $chr, $end, $fxn_type, $symbol, $study_cnt_total, $gwas_study_cnt, $gwas_sample_size)=$get_dz_snp->fetchrow_array()){
    $get_study_cnt->execute($broad_phenotype, $dbSNP, $genotype);
    $study_cnt = $get_study_cnt->fetchrow_array();
    $get_study_cnt->finish();
    if ($study_cnt==1){
	$get_LR_single_study->execute($broad_phenotype, $dbSNP, $genotype);
	($LR_final, $size_final)=$get_LR_single_study->fetchrow_array();
	$get_LR_single_study->finish();
	print OUT "$broad_phenotype\t$dbSNP\t$genotype\t$LR_final\t$study_cnt\t$size_final\t$chr\t$end\t$fxn_type\t$symbol\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
	next;
    }
    else {
	$get_study_id->execute($broad_phenotype, $dbSNP, $genotype);
	$LR_sum = 0;
	$size_sum = 0;
	$size_final = 0;
	$LR_final=0;
	while ($study_id=$get_study_id->fetchrow_array()){
	    print "$broad_phenotype, $dbSNP, $genotype, $study_id\n";
	    $get_LR_count->execute($broad_phenotype, $dbSNP, $study_id, $genotype);
	    $LR_cnt = $get_LR_count->fetchrow_array();
	    $get_LR_count->finish();
	    $LR1=0;
	    $size1=0;
	    # get $LR1 and $size1 for every dz/SNP/study_id
	    $get_LR->execute($broad_phenotype, $dbSNP, $study_id, $genotype);
	    if ($LR_cnt==1){
		($LR1, $size1)=$get_LR->fetchrow_array();
		$get_LR->finish();
	    }
	    else {
		$LR1_1 = 0;
		$size1_1 = 0;
		$LR1_sum = 0;
		$size1_sum = 0;
		$size1_max = 0;
		while (($LR1_1, $size1_1)=$get_LR->fetchrow_array()){
		    $LR1_sum += log($LR1_1)*sqrt($size1_1+1);
		    $size1_sum += sqrt($size1_1+1);
		    if ($size1_max < $size1_1){
			$size1_max=$size1_1;
		    }
		}
		$get_LR->finish();
		$LR1 = exp($LR1_sum/$size1_sum);
		$size1 = $size1_max;
	    }
	    print "size1: $size1, LR_sum: $LR_sum\n";
	    $LR_sum += log($LR1)*sqrt($size1);
	    $size_sum += sqrt($size1);
	    $size_final += $size1;
	    print "size1: $size1, size_sum: $size_sum, LR1: $LR1, LR_sum: $LR_sum\n";
	}
    }
    $get_study_id->finish();
    if ($size_sum>0){
	$LR_final = exp($LR_sum/$size_sum);
	print OUT "$broad_phenotype\t$dbSNP\t$genotype\t$LR_final\t$study_cnt\t$size_final\t$chr\t$end\t$fxn_type\t$symbol\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
    }
    else {
	print "size_sum= $size_sum in $broad_phenotype\t$dbSNP\t$genotype\t$LR_final\t$study_cnt\t$size_final\n";
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
