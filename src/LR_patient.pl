#!/usr/bin/perl

# Rong Chen, SMI
# This program is to report LR for a patient
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
getopts("h:d:u:p:i:j:o:g:s:t:", \%options);
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
$patient = $options{j};
$out = $options{o};

if (defined $options{g}){
    $gwas_study_cnt_cutoff = $options{g};
}
if (defined $options{s}){
    $gwas_sample_size_cutoff = $options{s};
}

$out_file = $current_dir."/".$out.".txt";
# Connect to the database.
$db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

$query = qq{drop table if exists $out};
$db->do($query);
$query = qq{create table $out(id int not NULL AUTO_INCREMENT primary key,
			      broad_phenotype varchar(128), index(broad_phenotype),
                              study_id int, index(study_id),
                              pubmed int, index(pubmed),
                              dbSNP bigint, index(dbSNP),
                              Narrow_phenotype varchar(128),
                              genotype varchar(64), index(genotype),
                              LR float, index(LR),
                              total int, index(total),
                              cases varchar(256), index(cases),
                              controls varchar(256), index(controls),
                              Model varchar(128),
			      GWAS varchar(10) DEFAULT NULL,
			      chromosome varchar(2) DEFAULT NULL,
			      end bigint(20) unsigned DEFAULT NULL,
			      fxn_type varchar(128) DEFAULT NULL,
			      symbol varchar(128) DEFAULT NULL,
			      population varchar(128),
			      study_cnt_total int(10) unsigned DEFAULT NULL,
			      gwas_study_cnt int(10) unsigned DEFAULT NULL,
			      gwas_sample_size int(10) unsigned DEFAULT NULL
                              )};
$db->do($query);
$query = qq{alter table $out disable keys};
$db->do($query);

$query = qq{select distinct dbSNP from $in join $patient using(dbSNP)};
$get_dbSNP = $db->prepare($query);

if (defined $options{t}){
    $Pmin = $options{t};
    if ((defined $options{g}) && (defined $options{s})){
	$query = qq{select id, genotype, strand from $in where dbSNP=? and gwas_study_cnt>=$gwas_study_cnt_cutoff and gwas_sample_size>=$gwas_sample_size_cutoff and Pmin<=$Pmin};
    }
    elsif (defined $options{g}){
	$query = qq{select id, genotype, strand from $in where dbSNP=? and gwas_study_cnt>=$gwas_study_cnt_cutoff and Pmin<=$Pmin};
    }
    elsif (defined $options{s}){
	$query = qq{select id, genotype, strand from $in where dbSNP=? and gwas_sample_size>=$gwas_sample_size_cutoff and Pmin<=$Pmin};
    }
    else{
	$query = qq{select id, genotype, strand from $in where dbSNP=? and Pmin<=$Pmin};
    }
}
else {
    if ((defined $options{g}) && (defined $options{s})){
	$query = qq{select id, genotype, strand from $in where dbSNP=? and gwas_study_cnt>=$gwas_study_cnt_cutoff and gwas_sample_size>=$gwas_sample_size_cutoff};
    }
    elsif (defined $options{g}){
	$query = qq{select id, genotype, strand from $in where dbSNP=? and gwas_study_cnt>=$gwas_study_cnt_cutoff};
    }
    elsif (defined $options{s}){
	$query = qq{select id, genotype, strand from $in where dbSNP=? and gwas_sample_size>=$gwas_sample_size_cutoff};
    }
    else{
	$query = qq{select id, genotype, strand from $in where dbSNP=?};
    } 
}
$get_strand = $db->prepare($query);

$query = qq{select genotype from $patient where dbSNP=?};
$get_genotype = $db->prepare($query);

$query = qq{select broad_phenotype, study_id, pubmed, dbSNP, Narrow_phenotype, LR, total, cases, controls, Model, GWAS, chromosome, end, fxn_type, symbol, population, study_cnt_total, gwas_study_cnt, gwas_sample_size from $in where id=?};
$get_data = $db->prepare($query);

$get_dbSNP->execute();
open(OUT, ">$out_file") || die "Error: $! when write $out_file\n";
while ($dbSNP=$get_dbSNP->fetchrow_array()){
    print "dbSNP: $dbSNP\n";
    $id_final = '';
    $get_genotype->execute($dbSNP);
    $genotype_p=$get_genotype->fetchrow_array();
    $get_genotype->finish();
    $get_strand->execute($dbSNP);
    while (($id, $genotype, $strand)=$get_strand->fetchrow_array()){
	$id_final ='';
	if ($strand eq '-'){
	    #print "Before flip: genotype: $genotype\n";
	    $genotype = &flip($genotype);
	    #print "After flip: genotype:$genotype\n";
	}
	if ($genotype eq $genotype_p){
	    $id_final = $id;
	    $get_data->execute($id_final);
	    ($phenotype, $study_id, $pubmed, $dbSNP, $Narrow_phenotype, $LR, $total, $cases, $controls, $Model, $GWAS, $chromosome, $end, $fxn_type, $symbol, $pop, $study_cnt_total, $gwas_study_cnt, $gwas_sample_size)=$get_data->fetchrow_array();
	    print OUT "$id_final\t$phenotype\t$study_id\t$pubmed\t$dbSNP\t$Narrow_phenotype\t$genotype_p\t$LR\t$total\t$cases\t$controls\t$Model\t$GWAS\t$chromosome\t$end\t$fxn_type\t$symbol\t$pop\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
	    $get_data->finish();
	    #print "Matched in 1,id:$id, flip negative strand genotype:$genotype\n";
	}
	else{
	    #print "before reverse: genotype $genotype\n";
	    $genotype2 = join("", reverse(split(//, $genotype)));
	    #print "After reverse: genotype $genotype2\n";
	    if (($genotype2 eq $genotype_p) && ($genotype ne $genotype2)){
		$id_final=$id;
		#print "original: $genotype, reverse: $genotype2, patient:$genotype_p\n";
		$get_data->execute($id_final);
		($phenotype, $study_id, $pubmed, $dbSNP, $Narrow_phenotype, $LR, $total, $cases, $controls, $Model, $GWAS, $chromosome, $end, $fxn_type, $symbol, $pop, $study_cnt_total, $gwas_study_cnt, $gwas_sample_size)=$get_data->fetchrow_array();
		print OUT "$id_final\t$phenotype\t$study_id\t$pubmed\t$dbSNP\t$Narrow_phenotype\t$genotype_p\t$LR\t$total\t$cases\t$controls\t$Model\t$GWAS\t$chromosome\t$end\t$fxn_type\t$symbol\t$pop\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
		$get_data->finish();
		
	    }
	    else {
		#print "Unmatched genotype at dbSNP $dbSNP, patient $genotype_p, disease $genotype\n";
	    }
	}
    }
    $get_strand->finish();
    if ($id_final eq ''){
	print "Unmatched genotype at dbSNP $dbSNP, patient $genotype_p, disease $genotype\n";
    }
}
$get_dbSNP->finish();
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
This program is to report match patient genotype with genotypes in the LR table
    Available arguments:
    -h mysql host name
    -d database name
    -u user name
    -p password
    -i input LR table name
    -j patient genotype input table
    -o output table name
    -g cutoff for the number of GWAS studies (>=, optional)
    -s cutoff for the total sample size of GWAS studies (>=, optional)
    -t cutoff for minimum P values (<=, optional)
    
  Usage: $0 -h bmir-db1 -d proj_dz_risk_westome -u rchen1 -p ** -i proj_patient_risk.LR_caucasian -j 20_genotype -o 20_out3 -g 2 -s 2000 -t 0.000001
END
}
