#!/usr/bin/perl

# This program is to report LR for a patient
# OUTPUT is a out_snps_all_tmp table
#
# The original code was built by Rong Chen, Buttelab, Stanford (08/19/2009).
# This is modified by Jieming Chen, Buttelab, UCSF (03/15/17), by
# 1) adding "use warnings" and "use strict"
# 2) declaring variables, and adding try and if statements to catch exceptions
# 3) adding comments to make code more readable
# 4) adding if statements to SKIP rows with genotype information 
# 5) TODO - maybe should skip rows with no strand info, but because previous version retains
#    them, and there are no rows with '-' strand, we assume EXPLICITLY that they are '+'
# 6) genotype data type in the output database is now char(2)
#    TODO - need to change this for non-SNPs

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
my $host = $options{h};
my $db2 = $options{d};
my $user = $options{u};
my $passwd = $options{p};
my $in = $options{i};
my $patient = $options{j};
my $out = $options{o};
my $gwas_study_cnt_cutoff; 
my $gwas_sample_size_cutoff;

$| = 1; ## no buffering for troubleshooting

if (defined $options{g}){
    $gwas_study_cnt_cutoff = $options{g};
}
if (defined $options{s}){
    $gwas_sample_size_cutoff = $options{s};
}

my $out_file = $current_dir."/".$out.".txt";

## Connect to the database.
my $db = DBI->connect("DBI:mysql:database=$db2;host=$host",
                        $user, $passwd, {'RaiseError' => 1});

## Create new patient_LR table
## Indexes are created
my $query = qq{drop table if exists $out};
$db->do($query);
$query = qq{create table $out(id int not NULL AUTO_INCREMENT primary key,
			      broad_phenotype varchar(128), index(broad_phenotype),
                              study_id int, index(study_id),
                              pubmed int, index(pubmed),
                              dbSNP bigint, index(dbSNP),
                              Narrow_phenotype varchar(128),
                              genotype char(2), index(genotype), # JC varchar2char
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

## LR_patient runs the longest
## find out which query runs the longest and see if we can optimize
#my $time1; # timedebug
#my $time2; # timedebug
#my $time3; # timedebug
#my $time4; # timedebug
#my $time5; # timedebug

## Join patient_LR table with patient_genotype table
$query = qq{select distinct dbSNP from $in join $patient using(dbSNP)};
my $get_dbSNP = $db->prepare($query);


## get_strand query (based on optional args -g -s -t)
if (defined $options{t}){
		    my $Pmin = $options{t};
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

my $get_strand = $db->prepare($query);

## get_genotype query
$query = qq{select genotype from $patient where dbSNP=?};
my $get_genotype = $db->prepare($query);

## get_data query
$query = qq{select broad_phenotype, study_id, pubmed, dbSNP, Narrow_phenotype, LR, total, cases, controls, Model, GWAS, chromosome, end, fxn_type, symbol, population, study_cnt_total, gwas_study_cnt, gwas_sample_size from $in where id=?};
my $get_data = $db->prepare($query);

## Based on each dbSNP rsID in the patient_genotype table, loop to get_data in varimed table
#$time1 = localtime(time); ## timedebug
$get_dbSNP->execute();
#$time2 = localtime(time); ## timedebug
#print "##JCtime## get_dbSNP takes:\n"; ## timedebug
#print "##JCtime## Started: $time1 ### Ended: $time2\n"; ## timedebug
my $phenotype;
my $study_id;
my $pubmed;
my $dbSNP;
my $Narrow_phenotype;
my $LR;
my $total;
my $cases;
my $controls;
my $Model;
my $GWAS;
my $chromosome;
my $end;
my $fxn_type;
my $symbol;
my $pop;
my $study_cnt_total;
my $gwas_study_cnt;
my $gwas_sample_size;
my $genotype;
my $genotype2;
my $id;
my $strand;
open(OUT, ">$out_file") || die "Error: $! when write $out_file\n";

#$time1 = localtime(time); ## timedebug
while (my $dbSNP=$get_dbSNP->fetchrow_array()){
    print "dbSNP: $dbSNP\n";
    my $id_final = '';
    
    $get_genotype->execute($dbSNP);
    my $genotype_p=$get_genotype->fetchrow_array();
    $get_genotype->finish();
#    $time2 = localtime(time); ## timedebug
		
    $get_strand->execute($dbSNP);
#    $time3 = localtime(time); ## timedebug
		
		
    while (($id, $genotype, $strand)=$get_strand->fetchrow_array()){
				$id_final ='';
				
				#print "JMJMJM id=$id; genotype=$genotype; strand=$strand;\n"; ## JC debug
				## JC $strand missing? assume '+' for now since no '-' and previous version let it pass
				if(!defined($strand) || ($strand eq ''))
				{
					$strand = '+';
					#print "Missing strand at dbSNP $dbSNP, patient $genotype_p, disease $genotype\n";
					#next;
				}
				
				## JC $genotype missing? skip rows with no genotypes
    		if((!defined($genotype)) || ($genotype eq ''))
    		{
					$genotype = 'NA';
					print "Missing genotype at dbSNP $dbSNP, patient $genotype_p, genotype $genotype\n";
					next;
				}
				
				if ($strand eq '-'){
	    			#print "Before flip: genotype: $genotype\n";
	    			$genotype = &flip($genotype);
	    			#print "After flip: genotype:$genotype\n";
				}
				
				#print "dbSNP $dbSNP, patient $genotype_p, disease $genotype\n"; ## JC debug
				## genotype == patient genotype
				if ($genotype eq $genotype_p){
			    	$id_final = $id;

			    	$get_data->execute($id_final);
			    	($phenotype, $study_id, $pubmed, $dbSNP, $Narrow_phenotype, $LR, $total, $cases, $controls, $Model, $GWAS, $chromosome, $end, $fxn_type, $symbol, $pop, $study_cnt_total, $gwas_study_cnt, $gwas_sample_size)=$get_data->fetchrow_array();
#			    	$time4 = localtime(time); ## timedebug
			    	
			    	## JC some of the information were NULL
			    	if(!defined($fxn_type)){ $fxn_type = 'NA'; }
			    	if(!defined($symbol)){ $symbol = 'NA'; }
			    	
			    	print OUT "$id_final\t$phenotype\t$study_id\t$pubmed\t$dbSNP\t$Narrow_phenotype\t$genotype_p\t$LR\t$total\t$cases\t$controls\t$Model\t$GWAS\t$chromosome\t$end\t$fxn_type\t$symbol\t$pop\t$study_cnt_total\t$gwas_study_cnt\t$gwas_sample_size\n";
			    	$get_data->finish();
			    	#print "Matched in 1,id:$id, flip negative strand genotype:$genotype\n";
				}
				else{ ## might be diff order in het genotype
	   		 		#print "before reverse: genotype $genotype\n";
	    			$genotype2 = join("", reverse(split(//, $genotype)));
	    			#print "After reverse: genotype $genotype2\n";
	    			if (($genotype2 eq $genotype_p) && ($genotype ne $genotype2)){
								$id_final=$id;
								#print "original: $genotype, reverse: $genotype2, patient:$genotype_p\n";

								$get_data->execute($id_final);
#								$time5 = localtime(time); ## timedebug
								
								($phenotype, $study_id, $pubmed, $dbSNP, $Narrow_phenotype, $LR, $total, $cases, $controls, $Model, $GWAS, $chromosome, $end, $fxn_type, $symbol, $pop, $study_cnt_total, $gwas_study_cnt, $gwas_sample_size)=$get_data->fetchrow_array();
								
								## JC some of the information were NULL
			    			if(!defined($fxn_type)){ $fxn_type = 'NA'; }
			    			if(!defined($symbol)){ $symbol = 'NA'; }
			    	
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
					#print "Unmatched genotype at dbSNP $dbSNP, patient $genotype_p, disease $genotype\n";
					## JC there is no genotype on those with missing data: either $id, $genotype, $strand, or options: $gwas_study_cnt_cutoff, $gwas_sample_size_cutoff, or $Pmin (the script was not able to execute get_strand)
					print "Missing genotype data at dbSNP $dbSNP, disease NA for $phenotype; patient $genotype_p \n";
    }
}
$get_dbSNP->finish();
close(OUT);

## Upload tmp file to database
$query = qq{load data local infile \'$out_file\' into table $out};
$db->do($query);
unlink($out_file);

## set up keys and index
$query = qq{alter table $out enable keys};
$db->do($query);
$query = qq{optimize table $out};
$db->do($query);

$end = localtime(time);

#print "##JCtime## Started while loop: $time1\n"; ## timedebug
#print "##JCtime## get_genotype: $time2\n"; ## timedebug
#print "##JCtime## get_strand  : $time3\n"; ## timedebug
#print "##JCtime## get_data    : $time4\n"; ## timedebug
#print "##JCtime## get_data    : $time5\n"; ## timedebug

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
This program is to report match patient genotype table with the LR table using rsID
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
