#!/usr/bin/env python
import pymysql, sys
from string import Template

## this script adapts Boris' original script to load a file (header-less)
## line by line into a user-defined but hard-coded mySQL database (host,user,passwd) 
## - currently on buttelab-s01
##
## inputs required:
## arg1: filename to be uploaded
## arg2: tablename to be created in database
## arg3: databasename

# Examples python ./load_tbl.py vcf2geno-NA12878-1kgp3-snps.txt NA12878_genotype var_ld_data_1000genomes_phase3

def load_data_from_file(data_file_name, target_table):

    # loads data from file data_file_name into table target_table
    connection = pymysql.connect(host, user, passwd, db, local_infile=True)
    cursor = connection.cursor()
    
    ## Create table as per requirement
    droptable = "DROP TABLE IF EXISTS " + target_table
    cursor.execute(droptable)
		
    createtable = "CREATE TABLE " + target_table + """ (
                 sample_ID  VARCHAR(50) NOT NULL,
                     dbSNP  BIGINT(20),
                     genotype CHAR(2) ); """
    
    cursor.execute(createtable)
    
    ## insert data
    sql = "LOAD DATA LOCAL INFILE \"" + data_file_name + \
          "\" INTO TABLE " + target_table + " CHARACTER SET UTF8 FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"

    cursor.execute(sql)

    ## create primary key and index
    pkey = "ALTER TABLE " + target_table + " ADD PRIMARY KEY(dbSNP)"
    cursor.execute(pkey)
    
    index = "CREATE INDEX genotype on " + target_table + " (genotype)"
    cursor.execute(index)

    ## commit
    connection.commit()
    cursor.close()
    connection.close()

    print("finished loading: ", data_file_name)
    sys.stdout.flush()

    return


def main():
		
    # Examples python ./load_tbl_1kgp3_snps.py vcf2geno-NA12878-1kgp3-snps.txt NA12878_genotype var_ld_data_1000genomes_phase3 
    global host, user, passwd, db, filename
    filename = sys.argv[1]
    tablename = sys.argv[2]
    tablename = tablename.replace("-","") ## for some reason, mySQL churn errors with "-" in table name
    host = "db address here"
    user = "username here"
    passwd = "password here"
    db = sys.argv[3]
    load_data_from_file(filename, tablename)

if __name__ == '__main__':
    main()
