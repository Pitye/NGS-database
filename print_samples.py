import pandas as pd
import os
import sys
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
import MySQLdb
import math
from datetime import datetime



directory_reports = os.path.normpath('C:/NGS_forensic_database/print_samples/prints')
path_samples_print = os.path.normpath('C:/NGS_forensic_database/print_samples/print_sample.txt')
path_locus_order = os.path.normpath('C:/NGS_forensic_database/print_samples/locus_order.txt')
table_list = ['AutoSTR', 'Y-STR']
selected_columns = ['allele', 'PubMed_ID','no_reads' ]
columns = ['allele', 'seq_name', 'PubMed_ID', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']
columnsY = ['allele', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']
sample_names = []
locus_list = []
markersAuto = ['D1S1656','TPOX' ,'D2S441' ,'D2S1338' ,'D3S1358' ,'D4S2408' ,'FGA' ,'D5S818' ,'CSF1PO' ,'D6S1043' ,'D7S820' ,'D8S1179' ,'D9S1122' ,'D10S1248' ,'TH01','vWA' ,'D12S391' ,'D13S317' ,'PentaE' ,'D16S539' ,'D17S1301' ,'D18S51' ,'D19S433' ,'D20S482' ,'D21S11' ,'PentaD' ,'D22S1045']
markersY = ['DYF387S1', 'DYS19', 'DYS385a-b', 'DYS389I', 'DYS389II', 'DYS390', 'DYS391', 'DYS392', 'DYS437', 'DYS438', 'DYS439', 'DYS448', 'DYS460', 'DYS481', 'DYS505', 'DYS522', 'DYS533', 'DYS549', 'DYS570', 'DYS576', 'DYS612', 'DYS635', 'DYS643', 'Y-GATA-H4']
markers = ['D1S1656','TPOX' ,'D2S441' ,'D2S1338' ,'D3S1358' ,'D4S2408' ,'FGA' ,'D5S818' ,'CSF1PO' ,'D6S1043' ,'D7S820' ,'D8S1179' ,'D9S1122' ,'D10S1248' ,'TH01','vWA' ,'D12S391' ,'D13S317' ,'PentaE' ,'D16S539' ,'D17S1301' ,'D18S51' ,'D19S433' ,'D20S482' ,'D21S11' ,'PentaD' ,'D22S1045', \
          'DYF387S1', 'DYS19', 'DYS385a-b', 'DYS389I', 'DYS389II', 'DYS390', 'DYS391', 'DYS392', 'DYS437', 'DYS438', 'DYS439', 'DYS448', 'DYS460', 'DYS481', 'DYS505', 'DYS522', 'DYS533', 'DYS549', 'DYS570', 'DYS576', 'DYS612', 'DYS635', 'DYS643', 'Y-GATA-H4']
 


### *** get locus order from file                    
for locus in csv.reader(open(path_locus_order)):
    if len(locus) == 1:
        locus_list.append(locus[0])
    else:
        print ('ERROR: wrong format ' + path_locus_order)    

### *** get samples from file
for sample in csv.reader(open(path_samples_print), delimiter=';'):
    if len(sample) == 1:
        sample_names.append(sample[0])
    else:
        print ('ERROR: wrong format ' + path_samples_print)


sample_records = {sample: {marker: {column: [] for column in columns} for marker in markers} for sample in sample_names}
        
### *** get data from from mySQL dtb - sample_records ***
db=MySQLdb.connect("localhost", "root", "Coufalka*1", "NGS_FORENSIC")
c = db.cursor()
if 'AutoSTR' in table_list:
    for sample in sample_names:
    
        sql_select_Query = "SELECT * FROM ngs_forensic.nomen_freq_autostrdata_flankingreg where sample_name = '%s'" % (sample)
        c.execute(sql_select_Query)
        records = c.fetchall()

        for row in records:
        #marker, row[1], allele, row[2], seq_name, row[3], PubMed_ID, row[4], sequence, row[5], no_reads, row[6], CE_validation, row[7], head_id, row[8], avg_no_reads, row[9], count_seq, row[10], frequency, row[11]

            column_index = 2
            for column in columns:
                sample_records[sample][row[1]][column].append(row[column_index])
                column_index += 1

    
if 'Y-STR' in table_list:        
    for sample in sample_names:
    
        sql_select_Query = "SELECT * FROM ngs_forensic.freq_y_strdata_flankingreg where sample_name = '%s'" % (sample)
        c.execute(sql_select_Query)
        records = c.fetchall()

        for row in records:
        #marker, row[1], allele, row[2], sequence, row[3], no_reads, row[4], CE_validation, row[5], head_id, row[6], avg_no_reads, row[7], count_seq, row[8], frequency, row[9]

            column_index = 2
            for column in columnsY:
                sample_records[sample][row[1]][column].append(row[column_index])
                column_index += 1    
        
db.close()         
#print (sample_records)        
for sample in sample_names:        
    for marker in markers:
        for column in columns:
            if sample_records[sample][marker][column] == []:
                sample_records[sample][marker][column] = ['null', 'null']
            if len(sample_records[sample][marker][column]) == 1:        
                sample_records[sample][marker][column].append('null')

### *** create report ***                 
col4print = ""
for column in selected_columns:
    col4print = col4print + "," + column + "," + column
now = datetime.now()
dt_string = now.strftime("%Y%m%d%H%M%S")
report_path_name = directory_reports + os.path.normpath('/') + 'report_print_' + dt_string + '.csv'
f = open (report_path_name, 'w+')
#f.writelines(["\nsample_name,marker" + col4print])
for sample in sample_names:
    f.writelines(["sample_name,marker" + col4print])
    print (sample + ' is printed')
    for marker in locus_list:
        row = ","
        for column in selected_columns:
            for i in range (2):
                value = str(sample_records[sample][marker][column][i]) 
                #print (value)
                if value == 'null':
                    row = row + ","
                else:
                    row = row + value + "," 
                #print (row)
        #print (sample + "," + marker + row)
        
        f.writelines(["\n" + sample + "," + marker + row])
    f.writelines(["\n\n"])
f.close()        
print ('done')        
