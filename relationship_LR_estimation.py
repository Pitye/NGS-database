import pandas as pd
import os
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
import MySQLdb
#import mysql.connector
#from mysql.connector import Error

### *** parameter settings ***
sample_names=['MS9-0,0065', 'AA11-KIGEN-SF26']

markers_dict = {'D1S1656': 0.5,'TPOX' : 0.5,'D2S441' : 0.5,'D2S1338' : 0.5,'D3S1358' : 0.5,'D4S2408' : 0.5,'FGA' : 0.5,'D5S818' : 0.5,'CSF1PO' : 0.5,'D6S1043' : 0.5,'D7S820' : 0.5,'D8S1179' : 0.5,'D9S1122' : 0.5,'D10S1248' : 0.5,'TH01': 0.5,'vWA' : 0.5,'D12S391' : 0.5,'D13S317' : 0.5,'PentaE' : 0.5,'D16S539' : 0.5,'D17S1301' : 0.5,'D18S51' : 0.5,'D19S433' : 0.5,'D20S482' : 0.5,'D21S11' : 0.5,'PentaD' : 0.5,'D22S1045' : 0.5}
markers = list(markers_dict.keys())
columns = ['allele', 'seq_name', 'PubMed_ID', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']

length_polymorphism_estimation = 'YES'
use_external_file_for_missing_markers = 'YES'
directory_path = os.path.normpath('C:/NGS_forensic_database/relationship_LR_estimation')


### *** create empty dict for data from mySQL dtb ***
#sample1_records = {marker: {column: [] for column in columns} for marker in markers}
#sample2_records = {marker: {column: [] for column in columns} for marker in markers}
sample_records = {sample: {marker: {column: [] for column in columns} for marker in markers} for sample in sample_names}
missed_markers = {sample: [] for sample in sample_names}
    
### *** get data from from mySQL dtb ***
db=MySQLdb.connect("localhost", "root", "Coufalka*1", "NGS_FORENSIC")
c = db.cursor()

for sample in sample_names:
    
    sql_select_Query = "SELECT * FROM ngs_forensic.nomen_freq_autostrdata_flankingreg where sample_name = '%s'" % (sample)
    c.execute(sql_select_Query)
    records = c.fetchall()

    print("Total number of rows for sample1: ", c.rowcount)

    print("\nPrinting sample " + sample )
    for row in records:
        #print("marker = ", row[1])
        #print("allele  = ", row[2]) 
        #print("seq_name  = ", row[3])
        #print("PubMed_ID  = ", row[4])
        #print("sequence  = ", row[5])
        #print("no_reads  = ", row[6])
        #print("CE_validation  = ", row[7])
        #print("head_id  = ", row[8])
        #print("avg_no_reads  = ", row[9])
        ##print("count_seq  = ", row[10])
        print("frequency  = ", row[11], "\n")
    
        column_index = 2
        for column in columns:
            sample_records[sample][row[1]][column].append(row[column_index])
            column_index += 1
db.close()       
#sql_select_Query2 = "SELECT * FROM ngs_forensic.nomen_freq_autostrdata_flankingreg where sample_name = '%s'" % (sample2_name)
#c.execute(sql_select_Query2)
#records2 = c.fetchall()

#print("Total number of rows for sample2: ", c.rowcount)

#print("\nPrinting sample " + sample2_name )
#for row in records2:
    #print("marker = ", row[1])
    #print("allele  = ", row[2])
    #print("seq_name  = ", row[3])
    #print("PubMed_ID  = ", row[4])
    #print("sequence  = ", row[5])
    #print("no_reads  = ", row[6])
    #print("CE_validation  = ", row[7])
    #print("head_id  = ", row[8])
    #print("avg_no_reads  = ", row[9])
    ##print("count_seq  = ", row[10])
    #print("frequency  = ", row[11], "\n")
    
    #column_index = 2
    #for column in columns:
        #sample2_records[row[1]][column].append(row[column_index])
        #column_index += 1
        
### *** check which markers don't have any data ***
for sample in sample_names:    
    for marker in markers:
        if not sample_records[sample][marker]['allele']:
            missed_markers[sample].append(marker)

#print (sample1_records['D1S1656'])
#if not sample1_records['D1S1656']['allele']:
    #print ('D1S1656 is empty')
    #del sample1_records['D1S1656']
#freq1 = sample1_records['D1S1656']['frequency'][0] + sample1_records['D1S1656']['frequency'][1]
for sample in sample_names:
    print (list(sample_records[sample].keys()))
    print (sample_records[sample]['D1S1656']) 
    #freq2 = sample_records[sample]['D1S1656']['frequency'][0] + sample2_records['D1S1656']['frequency'][1]
print (missed_markers)    
print('all done')
#print (freq1, freq2)

# Checking if 4 exists in list  
# using set() + in 
#test_list_set = set(test_list_set) 
#if 4 in test_list_set : 
    #print ("Element Exists") 
