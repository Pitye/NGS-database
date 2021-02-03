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


sample1_name='AA11-KIGEN-SF22'
sample2_name='AA11-KIGEN-SF26'
markers_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5,'TPOX' : 0.5,'D2S441' : 0.5,'D2S1338' : 0.5,'D3S1358' : 0.5,'D4S2408' : 0.5,'FGA' : 0.5,'D5S818' : 0.5,'CSF1PO' : 0.5,'D6S1043' : 0.5,'D7S820' : 0.5,'D8S1179' : 0.5,'D9S1122' : 0.5,'D10S1248' : 0.5,'TH01': 0.5,'vWA' : 0.5,'D12S391' : 0.5,'D13S317' : 0.5,'PentaE' : 0.5,'D16S539' : 0.5,'D17S1301' : 0.5,'D18S51' : 0.5,'D19S433' : 0.5,'D20S482' : 0.5,'D21S11' : 0.5,'PentaD' : 0.5,'D22S1045' : 0.5}
markers = list(markers_dict.keys())
columns = ['allele', 'seq_name', 'PubMed_ID', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']


sample1_records = {marker: {column: [] for column in columns} for marker in markers}
sample2_records = {marker: {column: [] for column in columns} for marker in markers}

db=MySQLdb.connect("localhost", "root", "Coufalka*1", "NGS_FORENSIC")
c = db.cursor()
sql_select_Query1 = "SELECT * FROM ngs_forensic.nomen_freq_autostrdata_flankingreg where sample_name = '%s'" % (sample1_name)
c.execute(sql_select_Query1)
records1 = c.fetchall()

print("Total number of rows for sample1: ", c.rowcount)

print("\nPrinting sample " + sample1_name )
for row in records1:
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
        sample1_records[row[1]][column].append(row[column_index])
        column_index += 1
        
sql_select_Query2 = "SELECT * FROM ngs_forensic.nomen_freq_autostrdata_flankingreg where sample_name = '%s'" % (sample2_name)
c.execute(sql_select_Query2)
records2 = c.fetchall()

print("Total number of rows for sample2: ", c.rowcount)

print("\nPrinting sample " + sample2_name )
for row in records2:
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
        sample2_records[row[1]][column].append(row[column_index])
        column_index += 1        

print (sample1_records['D1S1656']) 
freq1 = sample1_records['D1S1656']['frequency'][0] + sample1_records['D1S1656']['frequency'][1]

print (sample2_records['D1S1656']) 
freq2 = sample2_records['D1S1656']['frequency'][0] + sample2_records['D1S1656']['frequency'][1]
db.close()
print('all done')
print (freq1, freq2)