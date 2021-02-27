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



path_samples_delete = os.path.normpath('C:/NGS_forensic_database/delete_samples/delete_sample.txt')
table_list = ['AutoSTR', 'Y-STR']

sample_names = []

### *** get samples from file
for sample in csv.reader(open(path_samples_delete), delimiter=';'):
    if len(sample) == 1:
        sample_names.append(sample[0])
    else:
        print ('ERROR: wrong format ' + path_samples_print)

### *** delete data from from mySQL dtb  ***
db=MySQLdb.connect("localhost", "root", "Coufalka*1", "NGS_FORENSIC")
c = db.cursor()
for sample in sample_names:
    if 'AutoSTR' in table_list:
        select_head_id = "SELECT id FROM ngs_forensic_test.Heads_flankingReg WHERE sample_name = '%s'" % (sample)
        c.execute(select_head_id)
        row_count = c.rowcount
        if row_count == 0:
            print (sample, ' is not in AutoSTR database')
        else:
            sql_delete_Query = "DELETE FROM ngs_forensic_test.heads_flankingreg where sample_name = '%s'" % (sample)
            c.execute(sql_delete_Query)
            db.commit()
            print(sample + ' was deleted from AutoSTR database')
  
    if 'Y-STR' in table_list:
        select_head_id_Y = "SELECT id FROM ngs_forensic_test.Heads_y_flankingReg WHERE sample_name = '%s'" % (sample)
        c.execute(select_head_id)
        row_count = c.rowcount
        if row_count == 0:
            print (sample, ' is not in Y-STR database')
        else:
            sql_delete_Query_Y = "DELETE FROM ngs_forensic_test.heads_y_flankingreg where sample_name = '%s'" % (sample)
            c.execute(sql_delete_Query_Y)
            db.commit()
            print(sample + ' was deleted from Y_STR database')
        
db.close()  
print('done')
    