import os
import sys
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
#import MySQLdb
#import mysql.connector as MySQLdb
import math
from datetime import datetime
import sys

def systemLinux():
    if sys.platform == 'linux':
        return True

if systemLinux():
    import mysql.connector as MySQLdb
else:
    import MySQLdb

def main():
    path_samples_delete = os.path.normpath('C:/NGS_forensic_database/delete_samples/delete_sample.txt')
    table_list = ['AutoSTR', 'Y-STR']
    dbPass = 'XXX'
    delete_samples(path_samples_delete, table_list, dbPass)

def delete_samples(path_samples_delete, table_list, dbPass):
    sample_names = []
    
    ### *** get samples from file
    for sample in csv.reader(open(path_samples_delete), delimiter=';'):
        if len(sample) == 1:
            sample_names.append(sample[0])
        else:
            print ('ERROR: wrong format ' + path_samples_print)
    
    ### *** delete data from from mySQL dtb  ***
    if systemLinux():
        db = MySQLdb.connect(user="root", password=dbPass, database="ngs_forensic")
        c = db.cursor(buffered=True)
    else:
        db = MySQLdb.connect("localhost", "root", dbPass, "NGS_FORENSIC")
        c = db.cursor()

    for sample in sample_names:
        if 'AutoSTR' in table_list:
            select_head_id = "SELECT id FROM ngs_forensic.heads_flankingreg WHERE sample_name = '%s'" % (sample)
            c.execute(select_head_id)
            row_count = c.rowcount
            if row_count == 0:
                print (sample, ' is not in AutoSTR database')
            else:
                sql_delete_Query = "DELETE FROM ngs_forensic.heads_flankingreg where sample_name = '%s'" % (sample)
                c.execute(sql_delete_Query)
                db.commit()
                print(sample + ' was deleted from AutoSTR database')
    
        if 'Y-STR' in table_list:
            select_head_id_Y = "SELECT id FROM ngs_forensic.heads_y_flankingreg WHERE sample_name = '%s'" % (sample)
            c.execute(select_head_id_Y)
            row_count = c.rowcount
            if row_count == 0:
                print (sample, ' is not in Y-STR database')
            else:
                sql_delete_Query_Y = "DELETE FROM ngs_forensic.heads_y_flankingreg where sample_name = '%s'" % (sample)
                c.execute(sql_delete_Query_Y)
                db.commit()
                print(sample + ' was deleted from Y_STR database')
            
    db.close()  
    print('done')

if __name__ == '__main__':
    main()