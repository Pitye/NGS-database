import os
import csv
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
    query_table_dict = {}
    ### *** get samples from file
    for sample in csv.reader(open(path_samples_delete), delimiter=';'):
        if len(sample) == 1:
            sample_names.append(sample[0])
        else:
            print ('ERROR: wrong format ' + path_samples_delete)
    
    ### *** delete data from from mySQL dtb  ***
    if systemLinux():
        db = MySQLdb.connect(user="root", password=dbPass, database="ngs_forensic")
        c = db.cursor(buffered=True)
    else:
        db = MySQLdb.connect("localhost", "root", dbPass, "NGS_FORENSIC")
        c = db.cursor()

    if 'AutoSTR' in table_list:
        query_table_dict['heads_flankingreg'] = 'AutoSTR'

    if 'AutoSTR_Family-tree' in table_list:
        query_table_dict['heads_family_tree'] = 'AutoSTR_Family-tree'

    if 'Y-STR' in table_list:
        query_table_dict['heads_y_flankingreg'] = 'Y-STR'

    for table in query_table_dict.keys():
        for sample in sample_names:
            select_head_id = "SELECT id FROM ngs_forensic.%s WHERE sample_name = '%s'" % (table, sample)
            c.execute(select_head_id)
            row_count = c.rowcount
            if row_count == 0:
                print(sample, ' is not in ' + query_table_dict[table] + ' database')
            else:
                sql_delete_Query = "DELETE FROM ngs_forensic.%s WHERE sample_name = '%s'" % (table, sample)
                c.execute(sql_delete_Query)
                db.commit()
                print(sample + ' was deleted ' + query_table_dict[table] + ' database')
            
    db.close()  
    print('done')

if __name__ == '__main__':
    main()