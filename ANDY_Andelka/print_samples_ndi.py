import os
import csv

#import MySQLdb
#import mysql.connector as MySQLdb

from datetime import datetime
import xlwt

import sys

def systemLinux():
    if sys.platform == 'linux':
        return True

if systemLinux():
    import mysql.connector as MySQLdb
else:
    import MySQLdb

def main():
    directory_reports = os.path.normpath('C:/NGS_forensic_database/print_samples/prints')
    path_samples_print = os.path.normpath('C:/NGS_forensic_database/print_samples/print_sample.txt')
    path_locus_order = os.path.normpath('C:/NGS_forensic_database/print_samples/locus_order.txt')
    table_list = ['AutoSTR', 'AutoSTR_Family-tree', 'Y-STR']
    selected_columns = ['allele', 'PubMed_ID','no_reads' ]
    dbPass = 'XXX'
    print_samples (path_samples_print, path_locus_order, directory_reports, table_list, selected_columns, dbPass)

def print_samples (path_samples_print, path_locus_order, directory_reports, table_list, selected_columns, dbPass):
    columnsAuto = ['allele', 'seq_name', 'PubMed_ID', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']
    columnsY = ['allele', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']
    sample_names = []
    locus_list = []
    query_table_list = []
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

    sample_records = {sample: {marker: {column: [] for column in columnsAuto} for marker in markers} for sample in sample_names}

    ### *** get data from from mySQL dtb - sample_records ***
    if systemLinux():
        db = MySQLdb.connect(user="root", password=dbPass, database="ngs_forensic")
        c = db.cursor(buffered=True)
    else:
        db = MySQLdb.connect("localhost", "root", dbPass, "NGS_FORENSIC")
        c = db.cursor()

    if 'AutoSTR' in table_list:
        query_table_list.append('nomen_freq_autostrdata_flankingreg')

    if 'AutoSTR_Family-tree' in table_list:
        query_table_list.append('nomen_freq_autostrdata_family_tree')

    if 'Y-STR' in table_list:
        query_table_list.append('freq_y_strdata_flankingreg')



    for sample in sample_names:
        for table in query_table_list:
            sql_select_Query = "SELECT * FROM ngs_forensic.%s where sample_name = '%s'" % (table, sample)
            c.execute(sql_select_Query)
            records = c.fetchall()

            for row in records:
            #marker, row[1], allele, row[2], seq_name, row[3], PubMed_ID, row[4], sequence, row[5], no_reads, row[6], CE_validation, row[7], head_id, row[8], avg_no_reads, row[9], count_seq, row[10], frequency, row[11]
                column_index = 2
                if table in ['nomen_freq_autostrdata_flankingreg', 'nomen_freq_autostrdata_family_tree']:
                    columns = columnsAuto
                if table in ['freq_y_strdata_flankingreg']:
                    columns = columnsY
                for column in columns:
                    sample_records[sample][row[1]][column].append(row[column_index])
                    column_index += 1


    db.close()
    #print (sample_records)
    for sample in sample_names:
        for marker in markers:
            for column in columnsAuto:
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
    report_path_name_xls = directory_reports + os.path.normpath('/') + 'report_print_' + dt_string + '.xls'
    f = open (report_path_name, 'w+')
    f.writelines(["*** ANDY - NGS data interpreter - REPORT ***\n\n"])
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
    wb = xlwt.Workbook()
    sh = wb.add_sheet('report')
    with open(report_path_name, 'r') as fc:
        reader = csv.reader(fc)
        for r, row in enumerate(reader):
            for c, val in enumerate(row):
                sh.write(r, c, val)
    wb.save(report_path_name_xls)

    print ('done')

if __name__ == '__main__':
    main()
