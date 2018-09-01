import pandas as pd
import os
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
import pymysql
pymysql.install_as_MySQLdb()
import MySQLdb

in_directory = os.path.normpath('C:/NGS_forensic_database/xlsx_detail_reports') 
out_directory = os.path.normpath('C:/NGS_forensic_database/csv_output') 
xml_directory = os.path.normpath('C:/NGS_forensic_database/xml_CE') 
CSV_CE_directory = os.path.normpath('C:/NGS_forensic_database/csv_CE') 
sheets = ['Autosomal STRs', 'Y STRs', 'X STRs', 'iSNPs']
no_reads_for_validation_i = 30

#reading CSV files in directory for sheet[3] - 'iSNPs'
csv_list_3 = os.listdir(out_directory + os.path.normpath('/') + sheets[3] +  os.path.normpath('/'))

markers_heteroyzgozity_dict_i = {'rs1490413': 0.5,'rs560681' : 0.5,'rs1294331' : 0.5,'rs10495407' : 0.5,'rs891700' : 0.5,'rs1413212' : 0.5,'rs876724' : 0.5,'rs1109037' : 0.5,'rs993934' : 0.5,'rs12997453' : 0.5, \
                                 'rs907100' : 0.5,'rs1357617' : 0.5,'rs4364205' : 0.5,'rs2399332' : 0.5,'rs1355366': 0.5,'rs6444724' : 0.5,'rs2046361' : 0.5,'rs279844' : 0.5,'rs6811238' : 0.5,'rs1979255' : 0.5,\
                                 'rs717302' : 0.5,'rs159606' : 0.5,'rs13182883' : 0.5,'rs251934' : 0.5,'rs338882' : 0.5,'rs13218440' : 0.5,'rs1336071' : 0.5, 'rs214955': 0.5,'rs727811' : 0.5,'rs6955448' : 0.5,\
                                 'rs917118' : 0.5,'rs321198' : 0.5,'rs737681' : 0.5,'rs763869' : 0.5,'rs10092491' : 0.5,'rs2056277' : 0.5,'rs4606077' : 0.5,'rs1015250' : 0.5,'rs7041158' : 0.5,'rs1463729' : 0.5,\
                                 'rs1360288': 0.5,'rs10776839' : 0.5,'rs826472' : 0.5,'rs735155' : 0.5,'rs3780962' : 0.5,'rs740598' : 0.5,'rs964681' : 0.5,'rs1498553' : 0.5,'rs901398' : 0.5,'rs10488710' : 0.5,\
                                 'rs2076848' : 0.5,'rs2107612' : 0.5,'rs2269355' : 0.5,'rs2920816': 0.5,'rs2111980' : 0.5,'rs10773760' : 0.5,'rs1335873' : 0.5,'rs1886510' : 0.5,'rs1058083' : 0.5,'rs354439' : 0.5,\
                                 'rs1454361' : 0.5,'rs722290' : 0.5,'rs873196' : 0.5,'rs4530059' : 0.5,'rs1821380' : 0.5,'rs8037429' : 0.5,'rs1528460' : 0.5,'rs729172': 0.5,'rs2342747' : 0.5,'rs430046' : 0.5,\
                                 'rs1382387' : 0.5,'rs9905977' : 0.5,'rs740910' : 0.5,'rs938283' : 0.5,'rs8078417' : 0.5,'rs1493232' : 0.5,'rs9951171' : 0.5,'rs1736442' : 0.5,'rs1024116' : 0.5,'rs719366' : 0.5,\
                                 'rs576261': 0.5,'rs1031825' : 0.5,'rs445251' : 0.5,'rs1005533' : 0.5,'rs1523537' : 0.5,'rs722098' : 0.5,'rs2830795' : 0.5,'rs2831700' : 0.5,'rs914165' : 0.5,'rs221956' : 0.5,\
                                 'rs733164' : 0.5,'rs987640' : 0.5,'rs2040411' : 0.5,'rs1028528' : 0.5}
#

markers_i = list(markers_heteroyzgozity_dict_i.keys())

i_SNP_Data = {}


for file in csv_list_3:
    
    writing_row_marker = 0
    head_raw = []
    row_number = 0
    
    #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
    sample_raw_dict_i = {marker : [] for marker in markers_i}
    sample_dict_i = {marker : [] for marker in markers_i}
    sorted_sample_dict_i = {marker : [] for marker in markers_i}
    
    #reading csv files row by row
    for row in csv.reader(open(out_directory +  os.path.normpath('/') + sheets[3] +  os.path.normpath('/') + file), delimiter=','):
        
        row_number += 1
        
        if row [0] == 'Locus':
            writing_row_marker += 1
        
        if writing_row_marker == 0:
            head_raw.append(row)
        
        #writing to raw_dictionary for markers
        if writing_row_marker == 2:
            if row[0] in markers_i:
                sample_raw_dict_i[row[0]].append(row)
    
    #select true readings to sample_dict according to marker's heterozygozity              
    for marker in markers_i:
            
        sorted_markers_i = sorted(sample_raw_dict_i[marker], key = lambda x: int(x[3]),reverse = True)
        
        if len (sorted_markers_i) == 1:
            sample_dict_i[marker].append(sorted_markers_i[0])
            sample_dict_i[marker].append(sorted_markers_i[0])
            
        if len (sorted_markers_i) > 1 and int(sorted_markers_i[0][3]) > 0 :
            if int(sorted_markers_i[1][3])/int(sorted_markers_i[0][3]) > float(markers_heteroyzgozity_dict_i[marker]):
                sample_dict_i[marker].append(sorted_markers_i[0])
                sample_dict_i[marker].append(sorted_markers_i[1])
                
            else: 
                sample_dict_i[marker].append(sorted_markers_i[0])
                sample_dict_i[marker].append(sorted_markers_i[0])
                
        if len (sorted_markers_i) == 0 or int(sorted_markers_i[0][3]) == 0:
            sample_dict_i[marker].append(['Null', '0', 'Null', '0'])
            sample_dict_i[marker].append(['Null', '0', 'Null', '0'])
                
        
        sorted_sample_dict_i[marker] = sorted(sample_dict_i[marker], key = lambda x: x[1])
        

    sample_name = head_raw[2][1]
    
    # Auto_STR_Data is dictionary in formate {sample_name : {marker : [[marker, allele1, Yes/No, no.readings], [marker, allele2, Yes/No, no.readings]]} }
    i_SNP_Data[sample_name] = sorted_sample_dict_i
    #Auto_STR_Head[sample_name] = head_dict
    


samples_NGS_i = (list(i_SNP_Data.keys()))



#validating by no_reads
for sample_NGS_i in samples_NGS_i:
        
    for marker in markers_i:
                
        if int(i_SNP_Data[sample_NGS_i][marker][0][3]) >= no_reads_for_validation_i and int(i_SNP_Data[sample_NGS_i][marker][1][3]) >= no_reads_for_validation_i:
            i_SNP_Data[sample_NGS_i][marker][0] = i_SNP_Data[sample_NGS_i][marker][0] + ['validated_no_reads']
            i_SNP_Data[sample_NGS_i][marker][1] = i_SNP_Data[sample_NGS_i][marker][1] + ['validated_no_reads']
        else:
            i_SNP_Data[sample_NGS_i][marker][0] = i_SNP_Data[sample_NGS_i][marker][0] + ['not_validated_CE']
            i_SNP_Data[sample_NGS_i][marker][1] = i_SNP_Data[sample_NGS_i][marker][1] + ['not_validated_CE']
           
        

#print (i_SNP_Data)        
print('i_SNP_Data_ done')

#insert data from 'i_SNP_Data' to MySQL NGS_FORENSIC database (create dtbschema by querries in sql file)'
db=MySQLdb.connect("localhost", "root", "coufalka", "NGS_FORENSIC")
c = db.cursor()

for sample_name_i in samples_NGS_i:
    
    select_head_id_i = "SELECT id FROM Heads WHERE sample_name = '%s'" % (sample_name_i)
    c.execute(select_head_id_i)
    head_id_i = int(c.fetchone()[0])
    print (sample_name_i, ' inserted in database')
    for marker_i in markers_i:
        for x in range (2):
            al_i = i_SNP_Data[sample_name_i][marker_i][x][1]
            nr_i = int(i_SNP_Data[sample_name_i][marker_i][x][3])
            val_i = i_SNP_Data[sample_name_i][marker_i][x][4]
            insert_i_SNP = "INSERT INTO SNPdata (sample_name, marker, allele, no_reads, validation_by_no_reads, head_id) \
            VALUES ('%s', '%s', '%s', '%d', '%s', '%d')" % (sample_name_i, marker_i, al_i, nr_i, val_i, head_id_i)
            c.execute(insert_i_SNP)
            db.commit()
            #print (sample_name_i, ' ', marker_i, ' done')

            
   
            

db.close()
print('all done')
