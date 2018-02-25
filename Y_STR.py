import pandas as pd
import os
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
import MySQLdb

in_directory = '/Users/sweethome/Desktop/project_NGS/NGS/xlsx/'
out_directory = '/Users/sweethome/Desktop/project_NGS/NGS/csv/'
xml_directory = '/Users/sweethome/Desktop/project_NGS/NGS/xml_CE/'
CSV_CE_directory = '/Users/sweethome/Desktop/project_NGS/NGS/csv_CE/'
sheets = ['Autosomal STRs', 'Y STRs', 'X STRs', 'iSNPs']
no_reads_for_validation_Y = 150

#reading data from XML files to dictionary Data_CE for validating Auto_STR_Data
def read_XML_file(path):    
    tree = etree.parse(path)
    namespaces = { 'ns': 'urn:CODISImportFile-schema' }
    root = tree.getroot()

    XML_data = {}
    for specimen_element in root.findall("ns:SPECIMEN", namespaces):
        specimen_id = specimen_element.find("ns:SPECIMENID", namespaces).text
        genotype = {}
        for locus_element in specimen_element.findall("ns:LOCUS", namespaces):
            locus_name = locus_element.find("ns:LOCUSNAME", namespaces).text.replace(' ', '')
            allele_values = []

            for allele_value_element in locus_element.findall("*/ns:ALLELEVALUE", namespaces):
                allele_values.append(allele_value_element.text)
                
            if not locus_name.startswith('DYS') and not locus_name.startswith('YGAT') and len(allele_values) == 1:
                allele_values.append(allele_values[0])

            genotype[locus_name] = allele_values
            #dictionary[new_key] = dictionary.pop(old_key)
            if 'YGATAH4'in locus_name:
                genotype['Y-GATA-H4'] = genotype.pop('YGATAH4')
            if 'DYS385'in locus_name:    
                genotype['DYS385a-b'] = genotype.pop('DYS385')
            
            
        XML_data[specimen_id] = genotype

    return XML_data   

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


#Data_CE is dictionary in formate {sample_name : {markerA : [allele1, allele2], markerB: [allele1, allele2]}}
xml_list = os.listdir(xml_directory)
Data_CE = {}
for file in xml_list:
    if file.endswith('.xml'):
        path_xml = xml_directory + file
        Data_CE_part = read_XML_file(path_xml)   
        #print(list(Data_CE_part.keys()))
        Data_CE = merge_two_dicts(Data_CE, Data_CE_part)
            

#reading data from CSV_CE files to dictionary Data_CE for validating Auto_STR_Data
CSV_CE_list = os.listdir(CSV_CE_directory)
CSV_data_CE = {}
for file in CSV_CE_list:
    if file.endswith('.csv'):
        genotype = {}

        for row in csv.reader(open(CSV_CE_directory + file), delimiter=','):    
            specimen_id = row[9]
            locus_name = row[79].replace(' ', '')
            allele_values = row[80].split(", ")
            if not locus_name.startswith('DYS') and not locus_name.startswith('YGAT') and len(allele_values) == 1:
                allele_values.append(allele_values[0])
            genotype[locus_name] = allele_values
   
            
        CSV_data_CE[specimen_id] = genotype


Data_CE = merge_two_dicts(Data_CE, CSV_data_CE)            
print('Data_CE done') 

#reading CSV files in directory for sheet[1] - 'Y STRs'
csv_list_1 = os.listdir(out_directory + sheets[1] + '/')
#
#
markers_heteroyzgozity_dict_Y = {'DYS505': 0.5,'DYS570' : 0.5,'DYS576' : 0.5,'DYS522' : 0.5,'DYS481' : 0.5,'DYS19' : 0.5,'DYS391' : 0.5,'DYS635' : 0.5,'DYS437' : 0.5,'DYS439' : 0.5,\
                                 'DYS389I' : 0.5,'DYS389II' : 0.5,'DYS438' : 0.5,'DYS612' : 0.5,'DYS390': 0.5,'DYS643' : 0.5,'DYS533' : 0.5,'Y-GATA-H4' : 0.5,'DYS385a-b' : 0.5,'DYS460' : 0.5,\
                                 'DYS549' : 0.5,'DYS392' : 0.5,'DYS448' : 0.5,'DYF387S1' : 0.5}
#
markers_Y = list(markers_heteroyzgozity_dict_Y.keys())

Y_STR_Data = {}
Y_STR_Head = {}

for file in csv_list_1:
    if file.endswith('.csv'):
        writing_row_marker_Y = 0
        head_raw_Y = []
        row_number_Y = 0
    
    #create empty raw_dictionary_Y and dictionary with homo or hetero alleles for markers
        sample_raw_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
        sample_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
        sorted_sample_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
    
    #reading csv files row by row
        for row_Y in csv.reader(open(out_directory + sheets[1] + '/'+ file), delimiter=','):
        
            row_number_Y += 1
        
            if row_Y [0] == 'Locus':
                writing_row_marker_Y += 1
        
            if writing_row_marker_Y == 0:
                head_raw_Y.append(row_Y)
        
        #writing to raw_dictionary for markers
            if writing_row_marker_Y == 2:
                if row_Y[0] in markers_Y:
                    sample_raw_dict_Y[row_Y[0]].append(row_Y)
    
    #if gender is XY
        if head_raw_Y[6][1] == 'XY':
    
        #select true readings to sample_dict according to marker's heterozygozity              
            for marker_Y in markers_Y:
            
                sorted_markers_Y = sorted(sample_raw_dict_Y[marker_Y], key = lambda x: int(x[3]),reverse = True)
        
                if len (sorted_markers_Y) == 1:
                    sample_dict_Y[marker_Y].append(sorted_markers_Y[0])
                
            
                if len (sorted_markers_Y) > 1:
                    sample_dict_Y[marker_Y].append(sorted_markers_Y[0])
                
                    if int(sorted_markers_Y[1][3])/int(sorted_markers_Y[0][3]) > float(markers_heteroyzgozity_dict_Y[marker_Y]):
                        sample_dict_Y[marker_Y].append(sorted_markers_Y[1])
                    
                        if len (sorted_markers_Y) > 2:
                            if int(sorted_markers_Y[2][3])/int(sorted_markers_Y[1][3]) > float(markers_heteroyzgozity_dict_Y[marker_Y]):
                                sample_dict_Y[marker_Y].append(sorted_markers_Y[2])
                
               
                if len (sorted_markers_Y) == 0:
                    sample_dict_Y[marker_Y].append(['Null', '0', 'Null', '0', 'Null'])
                
                
                sorted_sample_dict_Y[marker_Y] = sorted(sample_dict_Y[marker_Y], key = lambda x: float(x[1]))

        
            head_dict_Y = {'Project' : head_raw_Y[3][1], 'Analysis': head_raw_Y[4][1], 'Run' : head_raw_Y[5][1], 'Gender' : head_raw_Y[6][1], 'Created' : head_raw_Y[7][1]}

            sample_name_Y = head_raw_Y[2][1]
    
        # Auto_STR_Data is dictionary in formate {sample_name : {marker : [[marker, allele1, Yes/No, no.readings, seq.], [marker, allele2, Yes/No, no.readings, seq.]]} }
            Y_STR_Data[sample_name_Y] = sorted_sample_dict_Y
            Y_STR_Head[sample_name_Y] = head_dict_Y
    



samples_NGS_Y = (list(Y_STR_Data.keys()))
samples_CE = (list(Data_CE.keys()))
genotype_NGS_Y = []

#validating 
for sample_NGS_Y in samples_NGS_Y:
    no_mismatches_Y = 0
    
    if sample_NGS_Y in samples_CE:
        
        for marker_Y in markers_Y:
            no_alleles_Y = len(Y_STR_Data[sample_NGS_Y][marker_Y])
            if marker_Y not in list(Data_CE[sample_NGS_Y].keys()):
                
                for i in range (no_alleles_Y):
                    if int(Y_STR_Data[sample_NGS_Y][marker_Y][i][3]) >= no_reads_for_validation_Y:
                        Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['validated_no_reads']
              
                    else:
                        Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['not_validated_CE']
                    
            
            if marker_Y in list(Data_CE[sample_NGS_Y].keys()):
                genotype_NGS_Y = []
                for i in range (no_alleles_Y):
                    genotype_NGS_Y = genotype_NGS_Y + [str(Y_STR_Data[sample_NGS_Y][marker_Y][i][1])]
                
                genotype_CE = Data_CE [sample_NGS_Y][marker_Y]
                #print (sample_NGS_Y, marker_Y, genotype_CE, genotype_NGS_Y )               
                if genotype_CE == genotype_NGS_Y:
                    for i in range (no_alleles_Y):
                        Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['validated_CE']
                    
                    
                else:
                    for i in range (no_alleles_Y):
                        Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['locus_mismatch_CE']
                    
                    no_mismatches_Y = no_mismatches_Y + 1
        
            Y_STR_Head[sample_NGS_Y]['No_mismatches_Y']= no_mismatches_Y
    else:
        print(sample_NGS_Y, ' is not in xml files')
        Y_STR_Head[sample_NGS_Y]['No_mismatches_Y']= no_mismatches_Y
        for marker_Y in markers_Y:
            no_alleles_Y = len(Y_STR_Data[sample_NGS_Y][marker_Y])
            for i in range (no_alleles_Y):
                if int(Y_STR_Data[sample_NGS_Y][marker_Y][i][3]) >= no_reads_for_validation_Y:
                    Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['validated_no_reads']
              
                else:
                    Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['not_validated_CE']
                    
       
        
print('Y_STR_Head done')
print('Y_STR_Data done') 
#print(Y_STR_Head)
#print(Y_STR_Data) 

#insert data from 'Auto_STR_Head' and 'Auto_STR_Data' to MySQL NGS_FORENSIC database (create dtbschema by querries in sql file)'
db=MySQLdb.connect("localhost", "root", "coufalka", "NGS_FORENSIC")
c = db.cursor()
sample_names_Y = (list(Y_STR_Data.keys()))
for sample_name_Y in sample_names_Y:
    pr = Y_STR_Head[sample_name_Y]['Project']
    an = Y_STR_Head[sample_name_Y]['Analysis']
    ru = Y_STR_Head[sample_name_Y]['Run']
    ge = Y_STR_Head[sample_name_Y]['Gender']
    cr = Y_STR_Head[sample_name_Y]['Created']
    nm = int(Y_STR_Head[sample_name_Y]['No_mismatches_Y'])
    
    insert_head_Y = "INSERT INTO Heads_Y (sample_name, project, analysis, run, gender, created, no_mismatches_Y) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%d')" % (sample_name_Y, pr, an, ru, ge, cr, nm) 
    c.execute(insert_head_Y)
    db.commit()
    print(sample_name_Y, " head_Y done")
    select_head_Y_id = "SELECT id FROM Heads_Y WHERE sample_name = '%s' AND \
                                                 project = '%s' AND \
                                                 analysis = '%s' AND \
                                                 run = '%s' AND \
                                                 gender = '%s' AND \
                                                 created = '%s' AND \
                                                 no_mismatches_Y = '%d' " % (sample_name_Y, pr, an, ru, ge, cr, nm)
    c.execute(select_head_Y_id)
    head_Y_id = int(c.fetchone()[0])
    print (head_Y_id)
    for marker_Y in markers_Y:
        no_alleles_Y = len(Y_STR_Data[sample_name_Y][marker_Y])
        for x in range (no_alleles_Y):
            al = Y_STR_Data[sample_name_Y][marker_Y][x][1]
            nr = int(Y_STR_Data[sample_name_Y][marker_Y][x][3])
            seq = Y_STR_Data[sample_name_Y][marker_Y][x][4]
            val = Y_STR_Data[sample_name_Y][marker_Y][x][5]
            insert_Y_STR = "INSERT INTO Y_STRdata (sample_name, marker, allele, sequence, no_reads, CE_validation, head_Y_id) \
            VALUES ('%s', '%s', '%s', '%s', '%d', '%s', '%d')" % (sample_name_Y, marker_Y, al, seq, nr, val, head_Y_id)
            c.execute(insert_Y_STR)
            db.commit()
            print (sample_name_Y, ' ', marker_Y, ' done')

            
   
            

db.close()
print('all done')    

