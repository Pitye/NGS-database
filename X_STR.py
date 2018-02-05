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
sheets = ['Autosomal STRs', 'Y STRs', 'X STRs', 'iSNPs']
no_reads_for_validation_X = 150



#reading CSV files in directory for sheet[1] - 'X STRs'
csv_list_2 = os.listdir(out_directory + sheets[2] + '/')
#
#
markers_heteroyzgozity_dict_X = {'DXS10135': 0.5,'DXS8378' : 0.5,'DXS7132' : 0.5,'DXS10074' : 0.5,'DXS10103' : 0.5,'HPRTB' : 0.5,'DXS7423' : 0.5}
#
markers_X = list(markers_heteroyzgozity_dict_X.keys())

X_STR_Data = {}
X_STR_Head = {}

for file in csv_list_2:
    
    writing_row_marker_X = 0
    head_raw_X = []
    row_number_X = 0
    
    #create empty raw_dictionary_X and dictionary with homo or hetero alleles for markers
    sample_raw_dict_X = {marker_X : [] for marker_X in markers_X}
    sample_dict_X = {marker_X : [] for marker_X in markers_X}
    sorted_sample_dict_X = {marker_X : [] for marker_X in markers_X}
    
    #reading csv files row by row
    for row_X in csv.reader(open(out_directory + sheets[2] + '/'+ file), delimiter=','):
        
        row_number_X += 1
        
        if row_X [0] == 'Locus':
            writing_row_marker_X += 1
        
        if writing_row_marker_X == 0:
            head_raw_X.append(row_X)
        
        #writing to raw_dictionary for markers
        if writing_row_marker_X == 2:
            if row_X[0] in markers_X:
                sample_raw_dict_X[row_X[0]].append(row_X)
    
    #if gender is XY
    if head_raw_X[6][1] == 'XY':
    
        #select true readings to sample_dict according to marker's heterozygozity              
        for marker_X in markers_X:
            
            sorted_markers_X = sorted(sample_raw_dict_X[marker_X], key = lambda x: int(x[3]),reverse = True)
        
            if len (sorted_markers_X) == 1:
                sample_dict_X[marker_X].append(sorted_markers_X[0])
                
            
            if len (sorted_markers_X) > 1:
                sample_dict_X[marker_X].append(sorted_markers_X[0])
                
                if int(sorted_markers_X[1][3])/int(sorted_markers_X[0][3]) > float(markers_heteroyzgozity_dict_X[marker_X]):
                    sample_dict_X[marker_X].append(sorted_markers_X[1])
                    
                    if len (sorted_markers_X) > 2:
                        if int(sorted_markers_X[2][3])/int(sorted_markers_X[1][3]) > float(markers_heteroyzgozity_dict_X[marker_X]):
                            sample_dict_X[marker_X].append(sorted_markers_X[2])
                
               
            if len (sorted_markers_X) == 0:
                sample_dict_X[marker_X].append(['Null', '0', 'Null', '0', 'Null'])
                
                
            sorted_sample_dict_X[marker_X] = sorted(sample_dict_X[marker_X], key = lambda x: float(x[1]))
    else:
        for marker_X in markers_X:
            
            sorted_markers_X = sorted(sample_raw_dict_X[marker_X], key = lambda x: int(x[3]),reverse = True)
        
            if len (sorted_markers_X) == 1:
                sample_dict_X[marker_X].append(sorted_markers_X[0])
                sample_dict_X[marker_X].append(sorted_markers_X[0])
            
            if len (sorted_markers_X) > 1:
                if int(sorted_markers_X[1][3])/int(sorted_markers_X[0][3]) > float(markers_heteroyzgozity_dict_X[marker_X]):
                    sample_dict_X[marker_X].append(sorted_markers_X[0])
                    sample_dict_X[marker_X].append(sorted_markers_X[1])
                
                else: 
                    sample_dict_X[marker_X].append(sorted_markers_X[0])
                    sample_dict_X[marker_X].append(sorted_markers_X[0])
                
            if len (sorted_markers_X) == 0:
                sample_dict_X[marker_X].append(['Null', '0', 'Null', '0', 'Null'])
                sample_dict_X[marker_X].append(['Null', '0', 'Null', '0', 'Null'])
                
        
            sorted_sample_dict_X[marker_X] = sorted(sample_dict_X[marker_X], key = lambda x: float(x[1]))
        
    
    
    head_dict_X = {'Project' : head_raw_X[3][1], 'Analysis': head_raw_X[4][1], 'Run' : head_raw_X[5][1], 'Gender' : head_raw_X[6][1], 'Created' : head_raw_X[7][1]}

    sample_name_X = head_raw_X[2][1]
    
    # X_STR_Data is dictionary in formate {sample_name : {marker : [[marker, allele1, Yes/No, no.readings, seq.], [marker, allele2, Yes/No, no.readings, seq.]]} }
    X_STR_Data[sample_name_X] = sorted_sample_dict_X
    X_STR_Head[sample_name_X] = head_dict_X
    
        #if gender is XX
    
    
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
            
print('Data_CE done') 


samples_NGS_X = (list(X_STR_Data.keys()))
samples_CE = (list(Data_CE.keys()))
genotype_NGS_X = []

#validating 
for sample_NGS_X in samples_NGS_X:
    no_mismatches_X = 0
    
    if sample_NGS_X in samples_CE:
        
        for marker_X in markers_X:
            no_alleles_X = len(X_STR_Data[sample_NGS_X][marker_X])
            if marker_X not in list(Data_CE[sample_NGS_X].keys()):
                
                for i in range (no_alleles_X):
                    if int(X_STR_Data[sample_NGS_X][marker_X][i][3]) >= no_reads_for_validation_X:
                        X_STR_Data[sample_NGS_X][marker_X][i] = X_STR_Data[sample_NGS_X][marker_X][i] + ['validated_no_reads']
              
                    else:
                        X_STR_Data[sample_NGS_X][marker_X][i] = X_STR_Data[sample_NGS_X][marker_X][i] + ['not_validated_CE']
                    
            
            if marker_X in list(Data_CE[sample_NGS_X].keys()):
                genotype_NGS_X = []
                for i in range (no_alleles_X):
                    genotype_NGS_X = genotype_NGS_X + [str(Y_STR_Data[sample_NGS_X][marker_X][i][1])]
                
                genotype_CE = Data_CE [sample_NGS_X][marker_X]
                #print (sample_NGS_Y, marker_Y, genotype_CE, genotype_NGS_Y )               
                if genotype_CE == genotype_NGS_X:
                    for i in range (no_alleles_X):
                        X_STR_Data[sample_NGS_X][marker_X][i] = X_STR_Data[sample_NGS_X][marker_X][i] + ['validated_CE']
                    
                    
                else:
                    for i in range (no_alleles_X):
                        X_STR_Data[sample_NGS_X][marker_X][i] = X_STR_Data[sample_NGS_X][marker_X][i] + ['locus_mismatch_CE']
                    
                    no_mismatches_X = no_mismatches_X + 1
        
            X_STR_Head[sample_NGS_X]['No_mismatches_X']= no_mismatches_X
    else:
        print(sample_NGS_X, ' is not in xml files')
        X_STR_Head[sample_NGS_X]['No_mismatches_X']= no_mismatches_X
        for marker_X in markers_X:
            no_alleles_X = len(X_STR_Data[sample_NGS_X][marker_X])
            for i in range (no_alleles_X):
                if int(X_STR_Data[sample_NGS_X][marker_X][i][3]) >= no_reads_for_validation_X:
                    X_STR_Data[sample_NGS_X][marker_X][i] = X_STR_Data[sample_NGS_X][marker_X][i] + ['validated_no_reads']
              
                else:
                    X_STR_Data[sample_NGS_X][marker_X][i] = X_STR_Data[sample_NGS_X][marker_X][i] + ['not_validated_CE']
                    
       
        
print('X_STR_Head done')
print('X_STR_Data done') 
#print(Y_STR_Head)
#print(Y_STR_Data) 

#insert data from 'Auto_STR_Head' and 'Auto_STR_Data' to MySQL NGS_FORENSIC database (create dtbschema by querries in sql file)'
db=MySQLdb.connect("localhost", "root", "coufalka", "NGS_FORENSIC")
c = db.cursor()
sample_names_X = (list(X_STR_Data.keys()))
for sample_name_X in sample_names_X:
    pr = X_STR_Head[sample_name_X]['Project']
    an = X_STR_Head[sample_name_X]['Analysis']
    ru = X_STR_Head[sample_name_X]['Run']
    ge = X_STR_Head[sample_name_X]['Gender']
    cr = X_STR_Head[sample_name_X]['Created']
    nm = int(X_STR_Head[sample_name_X]['No_mismatches_X'])
    
    insert_head_X = "INSERT INTO Heads_X (sample_name, project, analysis, run, gender, created, no_mismatches_X) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%d')" % (sample_name_X, pr, an, ru, ge, cr, nm) 
    c.execute(insert_head_X)
    db.commit()
    print(sample_name_X, " head_X done")
    select_head_X_id = "SELECT id FROM Heads_X WHERE sample_name = '%s' AND \
                                                 project = '%s' AND \
                                                 analysis = '%s' AND \
                                                 run = '%s' AND \
                                                 gender = '%s' AND \
                                                 created = '%s' AND \
                                                 no_mismatches_X = '%d' " % (sample_name_X, pr, an, ru, ge, cr, nm)
    c.execute(select_head_X_id)
    head_X_id = int(c.fetchone()[0])
    print (head_X_id)
    for marker_X in markers_X:
        no_alleles_X = len(X_STR_Data[sample_name_X][marker_X])
        for x in range (no_alleles_X):
            al = X_STR_Data[sample_name_X][marker_X][x][1]
            nr = int(X_STR_Data[sample_name_X][marker_X][x][3])
            seq = X_STR_Data[sample_name_X][marker_X][x][4]
            val = X_STR_Data[sample_name_X][marker_X][x][5]
            insert_X_STR = "INSERT INTO X_STRdata (sample_name, marker, allele, sequence, no_reads, CE_validation, head_X_id) \
            VALUES ('%s', '%s', '%s', '%s', '%d', '%s', '%d')" % (sample_name_X, marker_X, al, seq, nr, val, head_X_id)
            c.execute(insert_X_STR)
            db.commit()
            print (sample_name_X, ' ', marker_X, ' done')

            
   
            

db.close()
print('all done')    


