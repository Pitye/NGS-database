import pandas as pd
import os
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
#import MySQLdb
import mysql.connector as MySQLdb
from datetime import datetime
import xlwt

def main():

    in_directory = os.path.normpath('C:/NGS_forensic_database/xlsx_detail_reports') 
    out_directory = os.path.normpath('C:/NGS_forensic_database/csv_output') 
    xml_directory = os.path.normpath('C:/NGS_forensic_database/xml_CE') 
    CSV_CE_directory = os.path.normpath('C:/NGS_forensic_database/csv_CE') 
    reports_directory = os.path.normpath('C:/NGS_forensic_database/Y_STR_check_no_mismatches')
    no_reads_for_validation_Y = 150
    CheckIfSampleInDatabase = True
    dbPass = 'XXX'
    
    Y_STR_FR_check_no_mismatches(in_directory, out_directory, xml_directory, CSV_CE_directory, reports_directory, no_reads_for_validation_Y, CheckIfSampleInDatabase, dbPass)

def Y_STR_FR_check_no_mismatches(in_directory, out_directory, xml_directory, CSV_CE_directory, reports_directory, no_reads_for_validation_Y, CheckIfSampleInDatabase, dbPass):
    print('please wait...')
    sheets = ['Autosomal STR Coverage', 'Y STR Coverage', 'X STR Coverage', 'iSNP Coverage']
    #delete all in OUT DIRECTORY
    shutil.rmtree(out_directory)
    os.makedirs(out_directory)
    
    #create OUT DIRECTORIES_A,Y,X,i in CSV directory
    for sheet in sheets:
        os.makedirs(out_directory + os.path.normpath('/') + sheet)
    
    #create CSV files to OUT DIRECTORIES_A,Y,X
    xlsx_list = os.listdir(in_directory)
    for file in xlsx_list:
        if file.endswith('.xlsx'):
            for sheet in sheets:
                path_xlsx = in_directory + os.path.normpath('/') + file
                #path_csv = out_directory + file[0:-5] + '_' + sheet[0] + '.csv'
                path_csv = out_directory + os.path.normpath('/') + sheet + os.path.normpath('/') + file[0:-5] + '_' + sheet[0] + '.csv'
                df = pd.read_excel(path_xlsx, sheet, header=None)
                df.to_csv(path_csv, header=None, index=None)
    print('xlsx2csv done')

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
            path_xml = xml_directory +  os.path.normpath('/') + file
            Data_CE_part = read_XML_file(path_xml)   
            #print(list(Data_CE_part.keys()))
            Data_CE = merge_two_dicts(Data_CE, Data_CE_part)
                

    #reading data from CSV_CE files to dictionary Data_CE for validating Auto_STR_Data
    CSV_CE_list = os.listdir(CSV_CE_directory)
    CSV_data_CE = {}
    for file in CSV_CE_list:
        if file.endswith('.csv'):
            genotype = {}

            for row in csv.reader(open(CSV_CE_directory +  os.path.normpath('/') + file), delimiter=','):    
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
    csv_list_1 = os.listdir(out_directory +  os.path.normpath('/') + sheets[1] +  os.path.normpath('/'))
    #
    #
    markers_heteroyzgozity_dict_Y = {'DYS505': 0.5,'DYS570' : 0.5,'DYS576' : 0.5,'DYS522' : 0.5,'DYS481' : 0.5,'DYS19' : 0.5,'DYS391' : 0.5,'DYS635' : 0.5,'DYS437' : 0.5,'DYS439' : 0.5,\
                                     'DYS389I' : 0.5,'DYS389II' : 0.5,'DYS438' : 0.5,'DYS612' : 0.5,'DYS390': 0.5,'DYS643' : 0.5,'DYS533' : 0.5,'Y-GATA-H4' : 0.5,'DYS385a-b' : 0.5,'DYS460' : 0.5,\
                                     'DYS549' : 0.5,'DYS392' : 0.5,'DYS448' : 0.5,'DYF387S1' : 0.5}
    #
    markers_Y = list(markers_heteroyzgozity_dict_Y.keys())

    Y_STR_Data = {}
    Y_STR_Head = {}
    Y_STR_Data_raw = {}
    Y_STR_Data_unsorted = {}

    for file in csv_list_1:
        if file.endswith('.csv'):
            writing_row_marker_Y = 0
            head_raw_Y = []
            row_number_Y = 0
            project_name_Y = 'null'
            created_name_Y = 'null'
            user_name_Y = 'null'
            sample_name_Y = 'null'
        
            #create empty raw_dictionary_Y and dictionary with homo or hetero alleles for markers
            sample_raw_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
            sample_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
            sorted_sample_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
            row_selected_Y = []
            
            #reading csv files row by row
            for row_Y in csv.reader(open(out_directory +  os.path.normpath('/') + sheets[1] +  os.path.normpath('/')+ file), delimiter=','):
            
                row_number_Y += 1
            
                if row_number_Y == 3 and row_Y [0] == 'Project':
                    project_name_Y = row_Y [1]
                    print (project_name_Y)
                if row_number_Y == 4 and row_Y [0] == 'Created':    
                    created_name_Y = row_Y [1]
                if row_number_Y == 5 and row_Y [0] == 'User': 
                    user_name_Y = row_Y [1]
                if row_number_Y == 13 and row_Y [0] == 'Sample':
                    writing_row_marker_Y = 1
        
                #reading data for samples      
                if row_number_Y > 13 and writing_row_marker_Y == 1 :
                
                    #first sample starts
                    if sample_name_Y != row_Y [0] and sample_name_Y == 'null':
                        Y_STR_Head[row_Y[0]] = {'Project' : project_name_Y, 'Analysis': row_Y[1], 'Run' : 'user '+ user_name_Y, 'Gender' : 'null', 'Created' : created_name_Y}
                    #new sample starts - second, third, ... not first - new empty dictionary is created
                    if sample_name_Y != row_Y [0] and sample_name_Y != 'null':
                        Y_STR_Head[row_Y[0]] = {'Project' : project_name_Y, 'Analysis': row_Y[1], 'Run' : 'user '+ user_name_Y, 'Gender' : 'null', 'Created' : created_name_Y}
                        sample_raw_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
                    
                    #created Y_STR_Data_raw[sample_name_Y] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                    sample_name_Y = row_Y[0]
                    if row_Y[2] in markers_Y:
                        row_selected_Y = [row_Y[2], row_Y[3], 'N/A', row_Y[4], row_Y[5]]
                        sample_raw_dict_Y[row_Y[2]].append(row_selected_Y)
                        Y_STR_Data_raw[sample_name_Y] = sample_raw_dict_Y
            
            
    samples_Y = list(Y_STR_Data_raw.keys())
    for sample_Y in samples_Y:
        
            #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        sample_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
        sorted_sample_dict_Y = {marker_Y : [] for marker_Y in markers_Y}
        

        
        #select true readings to sample_dict according to marker's heterozygozity              
        for marker_Y in markers_Y:
                
            sorted_markers_Y = sorted(Y_STR_Data_raw[sample_Y][marker_Y], key = lambda x: int(x[3]),reverse = True)
            
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

            
        Y_STR_Data[sample_Y] = sorted_sample_dict_Y
        
        # Y_STR_Data is dictionary in formate {sample_name : {marker : [[marker, allele1, Yes/No, no.readings, seq.], [marker, allele2, Yes/No, no.readings, seq.]]} }
                
        

    print ('Y_STR_Head and Y_STR_Data done')

    samples_NGS_Y = (list(Y_STR_Data.keys()))
    samples_CE = (list(Data_CE.keys()))
    genotype_NGS_Y = []
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d%H%M%S")
    report_path_name = reports_directory + os.path.normpath('/') + 'report_Y_STR_check_no_mismatches_' + dt_string + '.csv'
    report_path_name_xls = reports_directory + os.path.normpath('/') + 'report_Y_STR_check_no_mismatches_' + dt_string + '.xls'
    f = open (report_path_name, 'w+')
    f.writelines(["*** ANDY - NGS data interpreter - REPORT ***\n"])
    
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
                    
                    #sorting for comparison
                    genotype_CE.sort()
                    genotype_NGS_Y.sort()
                    
                    if genotype_CE == genotype_NGS_Y:
                        for i in range (no_alleles_Y):
                            Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['validated_CE']
                        
                        
                    else:
                        for i in range (no_alleles_Y):
                            Y_STR_Data[sample_NGS_Y][marker_Y][i] = Y_STR_Data[sample_NGS_Y][marker_Y][i] + ['locus_mismatch_CE']
                            #print (sample_NGS_Y, marker_Y, 'genotype_CE is ', genotype_CE, 'genotype_NGS is', genotype_NGS_Y  )
                        genotype_CE_str = str(genotype_CE).replace('[','').replace(']','').replace('\'','')
                        genotype_NGS_str = str(genotype_NGS_Y).replace('[','').replace(']','').replace('\'','')
                        no_mismatches_Y = no_mismatches_Y + 1
                        f.writelines(["\n" + sample_NGS_Y + "," + marker_Y + "," + "genotype_CE" + "," + genotype_CE_str + "," + "genotype_NGS" + "," + genotype_NGS_str])

            Y_STR_Head[sample_NGS_Y]['No_mismatches_Y']= no_mismatches_Y
            #print(sample_NGS_Y, 'no_mismatches', no_mismatches_Y, "\n")
            f.writelines(["\n" + sample_NGS_Y + "," + "no_mismatches" + "," + str(no_mismatches_Y) +"\n"])
        else:
            #print(sample_NGS_Y, ' is not in xml or csv CE files', "\n")
            f.writelines(["\n" + sample_NGS_Y + "," + " is not in xml or csv CE files" + "\n"])
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
    
    
    
    # check if sample is already in database
    if CheckIfSampleInDatabase:
        f.writelines(["\n\n CHECK IF SAMPLES IN DATABASE"])
        #db=MySQLdb.connect("localhost", "root", dbPass, "NGS_FORENSIC")
        db=MySQLdb.connect(user="root", password=dbPass, database="ngs_forensic")
        #c = db.cursor()
        c = db.cursor(buffered=True)
        sample_names_Y = (list(Y_STR_Data.keys()))
        for sample_name_Y in sample_names_Y:
            select_head_id = "SELECT id FROM heads_y_flankingreg WHERE sample_name = '%s'" % (sample_name_Y)
            c.execute(select_head_id)
            row_count = c.rowcount
            if row_count == 0:
                f.writelines(["\n" + sample_name_Y + "," + " is not in database"])
            else:
                f.writelines(["\n" + sample_name_Y + "," + " entry already exists in database"])
                #print (sample_name_Y, ' entry already exists in Y_database' )
        db.close()
    f.close()
    wb = xlwt.Workbook()
    sh = wb.add_sheet('report')
    with open(report_path_name, 'r') as fc:
        reader = csv.reader(fc)
        for r, row in enumerate(reader):
            for c, val in enumerate(row):
                sh.write(r, c, val)
    wb.save(report_path_name_xls)
    
    print('report created')
    print('all done')
    
if __name__ == '__main__':
    main()
