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
    reports_directory = os.path.normpath('C:/NGS_forensic_database/AutoSTR_import2db')
    no_reads_for_validation = 150
    CheckIfSampleInDatabase = True
    dbPass = 'XXX'
    autoSTR_FR_import2db(in_directory, out_directory, xml_directory, CSV_CE_directory, reports_directory, no_reads_for_validation, CheckIfSampleInDatabase, dbPass)

def autoSTR_FR_import2db(in_directory, out_directory, xml_directory, CSV_CE_directory, reports_directory, no_reads_for_validation, CheckIfSampleInDatabase, dbPass):
    print('please wait...')
    sheets = ['Autosomal STR Coverage', 'X STR Coverage', 'Y STR Coverage', 'iSNP Coverage']

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
    #reading CSV files in directory for sheet[0] - 'Autosomal STRs'
    csv_list_0 = os.listdir(out_directory + os.path.normpath('/') + sheets[0] + os.path.normpath('/'))
    #markers_heteroyzgozity_dict = {'Amelogenin': 0.5,'D1S1656': 0.7,'TPOX' : 0.7,'D2S441' : 0.7,'D2S1338' : 0.7,'D3S1358' : 0.7,'D4S2408' : 0.7,'FGA' : 0.7,'D5S818' : 0.7,'CSF1PO' : 0.7,'D6S1043' : 0.7,'D7S820' : 0.7,'D8S1179' : 0.7,'D9S1122' : 0.7,'D10S1248' : 0.7,'TH01': 0.7,'vWA' : 0.7,'D12S391' : 0.7,'D13S317' : 0.7,'PentaE' : 0.7,'D16S539' : 0.7,'D17S1301' : 0.7,'D18S51' : 0.7,'D19S433' : 0.7,'D20S482' : 0.7,'D21S11' : 0.7,'PentaD' : 0.7,'D22S1045' : 0.7}
    markers_heteroyzgozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5,'TPOX' : 0.5,'D2S441' : 0.5,'D2S1338' : 0.5,'D3S1358' : 0.5,'D4S2408' : 0.5,'FGA' : 0.5,'D5S818' : 0.5,'CSF1PO' : 0.5,'D6S1043' : 0.5,'D7S820' : 0.5,'D8S1179' : 0.5,'D9S1122' : 0.5,'D10S1248' : 0.5,'TH01': 0.5,'vWA' : 0.5,'D12S391' : 0.5,'D13S317' : 0.5,'PentaE' : 0.5,'D16S539' : 0.5,'D17S1301' : 0.5,'D18S51' : 0.5,'D19S433' : 0.5,'D20S482' : 0.5,'D21S11' : 0.5,'PentaD' : 0.5,'D22S1045' : 0.5}
    markers = list(markers_heteroyzgozity_dict.keys())
    
    Auto_STR_Data = {}
    Auto_STR_Head = {}
    Auto_STR_Data_raw = {}
    Auto_STR_Data_unsorted = {}
    for file in csv_list_0:
        writing_row_marker = 0
        head_raw = []
        row_number = 0
        project_name = 'null'
        created_name = 'null'
        user_name = 'null'
        sample_name = 'null'
        #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        sample_raw_dict = {marker : [] for marker in markers}
        sample_dict = {marker : [] for marker in markers}
        sorted_sample_dict = {marker : [] for marker in markers}
        row_selected = []
        
        #reading csv files row by row
        for row in csv.reader(open(out_directory + os.path.normpath('/') + sheets[0] + os.path.normpath('/') + file), delimiter=','):
            row_number += 1
            #print (row_number, '', row [0], row [1])
            if row_number == 3 and row [0] == 'Project':
                project_name = row [1]
                print (project_name)
            if row_number == 4 and row [0] == 'Created':    
                created_name = row [1]
            if row_number == 5 and row [0] == 'User': 
                user_name = row [1]
            if row_number == 13 and row [0] == 'Sample':
                writing_row_marker = 1
            
            #reading data for samples      
            if row_number > 13 and writing_row_marker == 1 :
                
                #first sample starts
                if sample_name != row [0] and sample_name == 'null':
                    Auto_STR_Head[row[0]] = {'Project' : project_name, 'Analysis': row[1], 'Run' : 'user '+ user_name, 'Gender' : 'null', 'Created' : created_name}
                #new sample starts - second, third, ... not first - new empty dictionary is created
                if sample_name != row [0] and sample_name != 'null':
                    Auto_STR_Head[row[0]] = {'Project' : project_name, 'Analysis': row[1], 'Run' : 'user '+ user_name, 'Gender' : 'null', 'Created' : created_name}
                    sample_raw_dict = {marker : [] for marker in markers}
                    
                #created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                sample_name = row[0]
                if row[2] in markers:
                    row_selected = [row[2], row[3], 'N/A', row[4], row[5]]
                    sample_raw_dict[row[2]].append(row_selected)
                    Auto_STR_Data_raw[sample_name] = sample_raw_dict
            
            
    samples = list(Auto_STR_Data_raw.keys())
    for sample in samples:
        
        #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        sample_dict = {marker : [] for marker in markers}
        sorted_sample_dict = {marker : [] for marker in markers}
        
        #select true readings to sample_dict according to marker's heterozygozity 
        for marker in markers:
            sorted_markers = sorted(Auto_STR_Data_raw[sample][marker], key = lambda x: int(x[3]),reverse = True)
            
            if len (sorted_markers) == 1:
                sample_dict[marker].append(sorted_markers[0])
                sample_dict[marker].append(sorted_markers[0])
                
            if len (sorted_markers) > 1:
                if int(sorted_markers[1][3])/int(sorted_markers[0][3]) > float(markers_heteroyzgozity_dict[marker]):
                    sample_dict[marker].append(sorted_markers[0])
                    sample_dict[marker].append(sorted_markers[1])
                    
                else: 
                    sample_dict[marker].append(sorted_markers[0])
                    sample_dict[marker].append(sorted_markers[0])
                    
            if len (sorted_markers) == 0:
                sample_dict[marker].append(['Null', '0', 'Null', '0', 'Null'])
                sample_dict[marker].append(['Null', '0', 'Null', '0', 'Null'])
                    
            sorted_sample_dict[marker] = sorted(sample_dict[marker], key = lambda x: float(x[1]))
        
        Auto_STR_Data[sample] = sorted_sample_dict
        
        #determinate gender in HEAD from Amelogenin's alleles in Auto_STR_Data
        if Auto_STR_Data[sample]['Amelogenin'][0][1] == '0' and Auto_STR_Data[sample]['Amelogenin'][1][1] == '0':
            Auto_STR_Head[sample]['Gender'] = 'XX'
        else:
            if Auto_STR_Data[sample]['Amelogenin'][0][1] == '0' and Auto_STR_Data[sample]['Amelogenin'][1][1] == '6':
                Auto_STR_Head[sample]['Gender'] = 'XY' 
            else:
                Auto_STR_Head[sample]['Gender'] = 'N/A'
        if Auto_STR_Data[sample]['Amelogenin'][0][1] == '0':
            Auto_STR_Data[sample]['Amelogenin'][0][1] = 'X'
        if Auto_STR_Data[sample]['Amelogenin'][0][1] == '6':
            Auto_STR_Data[sample]['Amelogenin'][0][1] = 'Y'
        if Auto_STR_Data[sample]['Amelogenin'][1][1] == '0':
            Auto_STR_Data[sample]['Amelogenin'][1][1] = 'X'    
        if Auto_STR_Data[sample]['Amelogenin'][1][1] == '6':
            Auto_STR_Data[sample]['Amelogenin'][1][1] = 'Y'
            
    
    #for marker in markers:    
        #print(marker, Auto_STR_Data ['Cechova'][marker][0][1], Auto_STR_Data ['Cechova'][marker][1][1] )
        #print (sample, Auto_STR_Data[sample]['Amelogenin'][0][1], Auto_STR_Data[sample]['Amelogenin'][1][1], Auto_STR_Head[sample]['Gender'])
            
    print ('Auto_STR_Head and Auto_STR_Data done')
    
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
    
    
    #Data_CE is dictionary in formate {sample_name : {markerA : [allele1, 'allele2], markerB: [allele1, 'allele2]}}
    xml_list = os.listdir(xml_directory)
    Data_CE = {}
    for file in xml_list:
        if file.endswith('.xml'):
            path_xml = xml_directory + os.path.normpath('/') + file
            Data_CE_part = read_XML_file(path_xml)   
            #print(list(Data_CE_part.keys()))
            Data_CE = merge_two_dicts(Data_CE, Data_CE_part)
                
    
    #reading data from CSV_CE files to dictionary Data_CE for validating Auto_STR_Data
    CSV_CE_list = os.listdir(CSV_CE_directory)
    CSV_data_CE = {}
    for file in CSV_CE_list:
        if file.endswith('.csv'):
            genotype = {}
    
            for row in csv.reader(open(CSV_CE_directory + os.path.normpath('/') + file), delimiter=','):    
                specimen_id = row[9]
                locus_name = row[79].replace(' ', '')
                allele_values = row[80].split(", ")
                if not locus_name.startswith('DYS') and not locus_name.startswith('YGAT') and len(allele_values) == 1:
                    allele_values.append(allele_values[0])
                genotype[locus_name] = allele_values
    
                
            CSV_data_CE[specimen_id] = genotype
    
    
    Data_CE = merge_two_dicts(Data_CE, CSV_data_CE)
    print('Data_CE done') 
    
    
    samples_NGS = (list(Auto_STR_Data.keys()))
    samples_CE = (list(Data_CE.keys()))
    genotype_NGS = []
    
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d%H%M%S")
    report_path_name = reports_directory + os.path.normpath('/') + 'report_AutoSTR_import2db_' + dt_string + '.csv'
    report_path_name_xls = reports_directory + os.path.normpath('/') + 'report_AutoSTR_import2db_' + dt_string + '.xls'
    f = open (report_path_name, 'w+')
    f.writelines(["*** ANDY - NGS data interpreter - REPORT ***\n\n"])
    
    #validating 
    for sample_NGS in samples_NGS:
        no_mismatches = 0
        if sample_NGS in samples_CE:
            
            for marker in markers:
                if marker not in list(Data_CE[sample_NGS].keys()):
                    
                    if int(Auto_STR_Data[sample_NGS][marker][0][3]) >= no_reads_for_validation and int(Auto_STR_Data[sample_NGS][marker][1][3]) >= no_reads_for_validation:
                        Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['validated_no_reads']
                        Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['validated_no_reads']
                    else:
                        Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['not_validated_CE']
                        Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['not_validated_CE']
            
                
                if marker in list(Data_CE[sample_NGS].keys()):
                    genotype_NGS = [str(Auto_STR_Data[sample_NGS][marker][0][1])] + [str(Auto_STR_Data[sample_NGS][marker][1][1])]
                    genotype_CE = Data_CE [sample_NGS][marker]
                                
                    #sorting for comparison
                    genotype_CE.sort()
                    genotype_NGS.sort()
                    
                    if genotype_CE == genotype_NGS:
                    
                        Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['validated_CE']
                        Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['validated_CE']
                        
                    else:
                        
                        Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['locus_mismatch_CE']
                        Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['locus_mismatch_CE']
                        no_mismatches = no_mismatches + 1
            
            Auto_STR_Head[sample_NGS]['No_mismatches']= no_mismatches
        else:
            #print(sample_NGS, ' is not in xml or csv CE files')
            f.writelines([sample_NGS + "," + " is not in xml or csv CE files" + "\n"])
            Auto_STR_Head[sample_NGS]['No_mismatches']= no_mismatches
            for marker in markers:
                if int(Auto_STR_Data[sample_NGS][marker][0][3]) >= no_reads_for_validation and int(Auto_STR_Data[sample_NGS][marker][1][3]) >= no_reads_for_validation:
                    Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['validated_no_reads']
                    Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['validated_no_reads']
                else:
                    Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['not_validated_CE']
                    Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['not_validated_CE']
            
            
    print('Auto_STR_Head done')
    print('Auto_STR_Data done')
    
    #insert data from 'Auto_STR_Head' and 'Auto_STR_Data' to MySQL NGS_FORENSIC database (create dtbschema by querries in sql file)'
    #db=MySQLdb.connect("localhost", "root", dbPass, "NGS_FORENSIC")
    db=MySQLdb.connect(user="root", password=dbPass, database="NGS_FORENSIC")
    #c = db.cursor()
    c = db.cursor(buffered=True)
    sample_names = (list(Auto_STR_Data.keys()))
    for sample_name in sample_names:
        pr = Auto_STR_Head[sample_name]['Project']
        an = Auto_STR_Head[sample_name]['Analysis']
        ru = Auto_STR_Head[sample_name]['Run']
        ge = Auto_STR_Head[sample_name]['Gender']
        cr = Auto_STR_Head[sample_name]['Created']
        nm = int(Auto_STR_Head[sample_name]['No_mismatches'])
        
        # check if sample is already in database
        if CheckIfSampleInDatabase: 
            select_head_id = "SELECT id FROM Heads_flankingReg WHERE sample_name = '%s'" % (sample_name)
            c.execute(select_head_id)
            row_count = c.rowcount
        else:
            row_count = 0
        
        if row_count == 0:
            
            insert_head = "INSERT INTO Heads_flankingReg (sample_name, project, analysis, run, gender, created, no_mismatches) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%d')" % (sample_name, pr, an, ru, ge, cr, nm) 
            c.execute(insert_head)
            db.commit()
            #print(sample_name, " inserted in database")
            f.writelines(["\n" + sample_name + "," + " inserted in database"])
            select_head_id = "SELECT id FROM Heads_flankingReg WHERE sample_name = '%s' AND \
                                                    project = '%s' AND \
                                                    analysis = '%s' AND \
                                                    run = '%s' AND \
                                                    gender = '%s' AND \
                                                    created = '%s' AND \
                                                    no_mismatches = '%d' " % (sample_name, pr, an, ru, ge, cr, nm)
            c.execute(select_head_id)
            head_id = int(c.fetchone()[0])
            #print (head_id)
            for marker in markers:
                for x in range (2):
                    al = Auto_STR_Data[sample_name][marker][x][1]
                    nr = int(Auto_STR_Data[sample_name][marker][x][3])
                    seq = Auto_STR_Data[sample_name][marker][x][4]
                    val = Auto_STR_Data[sample_name][marker][x][5]
                    insert_Auto_STR = "INSERT INTO AutoSTRdata_flankingReg (sample_name, marker, allele, sequence, no_reads, CE_validation, head_id) \
                    VALUES ('%s', '%s', '%s', '%s', '%d', '%s', '%d')" % (sample_name, marker, al, seq, nr, val, head_id)
                    c.execute(insert_Auto_STR)
                    db.commit()
                    #print (sample_name, ' ', marker, ' done')
    
        else:
            #print (sample_name, 'was NOT inserted, entry already exists in database' )
            f.writelines(["\n" + sample_name + "," + " was NOT inserted - entry already exists in database"])
                
    
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
