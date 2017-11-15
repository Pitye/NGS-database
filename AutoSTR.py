import pandas as pd
import os
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter

in_directory = '/Users/sweethome/Desktop/project_NGS/NGS/xlsx/'
out_directory = '/Users/sweethome/Desktop/project_NGS/NGS/csv/'
xml_directory = '/Users/sweethome/Desktop/project_NGS/NGS/xml_CE/'
sheets = ['Autosomal STRs', 'Y STRs', 'X STRs', 'iSNPs']
no_reads_for_validation = 400

#delete all in OUT DIRECTORY
shutil.rmtree(out_directory)
os.makedirs(out_directory)

#create OUT DIRECTORIES_A,Y,X,i in CSV directory
for sheet in sheets:
    os.makedirs(out_directory + sheet)

#create CSV files to OUT DIRECTORIES_A,Y,X
xlsx_list = os.listdir(in_directory)
for file in xlsx_list:
    if file.endswith('.xlsx'):
        for sheet in sheets:
            path_xlsx = in_directory + file
            path_csv = out_directory + file[0:-5] + '_' + sheet[0] + '.csv'
            path_csv = out_directory + sheet + '/' + file[0:-5] + '_' + sheet[0] + '.csv'
            df = pd.read_excel(path_xlsx, sheet, header=None)
            df.to_csv(path_csv, header=None, index=None)
print('xlsx2csv done') 

#reading CSV files in directory for sheet[0] - 'Autosomal STRs'
csv_list_0 = os.listdir(out_directory + sheets[0] + '/')
#markers_heteroyzgozity_dict = {'D1S1656': 0.7,'TPOX' : 0.7,'D2S441' : 0.7,'D2S1338' : 0.7,'D3S1358' : 0.7,'D4S2408' : 0.7,'FGA' : 0.7,'D5S818' : 0.7,'CSF1PO' : 0.7,'D6S1043' : 0.7,'D7S820' : 0.7,'D8S1179' : 0.7,'D9S1122' : 0.7,'D10S1248' : 0.7,'TH01': 0.7,'vWA' : 0.7,'D12S391' : 0.7,'D13S317' : 0.7,'PentaE' : 0.7,'D16S539' : 0.7,'D17S1301' : 0.7,'D18S51' : 0.7,'D19S433' : 0.7,'D20S482' : 0.7,'D21S11' : 0.7,'PentaD' : 0.7,'D22S1045' : 0.7}
markers_heteroyzgozity_dict = {'D1S1656': 0.5,'TPOX' : 0.5,'D2S441' : 0.5,'D2S1338' : 0.5,'D3S1358' : 0.5,'D4S2408' : 0.5,'FGA' : 0.5,'D5S818' : 0.5,'CSF1PO' : 0.5,'D6S1043' : 0.5,'D7S820' : 0.5,'D8S1179' : 0.5,'D9S1122' : 0.5,'D10S1248' : 0.5,'TH01': 0.5,'vWA' : 0.5,'D12S391' : 0.5,'D13S317' : 0.5,'PentaE' : 0.5,'D16S539' : 0.5,'D17S1301' : 0.5,'D18S51' : 0.5,'D19S433' : 0.5,'D20S482' : 0.5,'D21S11' : 0.5,'PentaD' : 0.5,'D22S1045' : 0.5}
markers = list(markers_heteroyzgozity_dict.keys())

Auto_STR_Data = {}
Auto_STR_Head = {}

for file in csv_list_0:
    
    writing_row_marker = 0
    head_raw = []
    row_number = 0
    
    #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
    sample_raw_dict = {marker : [] for marker in markers}
    sample_dict = {marker : [] for marker in markers}
    sorted_sample_dict = {marker : [] for marker in markers}
    
    #reading csv files row by row
    for row in csv.reader(open(out_directory + sheets[0] + '/'+ file), delimiter=','):
        
        row_number += 1
        
        if row [0] == 'Locus':
            writing_row_marker += 1
        
        if writing_row_marker == 0:
            head_raw.append(row)
        
        #writing to raw_dictionary for markers
        if writing_row_marker == 2:
            if row[0] in markers:
                sample_raw_dict[row[0]].append(row)
    
    #select true readings to sample_dict according to marker's heterozygozity              
    for marker in markers:
            
        sorted_markers = sorted(sample_raw_dict[marker], key = lambda x: int(x[3]),reverse = True)
        
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
        
        
    head_dict = {'Project' : head_raw[3][1], 'Analysis': head_raw[4][1], 'Run' : head_raw[5][1], 'Gender' : head_raw[6][1], 'Created' : head_raw[7][1]}

    sample_name = head_raw[2][1]
    
    # Auto_STR_Data is dictionary in formate {sample_name : {marker : [[marker, allele1, Yes/No, no.readings, seq.], [marker, allele2, Yes/No, no.readings, seq.]]} }
    Auto_STR_Data[sample_name] = sorted_sample_dict
    Auto_STR_Head[sample_name] = head_dict
    
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
        path_xml = xml_directory + file
        Data_CE_part = read_XML_file(path_xml)   
        #print(list(Data_CE_part.keys()))
        Data_CE = merge_two_dicts(Data_CE, Data_CE_part)
            
print('Data_CE done') 


samples_NGS = (list(Auto_STR_Data.keys()))
samples_CE = (list(Data_CE.keys()))
genotype_NGS = []

#validating 
for sample_NGS in samples_NGS:
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
                               
                if genotype_CE == genotype_NGS:
                   
                    Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['validated_CE']
                    Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['validated_CE']
                    
                else:
                    
                    Auto_STR_Data[sample_NGS][marker][0] = Auto_STR_Data[sample_NGS][marker][0] + ['locus_mismatch_CE']
                    Auto_STR_Data[sample_NGS][marker][1] = Auto_STR_Data[sample_NGS][marker][1] + ['locus_mismatch_CE']
    else:
        print(sample_NGS, ' is not in xml files')

print('Auto_STR_Head done')
print('Auto_STR_Data done')
