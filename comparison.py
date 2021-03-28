import os
import shutil
import csv

in_GeneMarker_directory = os.path.normpath('C:/NGS_forensic_database/GeneMarker_reports') 
in_STRaitRazor_directory = os.path.normpath('C:/NGS_forensic_database/STRaitRazor_reports')
project_info = os.path.normpath('C:/NGS_forensic_database/GeneMarker_reports/project_info.txt')
GM_csv_list_0 = os.listdir(in_GeneMarker_directory)
SR_csv_list_0 = os.listdir(in_STRaitRazor_directory)
PowerSeq_markers_heteroyzgozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5,'TPOX' : 0.5,'D2S441' : 0.5,'D2S1338' : 0.5,'D3S1358' : 0.5,'FGA' : 0.5,'D5S818' : 0.5,'CSF1PO' : 0.5,'D7S820' : 0.5,'D8S1179' : 0.5,'D10S1248' : 0.5,'TH01': 0.5,'vWA' : 0.5,'D12S391' : 0.5,'D13S317' : 0.5,'PentaE' : 0.5,'D16S539' : 0.5,'D18S51' : 0.5,'D19S433' : 0.5,'D21S11' : 0.5,'PentaD' : 0.5,'D22S1045' : 0.5, 'DYS19' : 0.5, 'DYS385a/b' : 0.5}
PowerSeq_markers = list(PowerSeq_markers_heteroyzgozity_dict.keys())

GM_Auto_STR_Data = {}
GM_Auto_STR_Head = {}
GM_Auto_STR_Data_raw = {}
GM_Auto_STR_Data_unsorted = {}

SR_Auto_STR_Data = {}
SR_Auto_STR_Head = {}
SR_Auto_STR_Data_raw = {}
SR_Auto_STR_Data_unsorted = {}

PowerSeq_project_name = 'null'
PowerSeq_created_name = 'null'
PowerSeq_analysis_name = 'null'
PowerSeq_user_name = 'null'

for row in csv.reader(open(project_info), delimiter=','):
    if row[0].upper() =='PROJECT':
        PowerSeq_project_name = row [1]
    elif row[0].upper() =='CREATED':
        PowerSeq_created_name = row [1]
    elif row[0].upper() == 'ANALYSIS':
        PowerSeq_analysis_name = row[1]
    elif row[0].upper() == 'USER':
        PowerSeq_user_name = row[1]
    elif row[0].upper() == 'COMMENTS':
        continue
    else:
        print ('wrong parameter in GeneMarker project info', row[0])
        
print ('GM_project_name:', PowerSeq_project_name)
print ('GM_created_name:', PowerSeq_created_name)
print ('GM_analysis_name:', PowerSeq_analysis_name)
print ('GM_user_name:', PowerSeq_user_name)

for file in GM_csv_list_0:
    if file.endswith('.csv'):
        #GM_writing_row_marker = 0
        GM_head_raw = []
        GM_row_number = 0
    
        GM_sample_name = 'null'
        #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        GM_sample_raw_dict = {marker : [] for marker in PowerSeq_markers}
        GM_sample_dict = {marker : [] for marker in PowerSeq_markers}
        GM_sorted_sample_dict = {marker : [] for marker in PowerSeq_markers}
        GM_row_selected = []
    
        #reading csv files row by row
        for row in csv.reader(open(in_GeneMarker_directory + os.path.normpath('/') + file), delimiter=','):
        
            GM_row_number += 1
            #print (row_number, '', row [0], row [1])
            if GM_row_number == 1 and row [0] != '#Program: GeneMarkerHTS':
                print('wrong csv file or format', file )
                break
            if  GM_row_number == 3 and row[0].startswith('#Sample: '):
                if len(row[0][9:].split("_")[0]) >= 0 :
                    GM_sample_name = row[0][9:].split("_")[0]
                    
                    GM_Auto_STR_Head[GM_sample_name] = {'Project' : PowerSeq_project_name, 'Analysis': PowerSeq_analysis_name, 'Run' : 'user '+ PowerSeq_user_name, 'Gender' : 'null', 'Created' : PowerSeq_created_name}
          
        
            #reading data for samples      
            if  GM_row_number >= 6:
                
            #created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, ???Reads_ReverseSeq???, Reads_ForwardSeq, Repeat Sequence]}
                if row[0] in PowerSeq_markers:
                    if len(row)>17:
                        GM_row_selected = [row[0], row[1], row[9], row[8], row[17]]
                        #print(GM_row_selected)
                        GM_sample_raw_dict[row[0]].append(GM_row_selected)
                        GM_Auto_STR_Data_raw[GM_sample_name] = GM_sample_raw_dict
                    else:
                        print ('shortRow', GM_sample_name, row[0] )
#STRaitRazor
SR_sample_name = 'null'
for file in SR_csv_list_0:
    if file.endswith('.csv'):
        for row in csv.reader(open(in_STRaitRazor_directory + os.path.normpath('/') + file), delimiter=','):
            if row[0].split("_")[3] == 'R1':
                #new sample starts - new empty dictionary is created
                if SR_sample_name != row[0].split("_")[0]:
                    SR_Auto_STR_Head[row[0].split("_")[0]] = {'Project' : PowerSeq_project_name, 'Analysis': PowerSeq_analysis_name, 'Run' : 'user '+ PowerSeq_user_name, 'Gender' : 'null', 'Created' : PowerSeq_created_name}
                    SR_sample_raw_dict = {marker : [] for marker in PowerSeq_markers}
                #created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                SR_sample_name = row[0].split("_")[0]
                if row[1] in PowerSeq_markers:
                    SR_row_selected = [row[1], row[4], 'N/A', row[2], row[3]]
                    SR_sample_raw_dict[row[1]].append(SR_row_selected)
                    SR_Auto_STR_Data_raw[SR_sample_name] = SR_sample_raw_dict
        
#print(GM_Auto_STR_Data_raw['1-01'])
GM_samples = list(GM_Auto_STR_Data_raw.keys())
print ('GeneMarker samples: ', GM_samples)

SR_samples = list(SR_Auto_STR_Data_raw.keys())
print ('STRaitRazor samples: ', SR_samples)


def selectTrueAlleles (markers, samples, AutoSTR_Data_raw, heteroyzgozity_dict):
    AutoSTR_Data = {}
    for sample in samples:
        #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        sample_dict = {marker : [] for marker in markers}
        sorted_sample_dict = {marker : [] for marker in markers}
    
        #select true readings to sample_dict according to marker's heterozygozity 
        for marker in markers:
            sorted_markers = sorted(AutoSTR_Data_raw[sample][marker], key = lambda x: int(x[3]),reverse = True)
        
            if len (sorted_markers) == 1:
                sample_dict[marker].append(sorted_markers[0])
                sample_dict[marker].append(sorted_markers[0])
            
            if len (sorted_markers) > 1:
                if int(sorted_markers[1][3])/int(sorted_markers[0][3]) > float(heteroyzgozity_dict[marker]):
                    sample_dict[marker].append(sorted_markers[0])
                    sample_dict[marker].append(sorted_markers[1])
                
                else: 
                    sample_dict[marker].append(sorted_markers[0])
                    sample_dict[marker].append(sorted_markers[0])
                
            if len (sorted_markers) == 0:
                sample_dict[marker].append(['Null', '0', 'Null', '0', 'Null'])
                sample_dict[marker].append(['Null', '0', 'Null', '0', 'Null'])
                 
            if sample_dict['Amelogenin'][0][1] == 'chrX':
                sample_dict['Amelogenin'][0][1] = '0'         
            if sample_dict['Amelogenin'][0][1] == 'chrY' or sample_dict['Amelogenin'][0][1] == '1' :
                sample_dict['Amelogenin'][0][1] = '6'    
        
            if sample_dict['Amelogenin'][1][1] == 'chrX':
                sample_dict['Amelogenin'][1][1] = '0'
            if sample_dict['Amelogenin'][1][1] == 'chrY' or sample_dict['Amelogenin'][1][1] == '1': 
                sample_dict['Amelogenin'][1][1] = '6'
            
            if sample_dict[marker][0][1] == 'Unknown':
                sample_dict[marker][0][1] = '999'
            if sample_dict[marker][1][1] == 'Unknown':
                sample_dict[marker][1][1] = '999'
                
            sorted_sample_dict[marker] = sorted(sample_dict[marker], key = lambda x: float(x[1]))
    
        AutoSTR_Data[sample] = sorted_sample_dict
    return AutoSTR_Data

GM_Auto_STR_Data = selectTrueAlleles (PowerSeq_markers, GM_samples, GM_Auto_STR_Data_raw, PowerSeq_markers_heteroyzgozity_dict)
SR_Auto_STR_Data = selectTrueAlleles (PowerSeq_markers, SR_samples, SR_Auto_STR_Data_raw, PowerSeq_markers_heteroyzgozity_dict)

def markerAlleleList(AutoSTR_Data, sample, marker):
    AlleleList = [AutoSTR_Data[sample][marker][0][1],AutoSTR_Data[sample][marker][1][1]]
    return AlleleList
def markerNoReadsList(AutoSTR_Data, sample, marker):
    NoReadsList = [AutoSTR_Data[sample][marker][0][3],AutoSTR_Data[sample][marker][1][3]]
    return NoReadsList

for sample in GM_samples:
    if sample in SR_samples:
        for marker in PowerSeq_markers:
            GM_AlleleList = markerAlleleList(GM_Auto_STR_Data, sample, marker)
            SR_AlleleList = markerAlleleList(SR_Auto_STR_Data, sample, marker)
            GM_NoReadsList = markerNoReadsList(GM_Auto_STR_Data, sample, marker)
            SR_NoReadsList = markerNoReadsList(SR_Auto_STR_Data, sample, marker)
            if GM_AlleleList != SR_AlleleList:
                print (sample, marker, 'GeneMarker: ', GM_AlleleList,GM_NoReadsList, 'STRaitRazor: ', SR_AlleleList, SR_NoReadsList)
    else:
        print (sample, 'is not in STRaitRazor samples')

 