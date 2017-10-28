import numpy as np
import pandas as pd
import os
import shutil
import csv

in_directory = 'Desktop/NGS/xlsx/'
out_directory = 'Desktop/NGS/csv/'
sheets = ['Autosomal STRs', 'Y STRs', 'X STRs', 'iSNPs']

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
markers = ['D1S1656','TPOX','D2S441','D2S1338','D3S1358','D4S2408','FGA','D5S818','CSF1PO','D6S1043','D7S820','D8S1179','D9S1122','D10S1248','TH01','vWA','D12S391','D13S317','PentaE','D16S539','D17S1301','D18S51','D19S433','D20S482','D21S11','PentaD','D22S1045']

for file in csv_list_0:
    
    writing_row_marker = 0
    head_raw = []
    row_number = 0
    
    #create empty dictionary for markers
    sample_raw_dict = {marker : [] for marker in markers}
    
    #reading CSV file row by row
    for row in csv.reader(open(out_directory + sheets[0] + '/'+ file), delimiter=','):
        
        row_number += 1
        
        if row [0] == 'Locus':
            writing_row_marker += 1
        
        if writing_row_marker == 0:
            head_raw.append(row)
        
        if writing_row_marker == 2:
            
            #writing to dictionary for markers
            if row[0] in markers:
                sample_raw_dict[row[0]].append(row)   
            
             
    print (head_raw)
    print (sample_raw_dict) 
