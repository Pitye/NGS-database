import pandas as pd
import os
import shutil
import csv
#import operator

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
# D19 heterozygozyty 0.6!!!!
heterozygozity = [0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7 ]


for file in csv_list_0:
    
    writing_row_marker = 0
    head_raw = []
    row_number = 0
    
    #create empty raw_dictionary and dictionary with homo or hetero alleles for markers
    sample_raw_dict = {marker : [] for marker in markers}
    sample_dict = {marker : [] for marker in markers}
   
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
    
    #select true readings to sample_dict accord to marker's heterozygozity
    heterozygozity_place = 0            
    for marker in markers:
            
        sorted_markers = sorted(sample_raw_dict[marker], key = lambda x: int(x[3]))
        if len (sorted_markers) == 0:
            sample_dict[marker].append(['Null', 'Null', 'Null', 'Null', 'Null'])
            sample_dict[marker].append(['Null', 'Null', 'Null', 'Null', 'Null'])
            print('a')
        if len (sorted_markers) == 1:
            sample_dict[marker].append(sorted_markers[-1])
            sample_dict[marker].append(sorted_markers[-1])
            print('b')
        if len (sorted_markers) > 1:
            if int(sorted_markers[-2][3])/int(sorted_markers[-1][3]) > float(heterozygozity [heterozygozity_place]):
                sample_dict[marker].append(sorted_markers[-1])
                sample_dict[marker].append(sorted_markers[-2])
                print(int(sorted_markers[-2][3])/int(sorted_markers[-1][3]))
            else: 
                sample_dict[marker].append(sorted_markers[-1])
                sample_dict[marker].append(sorted_markers[-1])
                print(int(sorted_markers[-2][3])/int(sorted_markers[-1][3]))
        print(marker)
        print(sample_dict[marker][0][1])
        print(sample_dict[marker][1][1])
        
