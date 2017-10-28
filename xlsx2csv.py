import pandas as pd
import os
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
