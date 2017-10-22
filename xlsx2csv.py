# NGS-database
# xlsx2csv
import numpy as np
import pandas as pd
import os
in_directory = 'Desktop/NGS/xlsx/'
out_directory = 'Desktop/NGS/csv/'
sheets = ['Autosomal STRs', 'Y STRs', 'X STRs', 'iSNPs']
xlsx_list = os.listdir(in_directory)
for file in xlsx_list:
    for sheet in sheets:
        path_xlsx = in_directory + file
        path_csv = out_directory + file[0:-5] + sheet[0] + '.csv'
        df = pd.read_excel(path_xlsx, sheet, index_col=None, na_values=['NA'])
        df.to_csv(path_csv)
print ('hotovo') 
