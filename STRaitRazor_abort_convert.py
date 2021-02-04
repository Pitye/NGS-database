import os
import csv
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


input_file_path = os.path.normpath('C:/NGS_working_dir/STRaitRazor_abort.txt') 
output_file_path_AutoSTR = os.path.normpath('C:/NGS_working_dir/STRaitRazor_abort_converted_AutoSTR.txt')
output_file_path_YSTR = os.path.normpath('C:/NGS_working_dir/STRaitRazor_abort_converted_YSTR.txt')
output_file_path_XSTR = os.path.normpath('C:/NGS_working_dir/STRaitRazor_abort_NOTconverted_XSTR.txt')
output_file_path_AutoSTR_short = os.path.normpath('C:/NGS_working_dir/STRaitRazor_abort_converted_AutoSTR_short.txt')
output_file_path_YSTR_short = os.path.normpath('C:/NGS_working_dir/STRaitRazor_abort_converted_YSTR_short.txt')
 
for row in csv.reader(open (input_file_path), delimiter=','):
    
    sequence = row[3]
    
   #reverse complement of given markers 
    if row[1] in ['CSF1PO', 'D19S433', 'D1S1656', 'D2S1338', 'D5S818', 'D6S1043', 'D7S820', 'FGA', 'PentaE', 'vWA', \
                  'DYS19', 'DYS385', 'DYS389I', 'DYS389II', 'DYS390', 'DYS392', 'DYS460', 'DYS635', 'Y-GATA-H4' ]:
        my_dna = Seq(sequence)
        sequence = str(my_dna.reverse_complement())
   
    #rename
    if row[1] == 'DYS385':
        row[1]= 'DYS385a-b'
    
    
    #get shorter sequence
    if row[1] == 'PentaE':
        sequence = sequence [75:]
    if row[1] == 'TPOX':
        sequence = sequence [2:] 
    if row[1] == 'DYS392':
        sequence = sequence [104:] 
    if row[1] == 'DYS439':
        sequence = sequence [5:]
    if row[1] == 'DYS522':
        sequence = sequence [136:] 
    if row[1] == 'DYS570':
        sequence = sequence [32:]     
        
    #writing into file according to Auto,Y,X  marker   
    f= open(output_file_path_AutoSTR, "a")
    
    if str(row[1]) in ['D12S391', 'PentaD', 'FGA']:
        f= open(output_file_path_AutoSTR_short, "a")
    
    if 'DY' in str(row[1]) or 'GATA' in str(row[1]): 
        if str(row[1]) in ['DYF387S1', 'DYS385a-b', 'DYS389I', 'DYS390', 'DYS448', 'DYS460', 'DYS576', 'DYS612', 'DYS635', 'Y-GATA-H4']:
            f = open(output_file_path_YSTR_short, "a")
        else:    
            f = open(output_file_path_YSTR, "a")
    
    if 'DXS' in str(row[1]) or 'HPRTB' in str(row[1]):
        f = open(output_file_path_XSTR, "a")
        
    f.write (f'{row[0]},{"STRaitRazor_abortedRun"},{row[1]},{row[4]},{row[2]},{sequence}\n')
    f.close()

print('done')            
                




