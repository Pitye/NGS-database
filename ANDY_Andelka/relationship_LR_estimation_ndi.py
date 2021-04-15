import pandas as pd
import os
import sys
import shutil
import csv
import xml.etree.ElementTree as etree
import operator
from collections import Counter
#import MySQLdb
import mysql.connector as MySQLdb
import math
from datetime import datetime
import xlwt

#import mysql.connector
#from mysql.connector import Error

def main():

    ### *** parameter settings ***
    sample_names=['AA11-KIGEN-SF26', 'MK2_S2_L001_R1_001']
    length_polymorphism_estimation = True
    use_external_file_for_missing_markers = True
    directory_STR_lenght_profiles = os.path.normpath('C:/NGS_forensic_database/relationship_LR_estimation/STR_lenght_profiles')
    directory_STR_lenght_frequencies = os.path.normpath('C:/NGS_forensic_database/relationship_LR_estimation/STR_lenght_frequencies')
    directory_reports = os.path.normpath('C:/NGS_forensic_database/relationship_LR_estimation/reports')
    dbPass = 'XXX'
    
    relationship_LR_estimation (sample_names, length_polymorphism_estimation, use_external_file_for_missing_markers, directory_STR_lenght_profiles, directory_STR_lenght_frequencies, directory_reports, dbPass)
    
    
def relationship_LR_estimation (sample_names, length_polymorphism_estimation, use_external_file_for_missing_markers, directory_STR_lenght_profiles, directory_STR_lenght_frequencies, directory_reports, dbPass):
    print ('Please wait ...')
    if length_polymorphism_estimation == True:
        seq_lenght_types = ['lenght', 'seq']
    else:
        seq_lenght_types = ['seq']
    
    markers_dict = {'D1S1656': 0.5,'TPOX' : 0.5,'D2S441' : 0.5,'D2S1338' : 0.5,'D3S1358' : 0.5,'D4S2408' : 0.5,'FGA' : 0.5,'D5S818' : 0.5,'CSF1PO' : 0.5,'D6S1043' : 0.5,'D7S820' : 0.5,'D8S1179' : 0.5,'D9S1122' : 0.5,'D10S1248' : 0.5,'TH01': 0.5,'vWA' : 0.5,'D12S391' : 0.5,'D13S317' : 0.5,'PentaE' : 0.5,'D16S539' : 0.5,'D17S1301' : 0.5,'D18S51' : 0.5,'D19S433' : 0.5,'D20S482' : 0.5,'D21S11' : 0.5,'PentaD' : 0.5,'D22S1045' : 0.5}
    markers = list(markers_dict.keys())
    markers_evaluated = list(markers_dict.keys())
    columns = ['allele', 'seq_name', 'PubMed_ID', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']
    
    
    ### *** create empty dict for data from mySQL dtb and external files ***
    
    sample_records_original = {sample: {marker: {column: [] for column in columns} for marker in markers} for sample in sample_names}
    sample_records = {sample: {marker: {column: [] for column in columns} for marker in markers} for sample in sample_names}
    STR_profiles_records = {sample: {marker : {'STRallele' :[]} for marker in markers} for sample in sample_names}
    STR_frequencies = {marker:{} for marker in markers}
    STR_frequencies_database = {marker:{} for marker in markers}
    missed_markers = {sample: [] for sample in sample_names}
    result = {marker: {seq_lenght_type : {'alleles' : [], 'p/q_freq' : [], 'formulas': [], 'calculations' : []} for seq_lenght_type in seq_lenght_types } for marker in markers}
    result_PI = {seq_lenght_type : 1 for seq_lenght_type in seq_lenght_types}
    result_FSI = {seq_lenght_type : 1 for seq_lenght_type in seq_lenght_types}
    result_GI_AI_HIS = {seq_lenght_type : 1 for seq_lenght_type in seq_lenght_types}
    result_FirsCI = {seq_lenght_type : 1 for seq_lenght_type in seq_lenght_types}
    
    ### *** read external files - STR_Lenght_profiles ***
    if use_external_file_for_missing_markers == True:
        file_list_STR_profiles = os.listdir(directory_STR_lenght_profiles)
        for file in file_list_STR_profiles:
            line_number = 0
            for line in csv.reader(open(directory_STR_lenght_profiles + os.path.normpath('/') + file), delimiter=';'):
                if line_number == 0:
                    STRsample = line[0]
                    if STRsample not in sample_names:
                        print ('external STR profile ' + STRsample + ' will not be estimated, correct sample names: ')
                        print ( sample_names,"\n" )
                        break
                else:
                    if len(line) == 3 and line[0] in markers:
                        STR_profiles_records[STRsample][line[0]]['STRallele'].append(line[1])
                        STR_profiles_records[STRsample][line[0]]['STRallele'].append(line[2])
                    else:
                        print ('not correct structure in external STR profile ' + STRsample + ' line: ' + str(line_number+1) )
                line_number += 1    
    
    ### *** read external files - STR_Lenght_frequencies ***
        file_list_STR_frequencies = os.listdir(directory_STR_lenght_frequencies)
        for file in file_list_STR_frequencies:
            STRmarker = ''
            for line in csv.reader(open(directory_STR_lenght_frequencies + os.path.normpath('/') + file), delimiter=';'):
                if len(line) == 1:
                    STRmarker = line[0]
                    if STRmarker not in markers:
                        print ('external STR frequency for ' + STRmarker + ' will not be estimated, it is not in markers for databasing')
                if len(line) == 2 and STRmarker in markers:
                    dict1 = {line[0]:line[1]}
                    STR_frequencies[STRmarker].update(dict1)
    
    
    
    
    
    
    
    ### *** get data from from mySQL dtb - sample_records ***
    #db=MySQLdb.connect("localhost", "root", dbPass, "NGS_FORENSIC")
    db=MySQLdb.connect(user="root", password=dbPass, database="NGS_FORENSIC")
    #c = db.cursor()
    c = db.cursor(buffered=True)
    
    for sample in sample_names:
        
        sql_select_Query = "SELECT * FROM NGS_FORENSIC.nomen_freq_autostrdata_flankingreg where sample_name = '%s'" % (sample)
        c.execute(sql_select_Query)
        records = c.fetchall()
    
        #print("Total number of rows for sample1: ", c.rowcount)
    
        #print("\nPrinting sample " + sample )
        for row in records:
            #marker, row[1], allele, row[2], seq_name, row[3], PubMed_ID, row[4], sequence, row[5], no_reads, row[6], CE_validation, row[7], head_id, row[8], avg_no_reads, row[9], count_seq, row[10], frequency, row[11]
    
            column_index = 2
            for column in columns:
                sample_records[sample][row[1]][column].append(row[column_index])
                column_index += 1
    
        
        for marker in markers:
            for column in columns:
                sample_records_original[sample][marker][column] = sample_records[sample][marker][column]
                if sample_records_original[sample][marker][column] == []:
                    sample_records_original[sample][marker][column] = ['null','null']
                
    
    ### *** get data from from mySQL dtb - STR_frequencies_database ***
    
    sql_select_Query2 = "SELECT marker, allele, sum(frequency) FROM NGS_FORENSIC.nomenclature_markerautostrview_flankingreg group by marker, allele order by marker, allele"           
    c.execute(sql_select_Query2)
    records2 = c.fetchall()
    for row in records2:
        marker = row[0]
        dict_freq_database = {row[1]:row[2]}
        STR_frequencies_database[marker].update(dict_freq_database)
    
    db.close() 
    
    ### *** check which markers don't have any data and use external STR lenght profiles and STR lenght frecuencies ***
    
    for sample in sample_names:
        for marker in markers:
    
            if sample_records[sample][marker]['allele'] == []:
                missed_markers[sample].append(marker)   
                if use_external_file_for_missing_markers == True:
                    if STR_profiles_records[sample][marker]['STRallele']!=[]:
                        sample_records[sample][marker]['allele'].append(STR_profiles_records[sample][marker]['STRallele'][0])
                        sample_records[sample][marker]['allele'].append(STR_profiles_records[sample][marker]['STRallele'][1])
                        missed_markers[sample].remove(marker)
                    else:
                        print('for sample ' + sample + ' and marker ' + marker + ' sequential data are missing and no STR data provided')
    
    
        if missed_markers[sample]!=[]:
            for missed_marker in missed_markers[sample]:
                if missed_marker in markers_evaluated:
                    markers_evaluated.remove(missed_marker)
                    print ('missed markers: '+ missed_marker)
    
    
    ### *** find corresponding formula 
    def setColumn (param_column):
        if param_column == 'lenght':
            Fcolumn = 'allele'
        elif param_column == 'seq':
            Fcolumn = 'PubMed_ID'
        return Fcolumn
    
    def setAlleleABCD (letter, param_locus, param_seq_lenght_type):
        column = setColumn(param_seq_lenght_type)
        
        if letter == 'A':
            AlleleABCD = sample_records[sample_names[0]][param_locus][column][0]
        elif letter == 'B':
            AlleleABCD = sample_records[sample_names[0]][param_locus][column][1] 
        elif letter == 'C':
            AlleleABCD = sample_records[sample_names[1]][param_locus][column][0]
        elif letter == 'D':
            AlleleABCD = sample_records[sample_names[1]][param_locus][column][1]
        return AlleleABCD
        
    def isHomozygote(H_sample, H_locus, H_seq_lenght_type):
        column = setColumn(H_seq_lenght_type) 
        
        if sample_records[H_sample][H_locus][column][0] == sample_records[H_sample][H_locus][column][1]:
            return True
        else:
            return False
        
    def isAlleleShared(S_locus, S_allele, S_seq_lenght_type):
        column = setColumn(S_seq_lenght_type)
        checked_allele = setAlleleABCD(S_allele, S_locus, S_seq_lenght_type )
        if S_allele == 'A' or S_allele == 'B':
            checked_allele_list = sample_records[sample_names[1]][S_locus][column]
        elif S_allele == 'C' or S_allele == 'D':
            checked_allele_list = sample_records[sample_names[0]][S_locus][column] 
        else:
            print ('isAlleleShared wrong argument value allele')
            return
    
    
        if checked_allele in checked_allele_list:
            return True
        else:
            return False
    
    def setFrequency (p_q, locus, seq_lenght_type):
        if seq_lenght_type == 'seq' and sample_records[sample_names[0]][locus]['frequency'] != []:
            if p_q == 'p':
                frequency = sample_records[sample_names[0]][locus]['frequency'][0]
            if p_q == 'q':
                frequency = sample_records[sample_names[0]][locus]['frequency'][1]
        else:
    
            if p_q == 'p':
                allele = setAlleleABCD ('A', locus, 'lenght')
            if p_q == 'q':
                allele = setAlleleABCD ('B', locus, 'lenght')
    
            if allele in list(STR_frequencies[locus].keys()):   
                frequency = STR_frequencies[locus][allele]
            else:
                frequency = STR_frequencies_database[locus][allele]
    
        return frequency
    
    def round_up(n, decimals=0):
        multiplier = 10 ** decimals
        return math.ceil(n * multiplier) / multiplier
    
    def setFormulaCalculation(formula_type, locus, seq_lenght_type):
        p = float(setFrequency('p', locus, seq_lenght_type))
        q = float(setFrequency('q', locus, seq_lenght_type))
        if formula_type == '1':
            formula = ['1/p', '(p+1)(p+1)/(4*p*p)', '(p+1)/(2*p)', '(3*p+1)/(4*p)']
            result = [round_up(1/(p),15), round_up((p+1)*(p+1)/(4*p*p),15), round_up((p+1)/(2*p),15), round_up((3*p+1)/(4*p),15)]
    
        if formula_type == '2':
            formula = ['1/(2*p)','(p+1)/(4*p)', '(2*p+1)/(4*p)', '(6*p+1)/(8*p)']
            result = [round_up(1/(2*p),15), round_up((p+1)/(4*p),15), round_up((2*p+1)/(4*p),15), round_up((6*p+1)/(8*p),15)]
        if formula_type == '3':
            formula = ['1/(2*q)','(q+1)/(4*q)', '(2*q+1)/(4*q)', '(6*q+1)/(8*q)']
            result = [round_up(1/(2*q),15), round_up((q+1)/(4*q),15), round_up((2*q+1)/(4*q),15), round_up((6*q+1)/(8*q),15)]
        if formula_type == '4':
            formula = ['1/(4*p)', '(2*p+1)/(8*p)', '(4*p+1)/(8*p)', '(12*p+1)/(16*p)']
            result = [round_up(1/(4*p),15), round_up((2*p+1)/(8*p),15), round_up((4*p+1)/(8*p),15), round_up((12*p+1)/(16*p),15)]
        if formula_type == '5':
            formula = ['1/(4*q)', '(2*q+1)/(8*q)', '(4*q+1)/(8*q)', '(12*q+1)/(16*q)']
            result = [round_up(1/(4*q),15), round_up((2*q+1)/(8*q),15), round_up((4*q+1)/(8*q),15), round_up((12*q+1)/(16*q),15)]
        if formula_type == '6':
            formula = ['(p+q)/(4*p*q)', '(2*p*q+p+q+1)/(8*p*q)', '(4*p*q+p+q)/(8*p*q)', '(12*p*q+p+q+1)/(16*p*q)' ]
            result = [round_up((p+q)/(4*p*q),15), round_up((2*p*q+p+q+1)/(8*p*q),15), round_up((4*p*q+p+q)/(8*p*q),15), round_up((12*p*q+p+q+1)/(16*p*q),15)]
        if formula_type == '7':
            formula = ['*', '0.25', '0.5', '0.75']
            result = [0, 0.25, 0.5, 0.75]
        FormulaCalculation = {'formula' : formula, 'calculation': result}
    
        return FormulaCalculation
    
    ### *** find corresponding formula
    for seq_lenght_type in seq_lenght_types:
        column = setColumn(seq_lenght_type)
        for marker in markers_evaluated:
            if sample_records[sample_names[0]][marker][column] != [] and sample_records[sample_names[1]][marker][column] != []:
                if (isHomozygote(sample_names[0], marker, seq_lenght_type) == True) and (isHomozygote(sample_names[1], marker, seq_lenght_type) == True):
    
                    if setAlleleABCD ('A', marker, seq_lenght_type) == setAlleleABCD ('C', marker, seq_lenght_type):
                        formula_type = '1'
                    else:
                        formula_type = '7'
    
                elif (isHomozygote(sample_names[0], marker, seq_lenght_type) == True) and (isHomozygote(sample_names[1], marker, seq_lenght_type) == False):
    
                    if isAlleleShared(marker, "A", seq_lenght_type) == True:
                        formula_type = '2'
                    else:
                        formula_type = '7'
    
                elif isHomozygote(sample_names[0], marker, seq_lenght_type) == False and isHomozygote(sample_names[1], marker, seq_lenght_type) == True:
    
                    if setAlleleABCD ('A', marker, seq_lenght_type) == setAlleleABCD ('C', marker, seq_lenght_type):  
                        formula_type = '2'
                    elif setAlleleABCD ('B', marker, seq_lenght_type) == setAlleleABCD ('C', marker, seq_lenght_type):                                                       
                        formula_type = '3'
                    else:
                        formula_type = '7'
    
                elif isHomozygote(sample_names[0], marker, seq_lenght_type) == False and isHomozygote(sample_names[1], marker, seq_lenght_type) == False:
    
                    if isAlleleShared(marker, "A", seq_lenght_type) == True and isAlleleShared(marker, "B", seq_lenght_type) == True:
                        formula_type = '6' 
                    elif isAlleleShared(marker, "A", seq_lenght_type) == True and isAlleleShared(marker, "B", seq_lenght_type) == False:
                        formula_type = '4'
                    elif isAlleleShared(marker, "A", seq_lenght_type) == False and isAlleleShared(marker, "B", seq_lenght_type) == True:
                        formula_type = '5'
                    
                    else:
                        formula_type = '7' 
    
    
    
                ### *** result = {marker: {seq_lenght_type : {'alleles' : [], 'formulas': [], 'calculations' : []} for seq_lenght_type in seq_lenght_types } for marker in markers}
    
                result[marker][seq_lenght_type]['alleles'] = [setAlleleABCD ('A', marker, seq_lenght_type), setAlleleABCD ('B', marker, seq_lenght_type), setAlleleABCD ('C', marker, seq_lenght_type), setAlleleABCD ('D', marker, seq_lenght_type)]
                result[marker][seq_lenght_type]['p/q_freq'] = [setFrequency ('p', marker, seq_lenght_type), setFrequency ('q', marker, seq_lenght_type)]
                result[marker][seq_lenght_type]['formulas'] = setFormulaCalculation(formula_type, marker, seq_lenght_type)['formula']
                result[marker][seq_lenght_type]['calculations'] = setFormulaCalculation(formula_type, marker, seq_lenght_type)['calculation']
    
    ### *** if seq allele is missing, add lenght allele
    for marker in markers_evaluated:
        if result[marker]['seq']['alleles'] == [] and result[marker]['lenght']['alleles'] != []:
            
            result[marker]['seq']['alleles'] = result[marker]['lenght']['alleles']
            result[marker]['seq']['p/q_freq'] = result[marker]['lenght']['p/q_freq'] 
            result[marker]['seq']['formulas'] = result[marker]['lenght']['formulas']
            result[marker]['seq']['calculations'] = result[marker]['lenght']['calculations']
            
    for seq_lenght_type in seq_lenght_types:
        column = setColumn(seq_lenght_type)
        for marker in markers_evaluated:        
            
            if result[marker][seq_lenght_type]['calculations'] != []:
                result_PI[seq_lenght_type] =round_up( float(result_PI[seq_lenght_type]) * float(result[marker][seq_lenght_type]['calculations'][0]),20)
                result_FSI[seq_lenght_type] =round_up( float( result_FSI[seq_lenght_type]) * float(result[marker][seq_lenght_type]['calculations'][1]),20)
                result_GI_AI_HIS[seq_lenght_type] =round_up( float(result_GI_AI_HIS[seq_lenght_type]) * float(result[marker][seq_lenght_type]['calculations'][2]),20)
                result_FirsCI[seq_lenght_type] =round_up( float(result_FirsCI[seq_lenght_type]) * float(result[marker][seq_lenght_type]['calculations'][3]),20)
    
    for marker in markers_evaluated:
        for seq_lenght_type in seq_lenght_types:
            if result[marker][seq_lenght_type]['alleles'] == []:
                result[marker][seq_lenght_type]['alleles'] = ['null', 'null', 'null', 'null']
            if result[marker][seq_lenght_type]['p/q_freq'] == []:
                result[marker][seq_lenght_type]['p/q_freq'] = ['null', 'null']
            if result[marker][seq_lenght_type]['formulas'] == []:
                result[marker][seq_lenght_type]['formulas'] = ['null', 'null', 'null', 'null']
            if result[marker][seq_lenght_type]['calculations'] == []:
                result[marker][seq_lenght_type]['calculations'] = ['null', 'null', 'null', 'null']
    
    
    ### *** create report ***
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d%H%M%S")
    report_path_name = directory_reports + os.path.normpath('/') + 'report_' + sample_names[0] + '_' + sample_names[1] + '_' + dt_string + '.csv'
    report_path_name_xls = directory_reports + os.path.normpath('/') + 'report_' + sample_names[0] + '_' + sample_names[1] + '_' + dt_string + '.xls'
    f = open (report_path_name, 'w+')
    f.writelines([ "*** ANDY - NGS data interpreter - RELATIONSHIP_LR_ESTIMATION_REPORT ***", \
            "\n\n", \
            "DATA SELECTED FROM DATABASE", \
            "\n", \
            "\n" + sample_names[0] + ", allele_1, allele_2, PubMed_ID_1, PubMed_ID_2, no_reads_1, no_reads_2, CE_validation_1, CE_validation_2" + ",," + \
                sample_names[1] + ", allele_1, allele_2, PubMed_ID_1, PubMed_ID_2, no_reads_1, no_reads_2, CE_validation_1, CE_validation_2"])

    for marker in markers:
        f.writelines(["\n" + marker + "," + sample_records_original[sample_names[0]][marker]['allele'][0] + "," + sample_records_original[sample_names[0]][marker]['allele'][1] + "," + \
                sample_records_original[sample_names[0]][marker]['PubMed_ID'][0] + "," + sample_records_original[sample_names[0]][marker]['PubMed_ID'][1] + "," + \
                str(sample_records_original[sample_names[0]][marker]['no_reads'][0]) + "," + str(sample_records_original[sample_names[0]][marker]['no_reads'][1]) + "," + \
                sample_records_original[sample_names[0]][marker]['CE_validation'][0] + "," + sample_records_original[sample_names[0]][marker]['CE_validation'][1] +\
                ",,," + \
                sample_records_original[sample_names[1]][marker]['allele'][0] + "," + sample_records_original[sample_names[1]][marker]['allele'][1] + "," + \
                sample_records_original[sample_names[1]][marker]['PubMed_ID'][0] + "," + sample_records_original[sample_names[1]][marker]['PubMed_ID'][1] + "," + \
                str(sample_records_original[sample_names[1]][marker]['no_reads'][0]) + "," + str(sample_records_original[sample_names[1]][marker]['no_reads'][1]) + "," + \
                sample_records_original[sample_names[1]][marker]['CE_validation'][0] + "," + sample_records_original[sample_names[1]][marker]['CE_validation'][1] ])
    
    f.writelines([ "\n\n\nLR_ESTIMATION\n", \
                "\nlength_polymorphism_estimation: " + str(length_polymorphism_estimation), \
                "\nuse_external_file_for_missing_markers: " + str(use_external_file_for_missing_markers)])
    missing_markers = missed_markers[sample_names[0]] + missed_markers[sample_names[1]]
    missing_markers_set = set(missing_markers)
    for seq_lenght_type in seq_lenght_types:
        f.writelines([ "\n\n" + seq_lenght_type + "_ESTIMATION\n", \
                    "\nMarker, allele1, allele2, allele3, allele4, freq(p), freq(q), PI, FSI, GI/AI/HIS, First CI, LR-PI, LR-FSI, LR-GI/AI/HIS, LR-First CI" ])
        for marker in markers_evaluated:
            f.writelines(["\n" + marker + "," + result[marker][seq_lenght_type]['alleles'][0] + "," + result[marker][seq_lenght_type]['alleles'][1] + "," + \
                                                result[marker][seq_lenght_type]['alleles'][2] + "," + result[marker][seq_lenght_type]['alleles'][3] + "," + \
                                                str(result[marker][seq_lenght_type]['p/q_freq'][0]) + "," + str(result[marker][seq_lenght_type]['p/q_freq'][1]) + "," + \
                                                result[marker][seq_lenght_type]['formulas'][0] + "," + result[marker][seq_lenght_type]['formulas'][1] + "," + \
                                                result[marker][seq_lenght_type]['formulas'][2] + "," + result[marker][seq_lenght_type]['formulas'][3] + "," + \
                                                str(result[marker][seq_lenght_type]['calculations'][0]) + "," + str(result[marker][seq_lenght_type]['calculations'][1]) + "," + \
                                                str(result[marker][seq_lenght_type]['calculations'][2]) + "," + str(result[marker][seq_lenght_type]['calculations'][3]) ])
    
        f.writelines(["\nRESULTS,,,,,,,,,,," + str(result_PI[seq_lenght_type]) + "," + str(result_FSI[seq_lenght_type]) + "," +  str(result_GI_AI_HIS[seq_lenght_type]) + "," + str(result_FirsCI[seq_lenght_type]), \
                    "\n\nmissed markers:," + str(missing_markers_set) ]) 
    f.close()
    wb = xlwt.Workbook()
    sh = wb.add_sheet('report')
    with open(report_path_name, 'r') as fc:
        reader = csv.reader(fc)
        for r, row in enumerate(reader):
            for c, val in enumerate(row):
                sh.write(r, c, val)
    wb.save(report_path_name_xls)
    
    print ('markers evaluated for LR: ', markers_evaluated)
    print ('report exported')
    print('all done')

if __name__ == '__main__':
    main()

