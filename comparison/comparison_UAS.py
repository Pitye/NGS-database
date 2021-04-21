import os
import shutil
import csv
import pandas as pd
import sys


def systemLinux():
    if sys.platform == 'linux':
        return True


def getHome():
    if systemLinux():
        print('system is Linux')
        homePath = os.path.normpath('/home/pavla/')
    else:
        print('system is Windows')
        homePath = os.path.normpath('C:/')
    return homePath


def createPowerSeqProjectInfo(inProjectInfo):
    # project info from GeneMarker
    PowerSeq_project_name = 'null'
    PowerSeq_created_name = 'null'
    PowerSeq_analysis_name = 'null'
    PowerSeq_user_name = 'null'
    for row in csv.reader(open(inProjectInfo), delimiter=','):
        if row[0].upper() == 'PROJECT':
            PowerSeq_project_name = row[1]
        elif row[0].upper() == 'CREATED':
            PowerSeq_created_name = row[1]
        elif row[0].upper() == 'ANALYSIS':
            PowerSeq_analysis_name = row[1]
        elif row[0].upper() == 'USER':
            PowerSeq_user_name = row[1]
        elif row[0].upper() == 'COMMENTS':
            continue
        else:
            print('wrong parameter in PowerSeq project info', row[0])
    PowerSeq_project_info = [PowerSeq_project_name, PowerSeq_created_name, PowerSeq_analysis_name, PowerSeq_user_name]
    return PowerSeq_project_info


def createUasRaw(inUasDirectory, outUasDirectory, IlluminaMarkers):
    UAS_Auto_STR_Head = {}
    UAS_AutoSTRData_raw = {}
    UAS_sheets = ['Autosomal STR Coverage', 'X STR Coverage', 'Y STR Coverage', 'iSNP Coverage']

    # delete all in OUT DIRECTORY
    shutil.rmtree(outUasDirectory)
    os.makedirs(outUasDirectory)

    # create OUT DIRECTORIES_A,Y,X,i in CSV directory
    for sheet in UAS_sheets:
        os.makedirs(outUasDirectory + os.path.normpath('/') + sheet)

    # create CSV files to OUT DIRECTORIES_A,Y,X
    UAS_xlsx_list = os.listdir(inUasDirectory)
    for file in UAS_xlsx_list:
        if file.endswith('.xlsx'):
            for sheet in UAS_sheets:
                path_xlsx = inUasDirectory + os.path.normpath('/') + file
                # path_csv = out_directory + file[0:-5] + '_' + sheet[0] + '.csv'
                path_csv = outUasDirectory + os.path.normpath('/') + sheet + os.path.normpath('/') + file[
                                                                                                     0:-5] + '_' + sheet[0] + '.csv'
                df = pd.read_excel(path_xlsx, sheet, header=None)
                df.to_csv(path_csv, header=None, index=None)
    print('UAS: xlsx2csv done')

    UAS_csv_list_0 = os.listdir(outUasDirectory + os.path.normpath('/') + UAS_sheets[0] + os.path.normpath('/'))
    for file in UAS_csv_list_0:
        writing_row_marker = 0
        row_number = 0
        UAS_project_name = 'null'
        UAS_created_name = 'null'
        UAS_user_name = 'null'
        UAS_sample_name = 'null'
        # create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        UAS_sample_raw_dict = {marker: [] for marker in IlluminaMarkers}
        # UAS_row_selected = []

        # reading csv files row by row
        for row in csv.reader(
                open(outUasDirectory + os.path.normpath('/') + UAS_sheets[0] + os.path.normpath('/') + file),
                delimiter=','):
            row_number += 1
            # print (row_number, '', row [0], row [1])
            if row_number == 3 and row[0] == 'Project':
                UAS_project_name = row[1]
                print(UAS_project_name)
            if row_number == 4 and row[0] == 'Created':
                UAS_created_name = row[1]
            if row_number == 5 and row[0] == 'User':
                UAS_user_name = row[1]
            if row_number == 13 and row[0] == 'Sample':
                writing_row_marker = 1

            # reading data for samples
            if row_number > 13 and writing_row_marker == 1:

                # first sample starts
                if UAS_sample_name != row[0] and UAS_sample_name == 'null':
                    UAS_Auto_STR_Head[row[0]] = {'Project': UAS_project_name, 'Analysis': row[1],
                                                 'Run': 'user ' + UAS_user_name, 'Gender': 'null',
                                                 'Created': UAS_created_name}
                # new sample starts - second, third, ... not first - new empty dictionary is created
                if UAS_sample_name != row[0] and UAS_sample_name != 'null':
                    UAS_Auto_STR_Head[row[0]] = {'Project': UAS_project_name, 'Analysis': row[1],
                                                 'Run': 'user ' + UAS_user_name, 'Gender': 'null',
                                                 'Created': UAS_created_name}
                    # create empty raw_dictionary and dictionary with homo or hetero alleles for markers
                    UAS_sample_raw_dict = {marker: [] for marker in IlluminaMarkers}

                # created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                UAS_sample_name = row[0]
                if row[2] in IlluminaMarkers:
                    UAS_row_selected = [row[2], row[3], 'N/A', row[4], row[5]]
                    UAS_sample_raw_dict[row[2]].append(UAS_row_selected)
                    UAS_AutoSTRData_raw[UAS_sample_name] = UAS_sample_raw_dict
    UasRow = [UAS_Auto_STR_Head, UAS_AutoSTRData_raw]
    return UasRow


def createGeneMarkerRaw(inGeneMarkerDirectory, PowerSeqMarkers, PowerSeq_project_info, GMThreshold):
    # GeneMarker
    GM_csv_list_0 = os.listdir(inGeneMarkerDirectory)
    GM_Auto_STR_Head = {}
    GM_AutoSTRData_raw = {}
    ColumnReportIsPresent = False
    GM_sample_name = 'null'
    for file in GM_csv_list_0:
        if file.endswith('.csv'):
            # GM_writing_row_marker = 0
            # GM_head_raw = []
            GM_row_number = 0

            GM_sample_name = 'null'
            # create empty raw_dictionary and dictionary with homo or hetero alleles for markers
            GM_sample_raw_dict = {marker: [] for marker in PowerSeqMarkers}
            # GM_sample_dict = {marker: [] for marker in PowerSeqMarkers}
            # GM_sorted_sample_dict = {marker: [] for marker in PowerSeqMarkers}
            # GM_row_selected = []
            ColumnReportIsPresent = False
            # reading csv files row by row
            for row in csv.reader(open(inGeneMarkerDirectory + os.path.normpath('/') + file), delimiter=','):

                GM_row_number += 1
                # print (row_number, '', row [0], row [1])
                if GM_row_number == 1 and row[0] != '#Program: GeneMarkerHTS':
                    print('wrong csv file or format', file)
                    break
                if GM_row_number == 3 and row[0].startswith('#Sample: '):
                    if len(row[0][9:].split("_")[0]) >= 0:
                        GM_sample_name = row[0][9:].split("_")[0]

                        GM_Auto_STR_Head[GM_sample_name] = {'Project': PowerSeq_project_info[0],
                                                            'Analysis': PowerSeq_project_info[2],
                                                            'Run': 'user ' + PowerSeq_project_info[3], 'Gender': 'null',
                                                            'Created': PowerSeq_project_info[1]}
                if GM_row_number == 4 and row[2] == 'Report':
                    ColumnReportIsPresent = True

                # reading data for samples
                if GM_row_number >= 6:

                    # created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, ???Reads_ReverseSeq???, Reads_ForwardSeq, Repeat Sequence]}
                    if row[0] in PowerSeqMarkers:
                        if ColumnReportIsPresent:
                            if len(row) > 18:
                                if row[9].isnumeric():
                                    GM_row_selected = [row[0], row[1], row[10], row[9], row[18]]
                                    # print(GM_row_selected)
                                    if int(GM_row_selected[3]) >= GMThreshold:
                                        GM_sample_raw_dict[row[0]].append(GM_row_selected)
                                        GM_AutoSTRData_raw[GM_sample_name] = GM_sample_raw_dict

                            else:
                                print('shortRow', GM_sample_name, row[0])
                        else:
                            if len(row) > 17:
                                if row[8].isnumeric():
                                    GM_row_selected = [row[0], row[1], row[9], row[8], row[17]]
                                    # print(GM_row_selected)
                                    if int(GM_row_selected[3]) >= GMThreshold:
                                        GM_sample_raw_dict[row[0]].append(GM_row_selected)
                                        GM_AutoSTRData_raw[GM_sample_name] = GM_sample_raw_dict
                            else:
                                print('shortRow', GM_sample_name, row[0])
        if ColumnReportIsPresent:
            print(GM_sample_name, 'columnReportIsPresent')

    GeneMarkerRaw = [GM_Auto_STR_Head, GM_AutoSTRData_raw]
    return GeneMarkerRaw


def createSTRaitRazorRaw(inSTRaitRazorDirectory, PowerSeqMarkers, PowerSeq_project_info):
    # STRaitRazor
    SR_csv_list_0 = os.listdir(inSTRaitRazorDirectory)
    SR_Auto_STR_Head = {}
    SR_AutoSTRData_raw = {}
    SR_sample_name = 'null'
    SR_sample_raw_dict = {marker: [] for marker in PowerSeqMarkers}
    for file in SR_csv_list_0:
        if file.endswith('.csv'):
            for row in csv.reader(open(inSTRaitRazorDirectory + os.path.normpath('/') + file), delimiter=','):
                if row[0].split("_")[3] == 'R1':
                    # new sample starts - new empty dictionary is created
                    if SR_sample_name != row[0].split("_")[0]:
                        SR_Auto_STR_Head[row[0].split("_")[0]] = {'Project': PowerSeq_project_info[0],
                                                                  'Analysis': PowerSeq_project_info[2],
                                                                  'Run': 'user ' + PowerSeq_project_info[3],
                                                                  'Gender': 'null', 'Created': PowerSeq_project_info[1]}
                        SR_sample_raw_dict = {marker: [] for marker in PowerSeqMarkers}
                    # created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                    SR_sample_name = row[0].split("_")[0]
                    if row[1] in PowerSeqMarkers or row[1] == 'DYS385':
                        if row[1] == 'DYS385':
                            Locus = 'DYS385a/b'
                        else:
                            Locus = row[1]

                        SR_row_selected = [Locus, row[4], 'N/A', row[2], row[3]]
                        SR_sample_raw_dict[Locus].append(SR_row_selected)
                        SR_AutoSTRData_raw[SR_sample_name] = SR_sample_raw_dict
    STRaitRazorRaw = [SR_Auto_STR_Head, SR_AutoSTRData_raw]
    return STRaitRazorRaw


def getAuto_STR_Head(rawAnalysis):
    head = rawAnalysis[0]
    return head


def getAuto_STR_Data_raw(rawAnalysis):
    data = rawAnalysis[1]
    return data


def selectTrueAlleles(markers, samples, AutoSTR_Data_raw, heterozygozity_dict):
    AutoSTR_Data = {}
    for sample in samples:
        # create empty raw_dictionary and dictionary with homo or hetero alleles for markers
        sample_dict = {marker: [] for marker in markers}
        sorted_sample_dict = {marker: [] for marker in markers}

        # select true readings to sample_dict according to marker's heterozygozity
        for marker in markers:
            sorted_markers = sorted(AutoSTR_Data_raw[sample][marker], key=lambda x: int(x[3]), reverse=True)

            if len(sorted_markers) == 1:
                sample_dict[marker].append(sorted_markers[0])
                sample_dict[marker].append(sorted_markers[0])

            if len(sorted_markers) > 1:
                if int(sorted_markers[1][3]) / int(sorted_markers[0][3]) > float(heterozygozity_dict[marker]):
                    sample_dict[marker].append(sorted_markers[0])
                    sample_dict[marker].append(sorted_markers[1])

                else:
                    sample_dict[marker].append(sorted_markers[0])
                    sample_dict[marker].append(sorted_markers[0])

            if len(sorted_markers) == 0:
                sample_dict[marker].append(['Null', '0', 'Null', '0', 'Null'])
                sample_dict[marker].append(['Null', '0', 'Null', '0', 'Null'])

            if sample_dict['Amelogenin'][0][1] == 'chrX':
                sample_dict['Amelogenin'][0][1] = '0'
            if sample_dict['Amelogenin'][0][1] == 'chrY' or sample_dict['Amelogenin'][0][1] == '1':
                sample_dict['Amelogenin'][0][1] = '6'

            if sample_dict['Amelogenin'][1][1] == 'chrX':
                sample_dict['Amelogenin'][1][1] = '0'
            if sample_dict['Amelogenin'][1][1] == 'chrY' or sample_dict['Amelogenin'][1][1] == '1':
                sample_dict['Amelogenin'][1][1] = '6'

            if sample_dict[marker][0][1] == 'Unknown':
                sample_dict[marker][0][1] = '999'
            if sample_dict[marker][1][1] == 'Unknown':
                sample_dict[marker][1][1] = '999'

            sorted_sample_dict[marker] = sorted(sample_dict[marker], key=lambda x: float(x[1]))

        AutoSTR_Data[sample] = sorted_sample_dict
    return AutoSTR_Data


def markerAlleleList(AutoSTR_Data, sample, marker):
    AlleleList = [AutoSTR_Data[sample][marker][0][1], AutoSTR_Data[sample][marker][1][1]]
    return AlleleList


def markerNoReadsList(AutoSTR_Data, sample, marker):
    NoReadsList = [AutoSTR_Data[sample][marker][0][3], AutoSTR_Data[sample][marker][1][3]]
    return NoReadsList


reports_list = ['STRaitRazor', 'GeneMarker', 'UAS']
# reports_list = ['STRaitRazor',  'UAS']

home = getHome()

# setup of input directories
in_UAS_directory = home + os.path.normpath('/NGS_forensic_database/xlsx_detail_reports')
out_UAS_directory = home + os.path.normpath('/NGS_forensic_database/csv_output')
in_GeneMarker_directory = home + os.path.normpath('/NGS_forensic_database/GeneMarker_reports')
in_STRaitRazor_directory = home + os.path.normpath('/NGS_forensic_database/STRaitRazor_reports')
in_project_info = home + os.path.normpath('/NGS_forensic_database/GeneMarker_reports/project_info.txt')

# setup of markers
PowerSeq_markers_heterozygozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5, 'TPOX': 0.5, 'D2S441': 0.5, 'D2S1338': 0.5,
                                        'D3S1358': 0.5, 'FGA': 0.5, 'D5S818': 0.5, 'CSF1PO': 0.5, 'D7S820': 0.5,
                                        'D8S1179': 0.5, 'D10S1248': 0.5, 'TH01': 0.5, 'vWA': 0.5, 'D12S391': 0.5,
                                        'D13S317': 0.5, 'PentaE': 0.5, 'D16S539': 0.5, 'D18S51': 0.5, 'D19S433': 0.5,
                                        'D21S11': 0.5, 'PentaD': 0.5, 'D22S1045': 0.5, 'DYS19': 0.5, 'DYS385a/b': 0.5}
PowerSeq_markers = list(PowerSeq_markers_heterozygozity_dict.keys())
Illumina_markers_heterozygozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5, 'TPOX': 0.5, 'D2S441': 0.5, 'D2S1338': 0.5,
                                        'D3S1358': 0.5, 'D4S2408': 0.5, 'FGA': 0.5, 'D5S818': 0.5, 'CSF1PO': 0.5,
                                        'D6S1043': 0.5, 'D7S820': 0.5, 'D8S1179': 0.5, 'D9S1122': 0.5, 'D10S1248': 0.5,
                                        'TH01': 0.5, 'vWA': 0.5, 'D12S391': 0.5, 'D13S317': 0.5, 'PentaE': 0.5,
                                        'D16S539': 0.5, 'D17S1301': 0.5, 'D18S51': 0.5, 'D19S433': 0.5, 'D20S482': 0.5,
                                        'D21S11': 0.5, 'PentaD': 0.5, 'D22S1045': 0.5}
Illumina_markers = list(Illumina_markers_heterozygozity_dict.keys())

GM_threshold = 10
# create empty directories
# UAS_Auto_STR_Data = {}
# UAS_Auto_STR_Data_unsorted = {}
# GM_Auto_STR_Data = {}
# GM_Auto_STR_Data_unsorted = {}
# SR_Auto_STR_Data = {}
# SR_Auto_STR_Data_unsorted = {}

PowerSeq_ProjectInfo = createPowerSeqProjectInfo(in_project_info)

if 'UAS' in reports_list:
    UAS_raw = createUasRaw(in_UAS_directory, out_UAS_directory, Illumina_markers)
    UAS_Head = getAuto_STR_Head(UAS_raw)
    UAS_Auto_STR_Data_raw = getAuto_STR_Data_raw(UAS_raw)
    UAS_samples = list(UAS_Auto_STR_Data_raw.keys())
    UAS_Auto_STR_Data = selectTrueAlleles(Illumina_markers, UAS_samples, UAS_Auto_STR_Data_raw,
                                          Illumina_markers_heterozygozity_dict)
    print('UAS samples: ', UAS_samples)
if 'GeneMarker' in reports_list:
    GM_raw = createGeneMarkerRaw(in_GeneMarker_directory, PowerSeq_markers, PowerSeq_ProjectInfo, GM_threshold)
    GM_Head = getAuto_STR_Head(GM_raw)
    GM_Auto_STR_Data_raw = getAuto_STR_Data_raw(GM_raw)
    GM_samples = list(GM_Auto_STR_Data_raw.keys())
    GM_Auto_STR_Data = selectTrueAlleles(PowerSeq_markers, GM_samples, GM_Auto_STR_Data_raw,
                                         PowerSeq_markers_heterozygozity_dict)
    print('GeneMarker samples: ', GM_samples)
if 'STRaitRazor' in reports_list:
    SR_raw = createSTRaitRazorRaw(in_STRaitRazor_directory, PowerSeq_markers, PowerSeq_ProjectInfo)
    SR_Head = getAuto_STR_Head(SR_raw)
    SR_Auto_STR_Data_raw = getAuto_STR_Data_raw(SR_raw)
    SR_samples = list(SR_Auto_STR_Data_raw.keys())
    SR_Auto_STR_Data = selectTrueAlleles(PowerSeq_markers, SR_samples, SR_Auto_STR_Data_raw,
                                         PowerSeq_markers_heterozygozity_dict)
    print('STRaitRazor samples: ', SR_samples)

print('GM_project_name:', PowerSeq_ProjectInfo[0])
print('GM_created_name:', PowerSeq_ProjectInfo[1])
print('GM_analysis_name:', PowerSeq_ProjectInfo[2])
print('GM_user_name:', PowerSeq_ProjectInfo[3])

# for sample in GM_samples:
#     if sample in SR_samples:
#         for marker in PowerSeq_markers:
#             GM_AlleleList = markerAlleleList(GM_Auto_STR_Data, sample, marker)
#             SR_AlleleList = markerAlleleList(SR_Auto_STR_Data, sample, marker)
#             GM_NoReadsList = markerNoReadsList(GM_Auto_STR_Data, sample, marker)
#             SR_NoReadsList = markerNoReadsList(SR_Auto_STR_Data, sample, marker)
#             if GM_AlleleList != SR_AlleleList:
#                 print (sample, marker, 'GeneMarker: ', GM_AlleleList,GM_NoReadsList, 'STRaitRazor: ', SR_AlleleList, SR_NoReadsList)
#     else:
#         print (sample, 'is not in STRaitRazor samples')

# print (UAS_Auto_STR_Data)


# print(GM_Auto_STR_Data_raw['1-01'])
