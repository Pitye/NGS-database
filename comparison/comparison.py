import os
import shutil
import csv
import pandas as pd
import sys
import xml.etree.ElementTree as etree


def systemLinux():
    if sys.platform == 'linux':
        return True


def getHome():
    if systemLinux():
        print('system is Linux')
        homePath = os.path.normpath('/home/ku/')
    else:
        print('system is Windows')
        homePath = os.path.normpath('C:/')
    return homePath


def renameSamples(dict2rename, renameDict):
    newDict = {}
    for sample in list(dict2rename.keys()):
        if sample in list(renameDict.keys()):
            newDict[renameDict[sample]] = dict2rename[sample]
            del dict2rename[sample]
            print('SampleName', sample, 'changed to', renameDict[sample])
    updatedDict = merge2dicts(newDict, dict2rename)
    return updatedDict


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


def createUasRaw(inUasDirectory, outUasDirectory, IlluminaMarkers, version):
    UAS_Auto_STR_Head = {}
    UAS_AutoSTRData_raw = {}
    #UAS_sheets = ['Autosomal STR Coverage', 'X STR Coverage', 'Y STR Coverage', 'iSNP Coverage']
    UAS_sheets = ['Autosomal STR Coverage']
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
                                                                                                     0:-5] + '_' + \
                           sheet[0] + '.csv'
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
            if version == 'one report':
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
            if version == 'separate reports for samples':
                if row_number == 4 and row[0] == 'Project':
                    UAS_project_name = row[1]
                    print(UAS_project_name)
                if row_number == 5 and row[0] == 'Sample':
                    UAS_sample_name = row[1]
                if row_number == 7 and row[0] == 'Created':
                    UAS_created_name = row[1]
                if row_number == 8 and row[0] == 'User':
                    UAS_user_name = row[1]

                    UAS_Auto_STR_Head[row[0]] = {'Project': UAS_project_name, 'Analysis': row[1],
                                                'Run': 'user ' + UAS_user_name, 'Gender': 'null',
                                                'Created': UAS_created_name}

                    # created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                if row_number > 15:
                    if row[0] in IlluminaMarkers:
                        UAS_row_selected = [row[0], row[1], 'N/A', row[2], row[3]]
                        UAS_sample_raw_dict[row[0]].append(UAS_row_selected)
                        UAS_AutoSTRData_raw[UAS_sample_name] = UAS_sample_raw_dict
    if renameSamples_UAS:
        print('renamed samples UAS:')
        UAS_AutoSTRData_raw = renameSamples(UAS_AutoSTRData_raw, renameDict_UAS)
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
                        GM_sample_name = row[0][9:].split("_")[0].upper()

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
    if renameSamples_GM:
        print('renamed samples GeneMarker:')
        GM_AutoSTRData_raw = renameSamples(GM_AutoSTRData_raw, renameDict_GM)
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
                if len(row) > 0 and len(row[0].split("_")) > 3:
                    if row[0].split("_")[3] == 'R2':
                        # new sample starts - new empty dictionary is created
                        if SR_sample_name != row[0].split("_")[0].upper():
                            SR_Auto_STR_Head[row[0].split("_")[0].upper()] = {'Project': PowerSeq_project_info[0],
                                                                              'Analysis': PowerSeq_project_info[2],
                                                                              'Run': 'user ' + PowerSeq_project_info[3],
                                                                              'Gender': 'null',
                                                                              'Created': PowerSeq_project_info[1]}
                            SR_sample_raw_dict = {marker: [] for marker in PowerSeqMarkers}
                        # created Auto_STR_Data_raw[sample_name] = {Locus:[Locus, Allele Name, Typed Allele?, Reads, Repeat Sequence]}
                        SR_sample_name = row[0].split("_")[0].upper()

                        #if row[1] in PowerSeqMarkers or row[1] == 'DYS385':
                        #    if row[1] == 'DYS385':
                        #        Locus = 'DYS385a/b'
                        #    else:
                        #        Locus = row[1]
                        if row[1] in PowerSeqMarkers:
                            Locus = row[1]

                            SR_row_selected = [Locus, row[4], 'N/A', row[2], row[3]]
                            SR_sample_raw_dict[Locus].append(SR_row_selected)
                            SR_AutoSTRData_raw[SR_sample_name] = SR_sample_raw_dict
    if renameSamples_SR:
        print('renamed samples STRaitRazor:')
        SR_AutoSTRData_raw = renameSamples(SR_AutoSTRData_raw, renameDict_SR)
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

            if sample_dict['Amelogenin'][0][1] == 'chrX' or sample_dict['Amelogenin'][0][1] =='X':
                sample_dict['Amelogenin'][0][1] = '0'
            if sample_dict['Amelogenin'][0][1] == 'chrY' or sample_dict['Amelogenin'][0][1] == '1' or sample_dict['Amelogenin'][0][1] =='Y':
                sample_dict['Amelogenin'][0][1] = '6'

            if sample_dict['Amelogenin'][1][1] == 'chrX' or sample_dict['Amelogenin'][1][1] =='X':
                sample_dict['Amelogenin'][1][1] = '0'
            if sample_dict['Amelogenin'][1][1] == 'chrY' or sample_dict['Amelogenin'][1][1] == '1' or sample_dict['Amelogenin'][1][1] =='Y':
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
    if marker == 'Amelogenin':
        if AlleleList == ['0', '6']:
            AlleleList = ['X', 'Y']
        if AlleleList == ['0', '0']:
            AlleleList = ['X', 'X']
    return AlleleList


def markerNoReadsList(AutoSTR_Data, sample, marker):
    NoReadsList = [AutoSTR_Data[sample][marker][0][3], AutoSTR_Data[sample][marker][1][3]]
    return NoReadsList


def mergeLists2unique(list1, list2, list3):
    mergedList = list1 + list2 + list3
    uniqueList = list(dict.fromkeys(mergedList))
    return uniqueList


def commonMembers(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return list(set1.intersection(set2))


def merge2dicts(x, y):
    z = x.copy()  # start with x's keys and values
    z.update(y)  # modifies z with y's keys and values & returns None
    return z


def readCsvCE(CSV_CE_dir):
    CSV_CE_list = os.listdir(CSV_CE_dir)
    CSV_data_CE = {}
    if len(CSV_CE_list) > 0:
        # CSV_data_CE = {}
        for file in CSV_CE_list:
            if file.endswith('.csv'):
                genotype = {}

                for row in csv.reader(open(CSV_CE_dir + os.path.normpath('/') + file), delimiter=','):
                    specimen_id = row[9]
                    locus_name = row[79].replace(' ', '')
                    allele_values = row[80].split(", ")
                    if not locus_name.startswith('DYS') and not locus_name.startswith('YGAT') and len(
                            allele_values) == 1:
                        allele_values.append(allele_values[0])
                    genotype[locus_name] = allele_values

                CSV_data_CE[specimen_id] = genotype
    else:
        print('no CE csv files in', CSV_CE_dir)
    return CSV_data_CE


def readXmlCE_file(path):
    tree = etree.parse(path)
    namespaces = {'ns': 'urn:CODISImportFile-schema'}
    root = tree.getroot()

    XML_data = {}
    for specimen_element in root.findall("ns:SPECIMEN", namespaces):
        # specimen_id = specimen_element.find("ns:SPECIMENID", namespaces).text.split("-")[1].upper()
        specimen_id = specimen_element.find("ns:SPECIMENID", namespaces).text.upper()
        genotype = {}
        for locus_element in specimen_element.findall("ns:LOCUS", namespaces):
            locus_name = locus_element.find("ns:LOCUSNAME", namespaces).text.replace(' ', '')
            allele_values = []

            for allele_value_element in locus_element.findall("*/ns:ALLELEVALUE", namespaces):
                allele_values.append(allele_value_element.text)

            if not locus_name.startswith('DYS') and not locus_name.startswith('YGAT') and len(allele_values) == 1:
                allele_values.append(allele_values[0])

            genotype[locus_name] = allele_values
            # dictionary[new_key] = dictionary.pop(old_key)
            if 'YGATAH4' in locus_name:
                genotype['Y-GATA-H4'] = genotype.pop('YGATAH4')
            if 'DYS385' in locus_name:
                genotype['DYS385a-b'] = genotype.pop('DYS385')

        XML_data[specimen_id] = genotype
    return XML_data


def readXmlCE(XML_CE_dir):
    XML_CE_list = os.listdir(XML_CE_dir)
    XML_CE_data = {}
    if len(XML_CE_list) > 0:
        for file in XML_CE_list:
            if file.endswith('.xml'):
                path_xml = XML_CE_dir + os.path.normpath('/') + file
                XML_CE_data_part = readXmlCE_file(path_xml)
                XML_CE_data = merge2dicts(XML_CE_data, XML_CE_data_part)
    else:
        print('no CE xml files in', XML_CE_dir)
    return XML_CE_data


def removeFromList(originalList, removingList):
    for element in removingList:
        newList = originalList
        newList.remove(element)
    return newList


def sortListByValue(str_list):
    sorted_list1=[]
    if 'X' in str_list or 'Y' in str_list:
        sorted_list_str = str_list
    else:
        num_list = [float(i) for i in str_list]
        sorted_list = sorted(num_list)
        for i in sorted_list:
            if int(i) == float(i):
                sorted_list1.append(int(i))
            else:
                sorted_list1.append(float(i))

        sorted_list_str = [str(i) for i in sorted_list1]
    return sorted_list_str

# reports_list = ['STRaitRazor', 'GeneMarker', 'UAS']
# reports_list = ['STRaitRazor',  'UAS']
#reports_list = ['GeneMarker', 'UAS']
#reports_list = ['STRaitRazor', 'GeneMarker']
# reports_list = ['STRaitRazor']
# reports_list = ['GeneMarker']
reports_list = ['UAS']
reports_list.sort()

#renameDict_SR = {'5-25': '5-49', '5-33': 'JO23', '5-34': 'JO24', '5-36': 'MJ2', '5-37': 'MJ3', '5-38': 'MJ4',
               # '5-39': 'MJ5', '5-49': '5-25', 'JO23': '5-33', 'JO24': '5-34', 'MJ1': '5-35', 'MJ2': '5-36','MJ3': '5-37', 'MJ4': '5-38', 'MJ5': '5-39'}
# renameDict_SR = {'3499':'3491', '3500':'3492', '3501':'3493', '3502':'3494', '3503':'3495', '3504':'3496', '3505':'3497', '3506':'3498', '3507':'3499', '3508':'3500', '3509':'3501', '3510':'3502', '3511':'3503', '3512':'3504', '3513':'3505', '3514':'3506', '3515':'3507', '3516':'3508', '3517':'3509', '3518':'3510', '3519':'3511', '3520':'3512', '3521':'3513', '3522':'3514', '3523':'3515', '3524':'3516', '3525':'3517', '3526':'3518', '3527':'3519', '3528':'3520', '3529':'3521', '3530':'3522', '3531':'3523', '3532':'3524', '3533':'3525', '3534':'3526', '3535':'3527', '3536':'3528', '3537':'3529', '3539':'3530', '3540':'3531', '3541':'3532', '8349':'3533', '8350':'3534', '8351':'3535', '8352':'3536', '8353':'3537', '8354':'3539', '8355':'3540', '8356':'3541', '8357':'8349', '8358':'8350', '8359':'8351', '8360':'8352', '8361':'8353', '8362':'8354', '8363':'8355', '8364':'8356', '8365':'8357', '8366':'8358', '8367':'8359', '8368':'8360', '8369':'8361', '8370':'8362'}
# renameDict_GM = {'3499':'3491', '3500':'3492', '3501':'3493', '3502':'3494', '3503':'3495', '3504':'3496', '3505':'3497', '3506':'3498', '3507':'3499', '3508':'3500', '3509':'3501', '3510':'3502', '3511':'3503', '3512':'3504', '3513':'3505', '3514':'3506', '3515':'3507', '3516':'3508', '3517':'3509', '3518':'3510', '3519':'3511', '3520':'3512', '3521':'3513', '3522':'3514', '3523':'3515', '3524':'3516', '3525':'3517', '3526':'3518', '3527':'3519', '3528':'3520', '3529':'3521', '3530':'3522', '3531':'3523', '3532':'3524', '3533':'3525', '3534':'3526', '3535':'3527', '3536':'3528', '3537':'3529', '3539':'3530', '3540':'3531', '3541':'3532', '8349':'3533', '8350':'3534', '8351':'3535', '8352':'3536', '8353':'3537', '8354':'3539', '8355':'3540', '8356':'3541', '8357':'8349', '8358':'8350', '8359':'8351', '8360':'8352', '8361':'8353', '8362':'8354', '8363':'8355', '8364':'8356', '8365':'8357', '8366':'8358', '8367':'8359', '8368':'8360', '8369':'8361', '8370':'8362'}
#renameDict_SR = {'R701-A501': '1_38', 'R701-A502': '1_39', 'R701-A503': '1_40', 'R701-A504': '1_41', 'R701-A506': '1_42', 'R701-A507': '1_43', 'R701-A505': '1_44', 'R701-A508': '1_45', 'R702-A501': '1_46', 'R702-A502': '1_47', 'R702-A503': '1_48', 'R702-A504': '5_09', 'R702-A506': '5_20', 'R702-A507': '5_22', 'R702-A505': '5_23', 'R702-A508': '5_24', 'R703-A501': '5_25', 'R703-A502': '5_33', 'R703-A503': '5_34', 'R703-A504': '5_35', 'R703-A506': '5_36', 'R703-A507': '5_37', 'R703-A505': '5_38', 'R703-A508': '5_39', 'R704-A501': '5_49', 'R704-A502': 'JO23', 'R704-A503': 'JO24', 'R704-A504': 'MJ1', 'R704-A506': 'MJ2', 'R704-A507': 'MJ3', 'R704-A505': 'MJ4', 'R704-A508': 'MJ5', 'R705-A501': 'MJ6', 'R705-A502': 'MJ7', 'R705-A503': 'MJ8', 'R705-A504': 'MJ9', 'R705-A506': '1_01', 'R705-A507': '1_02', 'R705-A505': '1_08', 'R705-A508': '1_14', 'R706-A501': '1_18', 'R706-A502': '1_20', 'R706-A503': '1_23', 'R706-A504': '1_25', 'R706-A506': '1_31', 'R706-A507': '1_34', 'R706-A505': 'PK', 'R706-A508': 'NK', 'R707-A501': 'kost 5158', 'R707-A502': '4_01'}
# renameDict_GM = {'5-25': '5-49', '5-33': 'JO23', '5-34': 'JO24', '5-36': 'MJ2', '5-37': 'MJ3', '5-38': 'MJ4',
#                '5-39': 'MJ5', '5-49': '5-25', 'JO23': '5-33', 'JO24': '5-34', 'MJ1': '5-35', 'MJ2': '5-36',
#                 'MJ3': '5-37', 'MJ4': '5-38', 'MJ5': '5-39'}
renameDict_UAS = {'7_A':'7A', '8_A':'8A', '2_B':'2B', '2_A':'2A', '6_B':'6B', '6_A':'6A', '1_B':'1B', '9_A':'9A', '1_A':'1A', '3_A':'3A', '4_A':'4A', '4_B':'4B', '5_A':'5A', '8_B':'8B', '3_B':'3B', '9_B':'9B', '5_B':'5B', '7_B':'7B'}
samples_compar_dict = {'1A' :'1B', '2A':'2B', '3A':'3B', '4A':'4B', '5A':'5B', '6A':'6B', '7A':'7B', '8A':'8B', '9A':'9B', '10A':'10B', '11A':'11B', '12A':'12B', '13A':'13B', '14A':'14B', '15A':'15B', '16A':'16B', '17A':'17B', '18A':'18B', '19A':'19B', '20A':'20B', '21A':'21B', '22A':'22B', '23A':'23B'}
renameSamples_SR = False
renameSamples_GM = False
renameSamples_UAS = True

printMatches = True
CE_comparison = True
home = getHome()

# UAS version - 'one report' or 'separate reports for samples'
version = 'separate reports for samples'

# setup of input directories
in_UAS_directory = home + os.path.normpath('/NGS_forensic_database/xlsx_detail_reports')
out_UAS_directory = home + os.path.normpath('/NGS_forensic_database/csv_output')
in_GeneMarker_directory = home + os.path.normpath('/NGS_forensic_database/GeneMarker_reports')
in_STRaitRazor_directory = home + os.path.normpath('/NGS_forensic_database/STRaitRazor_reports')
in_project_info = home + os.path.normpath('/NGS_forensic_database/STRaitRazor_reports/project_info.txt')
in_CE_csv_directory = home + os.path.normpath('/NGS_forensic_database/csv_CE')
in_CE_xml_directory = home + os.path.normpath('/NGS_forensic_database/xml_CE')

# setup of markers
PowerSeq_markers_heterozygozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5, 'TPOX': 0.5, 'D2S441': 0.5, 'D2S1338': 0.5,
                                        'D3S1358': 0.5, 'D4S2408': 0.5, 'FGA': 0.5, 'D5S818': 0.5, 'CSF1PO': 0.5,
                                        'D6S1043': 0.5, 'D7S820': 0.5, 'D8S1179': 0.5, 'D9S1122': 0.5, 'D10S1248': 0.5,
                                        'TH01': 0.5, 'vWA': 0.5, 'D12S391': 0.5, 'D13S317': 0.5, 'PentaE': 0.5,
                                        'D16S539': 0.5, 'D17S1301': 0.5, 'D18S51': 0.5, 'D19S433': 0.5, 'D20S482': 0.5,
                                        'D21S11': 0.5, 'PentaD': 0.5, 'D22S1045': 0.5}
#PowerSeq_markers_heterozygozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5, 'TPOX': 0.5, 'D2S441': 0.5, 'D2S1338': 0.5,
#                                        'D3S1358': 0.5, 'FGA': 0.5, 'D5S818': 0.5, 'CSF1PO': 0.5, 'D7S820': 0.5,
#                                        'D8S1179': 0.5, 'D10S1248': 0.5, 'TH01': 0.5, 'vWA': 0.5, 'D12S391': 0.5,
#                                        'D13S317': 0.5, 'PentaE': 0.5, 'D16S539': 0.5, 'D18S51': 0.5, 'D19S433': 0.5,
#                                        'D21S11': 0.5, 'PentaD': 0.5, 'D22S1045': 0.5, 'DYS19': 0.5, 'DYS385a/b': 0.5}
PowerSeq_markers = list(PowerSeq_markers_heterozygozity_dict.keys())
Illumina_markers_heterozygozity_dict = {'Amelogenin': 0.5, 'D1S1656': 0.5, 'TPOX': 0.5, 'D2S441': 0.5, 'D2S1338': 0.5,
                                        'D3S1358': 0.5, 'D4S2408': 0.5, 'FGA': 0.5, 'D5S818': 0.5, 'CSF1PO': 0.5,
                                        'D6S1043': 0.5, 'D7S820': 0.5, 'D8S1179': 0.5, 'D9S1122': 0.5, 'D10S1248': 0.5,
                                        'TH01': 0.5, 'vWA': 0.5, 'D12S391': 0.5, 'D13S317': 0.5, 'PentaE': 0.5,
                                        'D16S539': 0.5, 'D17S1301': 0.5, 'D18S51': 0.5, 'D19S433': 0.5, 'D20S482': 0.5,
                                        'D21S11': 0.5, 'PentaD': 0.5, 'D22S1045': 0.5}
#Illumina_markers_heterozygozity_dict = {'Amelogenin': 0.3, 'D1S1656': 0.3, 'TPOX': 0.3, 'D2S441': 0.3, 'D2S1338': 0.3,
#                                       'D3S1358': 0.3, 'D4S2408': 0.3, 'FGA': 0.3, 'D5S818': 0.3, 'CSF1PO': 0.3,
#                                       'D6S1043': 0.3, 'D7S820': 0.3, 'D8S1179': 0.3, 'D9S1122': 0.3, 'D10S1248': 0.3,
#                                       'TH01': 0.3, 'vWA': 0.3, 'D12S391': 0.3, 'D13S317': 0.3, 'PentaE': 0.3,
#                                       'D16S539': 0.3, 'D17S1301': 0.3, 'D18S51': 0.3, 'D19S433': 0.3, 'D20S482': 0.3,
#                                       'D21S11': 0.3, 'PentaD': 0.3, 'D22S1045': 0.3}
Illumina_markers = list(Illumina_markers_heterozygozity_dict.keys())

GM_threshold = 10
# create empty directories
# UAS_Auto_STR_Data = {}
# UAS_Auto_STR_Data_unsorted = {}
# GM_Auto_STR_Data = {}
# GM_Auto_STR_Data_unsorted = {}
# SR_Auto_STR_Data = {}
# SR_Auto_STR_Data_unsorted = {}

# create empty lists
UAS_samples = []
GM_samples = []
SR_samples = []

if 'GeneMarker' in reports_list or 'STRaitRazor' in reports_list:
    PowerSeq_ProjectInfo = createPowerSeqProjectInfo(in_project_info)

if 'UAS' in reports_list:
    UAS_raw = createUasRaw(in_UAS_directory, out_UAS_directory, Illumina_markers, version)
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
    # print for comparison in GM
    for sample in GM_samples:
        for marker in PowerSeq_markers:
            #print(sample, marker, markerAlleleList(GM_Auto_STR_Data, sample, marker), markerNoReadsList(GM_Auto_STR_Data, sample, marker))
            print(sample, marker, markerAlleleList(GM_Auto_STR_Data, sample, marker))

    print('GeneMarker samples: ', GM_samples)
if 'STRaitRazor' in reports_list:
    SR_raw = createSTRaitRazorRaw(in_STRaitRazor_directory, PowerSeq_markers, PowerSeq_ProjectInfo)
    SR_Head = getAuto_STR_Head(SR_raw)
    SR_Auto_STR_Data_raw = getAuto_STR_Data_raw(SR_raw)
    SR_samples = list(SR_Auto_STR_Data_raw.keys())
    SR_Auto_STR_Data = selectTrueAlleles(PowerSeq_markers, SR_samples, SR_Auto_STR_Data_raw,
                                         PowerSeq_markers_heterozygozity_dict)
    print('STRaitRazor samples: ', SR_samples)


if CE_comparison:
    CE_csv_dict = readCsvCE(in_CE_csv_directory)
    CE_xml_dict = readXmlCE(in_CE_xml_directory)
    CE_data = merge2dicts(CE_csv_dict, CE_xml_dict)
    CE_samples = list(CE_data.keys())
    print('CE_samples: ', CE_samples)

if 'GeneMarker' in reports_list or 'STRaitRazor' in reports_list:
    print('GM_project_name:', PowerSeq_ProjectInfo[0])
    print('GM_created_name:', PowerSeq_ProjectInfo[1])
    print('GM_analysis_name:', PowerSeq_ProjectInfo[2])
    print('GM_user_name:', PowerSeq_ProjectInfo[3])

samples_all = mergeLists2unique(UAS_samples, GM_samples, SR_samples)
markers_all = mergeLists2unique(PowerSeq_markers, Illumina_markers, [])
# print('markers_all', markers_all)
markers_common = commonMembers(PowerSeq_markers, Illumina_markers)
# print('markers_common', markers_common)
markers_not_compared = removeFromList(markers_all, markers_common)


# markers_not_compared = markers_all - markers_common


def compare1analysisCE(analysisList, CE_compar):
    if CE_compar:
        if 'GeneMarker' in analysisList:
            samples0 = GM_samples
            data0 = GM_Auto_STR_Data
            markers = PowerSeq_markers
        if 'STRaitRazor' in analysisList:
            samples0 = SR_samples
            data0 = SR_Auto_STR_Data
            markers = PowerSeq_markers
        if 'UAS' in analysisList:
            samples0 = UAS_samples
            data0 = UAS_Auto_STR_Data
            markers = Illumina_markers
        for sample in samples0:
            if sample in CE_samples:
                markersMatch = []
                for marker in markers:
                    if marker in list(CE_data[sample].keys()):
                        CE_alleleList = sortListByValue(CE_data[sample][marker])
                        alleleList0 = markerAlleleList(data0, sample, marker)
                        noReadsList0 = markerNoReadsList(data0, sample, marker)
                        if alleleList0 != CE_alleleList:
                            print(sample, marker, analysisList[0] + ': ', alleleList0, noReadsList0, 'CE:',
                                  CE_alleleList)
                        if alleleList0 == CE_alleleList:
                            markersMatch.append(marker)
                if printMatches:
                    print(sample, "CE match: ", markersMatch, "\n")
            else:
                print(sample, "not in CE data\n")


def compare2analysis(analysisList, CE_compar):
    if 'UAS' not in analysisList:
        samples0 = GM_samples
        data0 = GM_Auto_STR_Data
        samples1 = SR_samples
        data1 = SR_Auto_STR_Data
        markers = PowerSeq_markers

    if 'GeneMarker' not in analysisList:
        samples0 = SR_samples
        data0 = SR_Auto_STR_Data
        samples1 = UAS_samples
        data1 = UAS_Auto_STR_Data
        markers = markers_common
        print('not compared markers', markers_not_compared)

    if 'STRaitRazor' not in analysisList:
        samples0 = GM_samples
        data0 = GM_Auto_STR_Data
        samples1 = UAS_samples
        data1 = UAS_Auto_STR_Data
        markers = markers_common
        print('not compared markers', markers_not_compared)

    for sample in samples_all:
        if sample not in samples0:
            print(sample, ' is not in ', analysisList[0])
        if sample not in samples1:
            print(sample, ' is not in ', analysisList[1])
        if sample in samples0 and sample in samples1:
            markersMatch = []
            for marker in markers:
                alleleList0 = markerAlleleList(data0, sample, marker)
                alleleList1 = markerAlleleList(data1, sample, marker)
                noReadsList0 = markerNoReadsList(data0, sample, marker)
                noReadsList1 = markerNoReadsList(data1, sample, marker)
                if CE_compar and sample in CE_samples and marker in list(CE_data[sample].keys()):
                    CE_alleleList = sortListByValue(CE_data[sample][marker])
                    if alleleList0 != alleleList1 or alleleList0 != CE_alleleList or alleleList1 != CE_alleleList:
                        print(sample, marker, analysisList[0] + ': ', alleleList0, noReadsList0,
                              analysisList[1] + ': ',
                              alleleList1, noReadsList1, 'CE:', CE_alleleList)
                    if alleleList0 == alleleList1 and alleleList0 == CE_alleleList and alleleList1 == CE_alleleList:
                        markersMatch.append(marker)

                else:
                    if alleleList0 != alleleList1:
                        print(sample, marker, analysisList[0] + ': ', alleleList0, noReadsList0, analysisList[1] + ': ',
                              alleleList1, noReadsList1, 'no CE data')
                    if alleleList0 == alleleList1:
                        markersMatch.append(marker)
            if printMatches:
                print(sample, analysisList[0], "-", analysisList[1], "match: ", markersMatch, "\n")


def compare3analysis(analysisList):
    if analysisList == ['GeneMarker', 'STRaitRazor', 'UAS']:
        samples0 = GM_samples
        data0 = GM_Auto_STR_Data
        samples1 = SR_samples
        data1 = SR_Auto_STR_Data
        samples2 = UAS_samples
        data2 = UAS_Auto_STR_Data
        markers = markers_common
    for sample in samples_all:
        if sample not in samples0:
            print(sample, ' is not in ', analysisList[0])
        if sample not in samples1:
            print(sample, ' is not in ', analysisList[1])
        if sample not in samples2:
            print(sample, ' is not in ', analysisList[2])

        if sample in samples0 and sample in samples1 and sample in samples2:
            for marker in markers:
                alleleList0 = markerAlleleList(data0, sample, marker)
                alleleList1 = markerAlleleList(data1, sample, marker)
                alleleList2 = markerAlleleList(data2, sample, marker)
                noReadsList0 = markerNoReadsList(data0, sample, marker)
                noReadsList1 = markerNoReadsList(data1, sample, marker)
                noReadsList2 = markerNoReadsList(data2, sample, marker)
                if alleleList0 != alleleList1 or alleleList1 != alleleList2 or alleleList0 != alleleList2:
                    print(sample, marker, analysisList[0] + ': ', alleleList0, noReadsList0, analysisList[1] + ': ',
                          alleleList1, noReadsList1, analysisList[2] + ': ', alleleList2, noReadsList2)


if len(reports_list) == 1:
    compare1analysisCE(reports_list, CE_comparison)
if len(reports_list) == 2:
    compare2analysis(reports_list, CE_comparison)
if len(reports_list) == 3:
    compare3analysis(reports_list)
print('comparison finished')

for s in UAS_samples:
    matchList = []
    if s.endswith('A') and s in samples_compar_dict.keys():
        s2 = samples_compar_dict[s]
        for m in Illumina_markers:
            if [UAS_Auto_STR_Data[s][m][0][1], UAS_Auto_STR_Data[s][m][1][1], UAS_Auto_STR_Data[s][m][0][4], UAS_Auto_STR_Data[s][m][1][4]] != [UAS_Auto_STR_Data[s2][m][0][1], UAS_Auto_STR_Data[s2][m][1][1], UAS_Auto_STR_Data[s2][m][0][4], UAS_Auto_STR_Data[s2][m][1][4]]:
                print ('no match:')
                print (s, m, [UAS_Auto_STR_Data[s][m][0][1], UAS_Auto_STR_Data[s][m][1][1], UAS_Auto_STR_Data[s][m][0][4], UAS_Auto_STR_Data[s][m][1][4]])
                print (s2, m, [UAS_Auto_STR_Data[s2][m][0][1], UAS_Auto_STR_Data[s2][m][1][1], UAS_Auto_STR_Data[s2][m][0][4], UAS_Auto_STR_Data[s2][m][1][4]])
            else:
                matchList.append(m)
        print(s, s2, 'match markers:', matchList)