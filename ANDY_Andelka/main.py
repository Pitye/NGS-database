from AutoSTR_FR_check_no_mismatches_ndi import autoSTR_FR_check_no_mismatches
from AutoSTR_FR_import2db_ndi import autoSTR_FR_import2db
from Y_STR_FR_check_no_mismatches_ndi import Y_STR_FR_check_no_mismatches
from Y_STR_FR_import2db_ndi import Y_STR_FR_import2db
from relationship_LR_estimation_ndi import relationship_LR_estimation
from print_samples_ndi import print_samples
from delete_samples_ndi import delete_samples
import os
import sys

from MySQLdb import connect

from functools import partial

from PyQt5.QtWidgets import QApplication

from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QGridLayout
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QStackedWidget

from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QCheckBox
from PyQt5.QtWidgets import QComboBox

from PyQt5.QtWidgets import QListWidget
from PyQt5.QtWidgets import QAbstractItemView

from PyQt5.QtWidgets import QFileDialog

from PyQt5.QtGui import QIntValidator

__version__ = '0.1'
__author__ = 'dh8'

class NdiUi(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('ANDY - NGS data interpreter')
        self.setFixedSize(700, 480)

        self.generalLayout = QVBoxLayout()
        self._centralWidget = QWidget(self)
        self.setCentralWidget(self._centralWidget)
        self._centralWidget.setLayout(self.generalLayout)

        self.buttonsModesLayout = QHBoxLayout()
        self.buttonMode1 = QPushButton('CHECK/IMPORT')
        self.buttonMode1.clicked.connect(self.setMode1)
        # self.buttonMode1 = QPushButton('Mode1', clicked=self.setMode1)
        self.buttonMode2 = QPushButton('RELATIONSHIP_ESTIMATION')
        self.buttonMode2.clicked.connect(self.setMode2)
        self.buttonMode3 = QPushButton('PRINT_SAMPLES')
        self.buttonMode3.clicked.connect(self.setMode3)
        self.buttonMode4 = QPushButton('DELETE_SAMPLES')
        self.buttonMode4.clicked.connect(self.setMode4)
        self.buttonsModesLayout.addWidget(self.buttonMode1)
        self.buttonsModesLayout.addWidget(self.buttonMode2)
        self.buttonsModesLayout.addWidget(self.buttonMode3)
        self.buttonsModesLayout.addWidget(self.buttonMode4)

        self.generalLayout.addLayout(self.buttonsModesLayout)

        self.modeWidget = QStackedWidget()
        self.generalLayout.addWidget(self.modeWidget)
        self.generateDbPasswWidget()

        self.mode1set = False
        self.mode2set = False
        self.mode3set = False
        self.mode4set = False

    def generateDbPasswWidget(self):
        self.dbPasswWidget = QWidget()

        self.dbPassw = QLineEdit(parent=self.dbPasswWidget)
        self.dbPassw.setEchoMode(QLineEdit.Password)
        self.dbPassw.setFixedSize(150, 20)
        self.dbPassw.setPlaceholderText("DB Password")

        self.dbTryConnection = QPushButton("Try DB connection", parent=self.dbPasswWidget)
        self.dbTryConnection.move(170, 0)
        self.dbTryConnection.clicked.connect(self.tryDbConnection)
        
        self.tryResult = QLabel(parent=self.dbPasswWidget)
        self.tryResult.setFixedSize(160, 20)
        self.tryResult.move(280, 2)

        self.modeWidget.addWidget(self.dbPasswWidget)
        self.modeWidget.setCurrentWidget(self.dbPasswWidget)
        
    def tryDbConnection(self):
        try:
            db = connect("localhost", "root", self.dbPassw.text(), "NGS_FORENSIC")
        except:
            self.tryResult.setText('Connection failed')
            return
        self.tryResult.setText('Connection OK')
        db.close()

    def getBrowseDirName(self, label):
        dir = QFileDialog.getExistingDirectory()
        label.setText(dir)
        
###todo - open file ad to mode 3
    def getBrowseFileName(self, label):
        file = QFileDialog.getOpenFileName()
        label.setText(file[0])
        
    def setDefaultModeButtonText(self):
        self.buttonMode1.setText('CHECK/IMPORT')
        self.buttonMode2.setText('RELATIONSHIP_ESTIMATION')
        self.buttonMode3.setText('PRINT_SAMPLES')
        self.buttonMode4.setText('DELETE_SAMPLES')

    def setMode1(self):
        if not self.mode1set:
            self.generateMode1Widget()
            self.modeWidget.addWidget(self.mode1Widget)
            self.mode1set = True
        self.modeWidget.setCurrentWidget(self.mode1Widget)
        self.setDefaultModeButtonText()
        self.buttonMode1.setText('<CHECK/IMPORT>')

    def generateMode1Widget(self):
        self.mode1Widget = QWidget()
        self.mode1Layout = QGridLayout()

        self.mode1buttonDir1 = QPushButton("in_directory Illumina")
        self.mode1labelDir1 = QLabel(os.path.normpath('C:/NGS_forensic_database/xlsx_detail_reports'))
        self.mode1buttonDir1.clicked.connect(partial(self.getBrowseDirName, self.mode1labelDir1))
        # self.mode1buttonDir1.clicked.connect(lambda: self.getBrowseDirName(self.mode1labelDir1))

        self.mode1buttonDir2 = QPushButton("out_directory")
        self.mode1labelDir2 = QLabel(os.path.normpath('C:/NGS_forensic_database/csv_output'))
        self.mode1buttonDir2.clicked.connect(partial(self.getBrowseDirName, self.mode1labelDir2))

        self.mode1buttonDir3 = QPushButton("xml_CE_directory")
        self.mode1labelDir3 = QLabel(os.path.normpath('C:/NGS_forensic_database/xml_CE'))
        self.mode1buttonDir3.clicked.connect(partial(self.getBrowseDirName, self.mode1labelDir3))

        self.mode1buttonDir4 = QPushButton("CSV_CE_directory")
        self.mode1labelDir4 = QLabel(os.path.normpath('C:/NGS_forensic_database/csv_CE'))
        self.mode1buttonDir4.clicked.connect(partial(self.getBrowseDirName, self.mode1labelDir4))

        self.mode1buttonDir5 = QPushButton("report_directory")
        self.mode1labelDir5 = QLabel()
        self.mode1buttonDir5.clicked.connect(partial(self.getBrowseDirName, self.mode1labelDir5))

        self.mode1labelNumber = QLabel("no_reads_for_validation:")
        self.mode1inputNumber = QLineEdit("150")
        self.onlyInt = QIntValidator()
        self.mode1inputNumber.setValidator(self.onlyInt)

        self.mode1booleanCheckBox = QCheckBox("CheckIfSampleInDatabase")
        self.mode1booleanCheckBox.setChecked(True)

        self.mode1scriptSelector = QComboBox()
        self.mode1scriptSelector.addItems(["AutoSTR_check_no_mismatches", "AutoSTR_import2db", "Y-STR_check_no_mismatches", "Y-STR_import2db"])

        self.mode1buttonRun = QPushButton("run")
        self.mode1buttonRun.clicked.connect(self.runMode1)
        
        self.mode1labelState = QLabel("")

        self.mode1Layout.addWidget(self.mode1buttonDir1, 0, 0)
        self.mode1Layout.addWidget(self.mode1labelDir1, 0, 1)
        self.mode1Layout.addWidget(self.mode1buttonDir2, 1, 0)
        self.mode1Layout.addWidget(self.mode1labelDir2, 1, 1)
        self.mode1Layout.addWidget(self.mode1buttonDir3, 2, 0)
        self.mode1Layout.addWidget(self.mode1labelDir3, 2, 1)
        self.mode1Layout.addWidget(self.mode1buttonDir4, 3, 0)
        self.mode1Layout.addWidget(self.mode1labelDir4, 3, 1)
        self.mode1Layout.addWidget(self.mode1buttonDir5, 4, 0)
        self.mode1Layout.addWidget(self.mode1labelDir5, 4, 1)
        self.mode1Layout.addWidget(self.mode1labelNumber, 5, 0)
        self.mode1Layout.addWidget(self.mode1inputNumber, 5, 1)
        self.mode1Layout.addWidget(self.mode1booleanCheckBox, 6, 0)
        self.mode1Layout.addWidget(self.mode1scriptSelector, 7, 0)
        self.mode1Layout.addWidget(self.mode1buttonRun, 7, 1)
        self.mode1Layout.addWidget(self.mode1labelState, 8, 1)
        self.mode1Widget.setLayout(self.mode1Layout)
        
    def runMode1(self):
        
        self.mode1labelState.setText("running")
        self.mode1labelState.repaint()
        
        dir1 = self.mode1labelDir1.text()
        dir2 = self.mode1labelDir2.text()
        dir3 = self.mode1labelDir3.text()
        dir4 = self.mode1labelDir4.text()
        dir5 = self.mode1labelDir5.text()
        number = int(self.mode1inputNumber.text())
        boolean = self.mode1booleanCheckBox.isChecked()
        script = self.mode1scriptSelector.currentIndex()
        
        #in_directory, out_directory, xml_directory, CSV_CE_directory, no_reads_for_validation, CheckIfSampleInDatabase
        if script == 0:
            autoSTR_FR_check_no_mismatches(dir1, dir2, dir3, dir4, dir5, number, boolean, self.dbPassw.text())
            self.mode1labelState.setText("finished")
        elif script == 1:
            autoSTR_FR_import2db(dir1, dir2, dir3, dir4, dir5, number, boolean, self.dbPassw.text())
            self.mode1labelState.setText("finished")
        elif script == 2:
            Y_STR_FR_check_no_mismatches(dir1, dir2, dir3, dir4, dir5, number, boolean, self.dbPassw.text())
            self.mode1labelState.setText("finished")
        elif script == 3:
            Y_STR_FR_import2db(dir1, dir2, dir3, dir4, dir5, number, boolean, self.dbPassw.text())
            self.mode1labelState.setText("finished")

    def setMode2(self):
        if not self.mode2set:
            self.generateMode2Widget()
            self.modeWidget.addWidget(self.mode2Widget)
            self.mode2set = True
        self.modeWidget.setCurrentWidget(self.mode2Widget)
        self.setDefaultModeButtonText()
        self.buttonMode2.setText('<RELATIONSHIP_ESTIMATION>')

    def generateMode2Widget(self):
        self.mode2Widget = QWidget()
        self.mode2Layout = QGridLayout()

        self.mode2sampleName1Label = QLabel("Sample name 1:")
        self.mode2sampleName1Input = QLineEdit()

        self.mode2sampleName2Label = QLabel("Sample name 2:")
        self.mode2sampleName2Input = QLineEdit()

        self.mode2boolean1CheckBox = QCheckBox("length_polymorphism_estimation")

        self.mode2boolean2CheckBox = QCheckBox("use_external_file_for_missing_markers")

        self.mode2buttonDir1 = QPushButton("CE_profiles_directory")
        self.mode2labelDir1 = QLabel("C:/NGS_forensic_database/relationship_LR_estimation/STR_lenght_profiles")
        self.mode2buttonDir1.clicked.connect(partial(self.getBrowseDirName, self.mode2labelDir1))

        self.mode2buttonDir2 = QPushButton("CE_frequencies_directory")
        self.mode2labelDir2 = QLabel("C:/NGS_forensic_database/relationship_LR_estimation/STR_lenght_frequencies")
        self.mode2buttonDir2.clicked.connect(partial(self.getBrowseDirName, self.mode2labelDir2))

        self.mode2buttonDir3 = QPushButton("reports_directory")
        self.mode2labelDir3 = QLabel("C:/NGS_forensic_database/relationship_LR_estimation/reports")
        self.mode2buttonDir3.clicked.connect(partial(self.getBrowseDirName, self.mode2labelDir3))

        self.mode2buttonRun = QPushButton("run")
        self.mode2buttonRun.clicked.connect(self.runMode2)
        
        self.mode2labelState = QLabel("")

        self.mode2Layout.addWidget(self.mode2sampleName1Label, 0, 0)
        self.mode2Layout.addWidget(self.mode2sampleName1Input, 0, 1)
        self.mode2Layout.addWidget(self.mode2sampleName2Label, 1, 0)
        self.mode2Layout.addWidget(self.mode2sampleName2Input, 1, 1)
        self.mode2Layout.addWidget(self.mode2boolean1CheckBox, 2, 0)
        self.mode2Layout.addWidget(self.mode2boolean2CheckBox, 3, 0)
        self.mode2Layout.addWidget(self.mode2buttonDir1, 4, 0)
        self.mode2Layout.addWidget(self.mode2labelDir1, 4, 1)
        self.mode2Layout.addWidget(self.mode2buttonDir2, 5, 0)
        self.mode2Layout.addWidget(self.mode2labelDir2, 5, 1)
        self.mode2Layout.addWidget(self.mode2buttonDir3, 6, 0)
        self.mode2Layout.addWidget(self.mode2labelDir3, 6, 1)
        self.mode2Layout.addWidget(self.mode2buttonRun, 7, 0)
        self.mode2Layout.addWidget(self.mode2labelState, 8, 1)

        self.mode2Widget.setLayout(self.mode2Layout)

    def addSamples(self):
        sample1 = self.mode2sampleName1Input.text()
        sample2 = self.mode2sampleName2Input.text()
        sampleList = list()
        sampleList.append(sample1)
        sampleList.append(sample2)
        return sampleList

    
    def runMode2 (self):
        self.mode2labelState.setText("running")
        self.mode2labelState.repaint()
        
        samples = self.addSamples()
        dir1 = self.mode2labelDir1.text()
        dir2 = self.mode2labelDir2.text()
        dir3 = self.mode2labelDir3.text()
        boolean1 = self.mode2boolean1CheckBox.isChecked()
        boolean2 = self.mode2boolean2CheckBox.isChecked()
        relationship_LR_estimation(samples, boolean1, boolean2, dir1, dir2, dir3, self.dbPassw.text())
        self.mode2labelState.setText("finished")
        
    def setMode3(self):
        if not self.mode3set:
            self.generateMode3Widget()
            self.modeWidget.addWidget(self.mode3Widget)
            self.mode3set = True
        self.modeWidget.setCurrentWidget(self.mode3Widget)
        self.setDefaultModeButtonText()
        self.buttonMode3.setText('<PRINT_SAMPLES>')

    def generateMode3Widget(self):
        self.mode3Widget = QWidget()
        self.mode3Layout = QGridLayout()

        self.mode3buttonFile1 = QPushButton("sample_list")
        self.mode3labelFile1 = QLabel("C:/NGS_forensic_database/print_samples/print_sample.txt")
        self.mode3buttonFile1.clicked.connect(partial(self.getBrowseFileName, self.mode3labelFile1))

        self.mode3buttonFile2 = QPushButton("locus_order")
        self.mode3labelFile2 = QLabel("C:/NGS_forensic_database/print_samples/locus_order.txt")
        self.mode3buttonFile2.clicked.connect(partial(self.getBrowseFileName, self.mode3labelFile2))
        
        self.mode3buttonDir3 = QPushButton("reports_directory")
        self.mode3labelDir3 = QLabel("C:/NGS_forensic_database/print_samples/prints")
        self.mode3buttonDir3.clicked.connect(partial(self.getBrowseDirName, self.mode3labelDir3))

        self.mode3categoriesList = QListWidget()
        self.mode3categoriesList.addItems(["AutoSTR", "Y-STR"])
        self.mode3categoriesList.setSelectionMode(QAbstractItemView.MultiSelection)

        #columns = ['allele', 'seq_name', 'PubMed_ID', 'sequence', 'no_reads', 'CE_validation', 'head_id', 'avg_no_reads', 'count_seq', 'frequency']
        self.mode3columnsList = QListWidget()
        self.mode3columnsList.addItems(["allele", "seq_name", "PubMed_ID", "sequence", "no_reads", "CE_validation", "head_id", "avg_no_reads", "count_seq", "frequency"])
        self.mode3columnsList.setSelectionMode(QAbstractItemView.MultiSelection)

        self.mode3buttonRun = QPushButton("run")
        self.mode3buttonRun.clicked.connect(self.runMode3)

        self.mode3labelState = QLabel("")

        self.mode3Layout.addWidget(self.mode3buttonFile1, 0, 0)
        self.mode3Layout.addWidget(self.mode3labelFile1, 0, 1)
        self.mode3Layout.addWidget(self.mode3buttonFile2, 1, 0)
        self.mode3Layout.addWidget(self.mode3labelFile2, 1, 1)
        self.mode3Layout.addWidget(self.mode3buttonDir3, 2, 0)
        self.mode3Layout.addWidget(self.mode3labelDir3, 2, 1)
        self.mode3Layout.addWidget(self.mode3categoriesList, 3, 0)
        self.mode3Layout.addWidget(self.mode3columnsList, 3, 1)
        self.mode3Layout.addWidget(self.mode3buttonRun, 4, 0)
        self.mode3Layout.addWidget(self.mode3labelState, 5, 1)

        self.mode3Widget.setLayout(self.mode3Layout)

    def SelectedItems(self,item):
        if item == "columns":
            selectedItems = self.mode3columnsList.selectedItems()
        if item == "category_mode3":
            selectedItems = self.mode3categoriesList.selectedItems()
        if item == "category_mode4":
            selectedItems = self.mode4categoriesList.selectedItems()
        selectedItemsTextList = list()
        for selectedItem in selectedItems:
            selectedItemsTextList.append(selectedItem.text())
        return selectedItemsTextList

    def runMode3 (self):
        self.mode3labelState.setText("running")
        self.mode3labelState.repaint()
        file1 = self.mode3labelFile1.text()
        file2 = self.mode3labelFile2.text()
        dir3 = self.mode3labelDir3.text()
        categoriesList = self.SelectedItems ("category_mode3")
        columnsList = self.SelectedItems ("columns")
        print_samples(file1, file2, dir3, categoriesList, columnsList, self.dbPassw.text())
        self.mode3labelState.setText("finished")

    def setMode4(self):
        if not self.mode4set:
            self.generateMode4Widget()
            self.modeWidget.addWidget(self.mode4Widget)
            self.mode4set = True
        self.modeWidget.setCurrentWidget(self.mode4Widget)
        self.setDefaultModeButtonText()
        self.buttonMode4.setText('<DELETE_SAMPLES>')

    def generateMode4Widget(self):
        self.mode4Widget = QWidget()
        self.mode4Layout = QGridLayout()

        self.mode4buttonFile1 = QPushButton("sample_list")
        self.mode4labelFile1 = QLabel("C:/NGS_forensic_database/delete_samples/delete_sample.txt")
        self.mode4buttonFile1.clicked.connect(partial(self.getBrowseFileName, self.mode4labelFile1))

        self.mode4categoriesList = QListWidget()
        self.mode4categoriesList.addItems(["AutoSTR", "Y-STR"])
        self.mode4categoriesList.setSelectionMode(QAbstractItemView.MultiSelection)

        self.mode4buttonRun = QPushButton("run")
        self.mode4buttonRun.clicked.connect(self.runMode4)
        
        self.mode4labelState = QLabel("")

        self.mode4Layout.addWidget(self.mode4buttonFile1, 0, 0)
        self.mode4Layout.addWidget(self.mode4labelFile1, 0, 1)
        self.mode4Layout.addWidget(self.mode4categoriesList, 1, 0)
        self.mode4Layout.addWidget(self.mode4buttonRun, 2, 0)
        self.mode4Layout.addWidget(self.mode4labelState, 1, 1)

        self.mode4Widget.setLayout(self.mode4Layout)

    def runMode4 (self):
        self.mode4labelState.setText("running")
        self.mode4labelState.repaint()
        file1 = self.mode4labelFile1.text()
        categoriesList = self.SelectedItems ("category_mode4")
        delete_samples(file1, categoriesList, self.dbPassw.text())
        self.mode4labelState.setText("finished")

def main():
    ndi = QApplication([])
    view = NdiUi()
    view.show()
    sys.exit(ndi.exec())

if __name__ == '__main__':
    main()
