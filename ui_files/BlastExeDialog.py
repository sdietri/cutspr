# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'BlastExeDialog.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_BlastExecutableSelection(object):
    def setupUi(self, BlastExecutableSelection):
        BlastExecutableSelection.setObjectName("BlastExecutableSelection")
        BlastExecutableSelection.resize(849, 351)
        BlastExecutableSelection.setSizeGripEnabled(False)
        self.verticalLayoutWidget = QtWidgets.QWidget(BlastExecutableSelection)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 848, 351))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetMaximumSize)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.openMessage = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.openMessage.setText("")
        self.openMessage.setObjectName("openMessage")
        self.verticalLayout.addWidget(self.openMessage)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.label_3 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.label_4 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_4.setOpenExternalLinks(True)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.label_5 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.verticalLayout.addWidget(self.label_5)
        self.label_6 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_6.setObjectName("label_6")
        self.verticalLayout.addWidget(self.label_6)
        spacerItem1 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.label = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.blastnExeFile = CQLineEdit(self.verticalLayoutWidget)
        self.blastnExeFile.setObjectName("blastnExeFile")
        self.verticalLayout.addWidget(self.blastnExeFile)
        self.label_2 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.makeblastDBExeFile = CQLineEdit(self.verticalLayoutWidget)
        self.makeblastDBExeFile.setObjectName("makeblastDBExeFile")
        self.verticalLayout.addWidget(self.makeblastDBExeFile)
        self.savButton = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.savButton.setMaximumSize(QtCore.QSize(80, 16777215))
        self.savButton.setObjectName("savButton")
        self.verticalLayout.addWidget(self.savButton)

        self.retranslateUi(BlastExecutableSelection)
        QtCore.QMetaObject.connectSlotsByName(BlastExecutableSelection)

    def retranslateUi(self, BlastExecutableSelection):
        _translate = QtCore.QCoreApplication.translate
        BlastExecutableSelection.setWindowTitle(_translate("BlastExecutableSelection", "Select Blast Executables"))
        self.label_3.setText(_translate("BlastExecutableSelection", "To obtain the latest NCBI-Blast Executables, go to:"))
        self.label_4.setText(_translate("BlastExecutableSelection", "<a href=\"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download\">NCBI Blast Download</a>"))
        self.label_5.setText(_translate("BlastExecutableSelection", "and install the version for your OS and restart CutSPR. In case this dialog shows up again, please select the corresponding executables"))
        self.label_6.setText(_translate("BlastExecutableSelection", "with the dialog fields below:"))
        self.label.setText(_translate("BlastExecutableSelection", "Blastn Executable:"))
        self.label_2.setText(_translate("BlastExecutableSelection", "MakeBlastDB Executable:"))
        self.savButton.setText(_translate("BlastExecutableSelection", "Save"))

from .cqlineedit import CQLineEdit