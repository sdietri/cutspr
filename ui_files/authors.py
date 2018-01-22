# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Authors.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Authors(object):
    def setupUi(self, Authors):
        Authors.setObjectName("Authors")
        Authors.resize(301, 162)
        self.label = QtWidgets.QLabel(Authors)
        self.label.setGeometry(QtCore.QRect(10, 10, 59, 14))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Authors)
        self.label_2.setGeometry(QtCore.QRect(10, 30, 281, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(Authors)
        self.label_3.setGeometry(QtCore.QRect(20, 50, 221, 51))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(Authors)
        self.label_4.setGeometry(QtCore.QRect(180, 50, 59, 14))
        self.label_4.setOpenExternalLinks(True)
        self.label_4.setObjectName("label_4")
        self.label_8 = QtWidgets.QLabel(Authors)
        self.label_8.setGeometry(QtCore.QRect(10, 110, 281, 16))
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(Authors)
        self.label_9.setGeometry(QtCore.QRect(20, 130, 131, 16))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(Authors)
        self.label_10.setGeometry(QtCore.QRect(180, 130, 59, 14))
        self.label_10.setOpenExternalLinks(True)
        self.label_10.setObjectName("label_10")

        self.retranslateUi(Authors)
        QtCore.QMetaObject.connectSlotsByName(Authors)

    def retranslateUi(self, Authors):
        _translate = QtCore.QCoreApplication.translate
        Authors.setWindowTitle(_translate("Authors", "Dialog"))
        self.label.setText(_translate("Authors", "Authors:"))
        self.label_2.setText(_translate("Authors", "Corresponding Author, Method development:"))
        self.label_3.setText(_translate("Authors", "<html><head/><body><p>Dr. Robert Hertel, </p><p>Tobias Schilling</p></body></html>"))
        self.label_4.setText(_translate("Authors", "<a href=\"mailto:rhertel@gwdg.de?Subject=CutSPR%20Question\">Send Mail</a>"))
        self.label_8.setText(_translate("Authors", "CutSPR development:"))
        self.label_9.setText(_translate("Authors", "<html><head/><body><p>Dr. Sascha Dietrich</p></body></html>"))
        self.label_10.setText(_translate("Authors", "<a href=\"mailto:sdietri@gwdg.de?Subject=CutSPR%20Question\">Send Mail</a>"))

