# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'VectorPrep.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_VectorPrep(object):
    def setupUi(self, VectorPrep):
        VectorPrep.setObjectName("VectorPrep")
        VectorPrep.resize(388, 69)
        self.tagForw = QtWidgets.QLineEdit(VectorPrep)
        self.tagForw.setGeometry(QtCore.QRect(10, 10, 241, 20))
        self.tagForw.setObjectName("tagForw")
        self.tagRev = QtWidgets.QLineEdit(VectorPrep)
        self.tagRev.setGeometry(QtCore.QRect(10, 40, 241, 20))
        self.tagRev.setObjectName("tagRev")
        self.label = QtWidgets.QLabel(VectorPrep)
        self.label.setGeometry(QtCore.QRect(260, 10, 241, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(VectorPrep)
        self.label_2.setGeometry(QtCore.QRect(260, 40, 241, 16))
        self.label_2.setObjectName("label_2")

        self.retranslateUi(VectorPrep)
        QtCore.QMetaObject.connectSlotsByName(VectorPrep)

    def retranslateUi(self, VectorPrep):
        _translate = QtCore.QCoreApplication.translate
        VectorPrep.setWindowTitle(_translate("VectorPrep", "Vector Preparation"))
        self.tagForw.setText(_translate("VectorPrep", "TACG"))
        self.tagRev.setText(_translate("VectorPrep", "AAAC"))
        self.label.setText(_translate("VectorPrep", "Restriction site one"))
        self.label_2.setText(_translate("VectorPrep", "Restriction site two"))

