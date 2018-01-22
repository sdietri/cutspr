# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import QLineEdit
from PyQt5.Qt import pyqtSignal
from PyQt5.QtCore import Qt

class CQLineEdit(QLineEdit):
    clicked = pyqtSignal()
    
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton: 
            self.clicked.emit()
        else:
            super().mousePressEvent(event)