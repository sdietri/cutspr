# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import QLabel
from PyQt5.QtCore import Qt, pyqtSignal

class CQlabel(QLabel):
    clicked= pyqtSignal()
    def __init__(self, title):
        super().__init__(title)
        #self.setAcceptDrops(True)
    
    def mouseReleaseEvent(self, e):
        self.clicked.emit()
        