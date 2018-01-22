#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:24:51 2017

Copyright 2017 Sascha Dietrich

 This file is part of CutSPR.

    CutSPR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CutSPR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CutSPR.  If not, see <http://www.gnu.org/licenses/>.

"""
from ui_files.main import Ui_CutSPR
from ui_files.VectorPrep import Ui_VectorPrep
from ui_files.BlastExeDialog import Ui_BlastExecutableSelection
from ui_files.authors import Ui_Authors
from PyQt5.QtWidgets import QMainWindow, QApplication, QStyleFactory, QFileDialog, QLineEdit, \
                            QLabel, QListWidgetItem, QHBoxLayout,QFrame, QDialog, QSpacerItem, QSizePolicy
from PyQt5.QtGui import QFont, QColor

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import MeltingTemp as mt

from Bio.Blast import NCBIXML

from lib.seqHandling import getRegionsForExtPrimer, findExternalPrimers, findInternalPrimers, Primer

import os,sys,re,subprocess,tempfile,shutil

regexTranslationTable = {}
regexTranslationTable["A"] = "A"
regexTranslationTable["T"] = "T"
regexTranslationTable["C"] = "C"
regexTranslationTable["G"] = "G"
regexTranslationTable["R"] = "[AG]"
regexTranslationTable["Y"] = "[CT]"
regexTranslationTable["S"] = "[CG]"
regexTranslationTable["W"] = "[AT]"
regexTranslationTable["K"] = "[GT]"
regexTranslationTable["M"] = "[AC]"
regexTranslationTable["B"] = "[CGT]"
regexTranslationTable["D"] = "[AGT]"
regexTranslationTable["H"] = "[ACT]"
regexTranslationTable["V"] = "[ACG]"
regexTranslationTable["N"] = "[ACGT]"

class PAMhit:
    def __init__(self, seq, start, strand, stop):
        self.start = start
        self.stop = stop
        
        self.ref = ""
        self.refStart = 0
        self.refStop = 0
        
        self.seq = seq
        self.strand = strand
        
        self.genRefCount = 0
        self.secondHitPerc = 0.0
        
class BlastExeDialog(QDialog, Ui_BlastExecutableSelection):
    def __init__(self):
        super(BlastExeDialog,self).__init__()
        self.setupUi(self)
        
    def writeIni(self, blastnExe, makeblastdbExe):
        f = open('cutspr.ini', "w")
        f.write("blastn=" + blastnExe + os.linesep)
        f.write("makeblastdb=" + makeblastdbExe + os.linesep)
        f.flush()
        f.close()
        
        self.close()
        
    def fileSelectDialog(self, target):
        files_types = "Executable (*)"
        exeFile = QFileDialog.getOpenFileName(self, 'Exe File', '' ,files_types)
        #print(exeFile[0])
        if not exeFile[0] == '':
            target.setText(exeFile[0])
    
class VectorPrepDialog(QDialog, Ui_VectorPrep):
    def __init__(self, overhangA, overhangB, vecType):
        super(VectorPrepDialog, self).__init__()
        
        self.setupUi(self)
        self.tagForw.setText(overhangA)
        self.tagRev.setText(overhangB)
        
class AuthorsDialog(QDialog, Ui_Authors):
    def __init__(self):
        super(AuthorsDialog, self).__init__()
        self.setupUi(self)
    
class MainDialog(QMainWindow, Ui_CutSPR):
    def __init__(self,scrH,scrW):
        super(MainDialog, self).__init__()
        
        self.setupUi(self)
        
        self.fileSelectionButton.clicked.connect(self.openFileSelectionDialog)
        self.cloningSitesgRNA.clicked.connect(lambda: self.openVectorPrep("sgRNA"))
        self.cloningSiteRecomb.clicked.connect(lambda: self.openVectorPrep("cassette"))
        self.resultsList.itemClicked.connect(self.sgRNAPrimerSelect)
        self.actionBLAST_executables.triggered.connect(lambda: self.openExeSelDialog())
        self.actionContact.triggered.connect(lambda: self.openAuthorsDialog())
        self.actionQuit.triggered.connect(self.close)
        self.cutRun.clicked.connect(self.runCut)
        self.flankASize.textChanged.connect(self.checkFlanks)
        self.flankBSize.textChanged.connect(self.checkFlanks)
        self.circularTargetCheck.toggled.connect(self.checkFlanks)
        self.sgRNAHeader.clicked.connect(lambda: self.toggleTextBlock(self.sgRNAdesc))
        self.indelCassHeader.clicked.connect(lambda: self.toggleTextBlock(self.indelCassDesc))
        self.vecConstHeader.clicked.connect(lambda: self.toggleTextBlock(self.vecConstDesc))
        
        self.sgRNAdesc.hide()
        self.indelCassDesc.hide()
        self.vecConstDesc.hide()
        
        self.foundTarget = False
        self.targetSelected = False
        self.genomeRefFasta = None
        
        
        self.insertSeq.textChanged.connect(self.insertSeqChanged)
        self.insert = ""
        
        self.sequences = {}
        
        self.overhangAsgRNA = "tacg"
        self.overhangBsgRNA = "aaac"
        self.overhangAcassette = "taggatccggccaacgaggcc"
        self.overhangBcassette = "taggatccggccttattggcc"
        
        self.blastnExe = None
        self.makeBlastDBExe = None
        
        self.gotBlastExes = False
        
        
        #hide/deactivate inactive elements until they're needed
        
        self.cassetteMsgBox.hide()
        self.tabWidget.setTabEnabled(2,False)
        self.tabWidget.setTabEnabled(3,False)
        
        #check for blast executables
        #Are we on windows, MacOS or linux?
        
        self.blastnExe = shutil.which("blastn")
        self.makeBlastDBExe = shutil.which("makeblastdb")
        
        if self.blastnExe == None or self.makeBlastDBExe == None:
            self.gotBlastExes = False
        else:
            self.gotBlastExes = True
            
        #do we have a cutspr.ini?
        #self.gotBlastExes = False
        
        if not self.gotBlastExes and os.path.isfile("cutspr.ini"):
            f = open("cutspr.ini")
            lines = f.readlines()
            
            gotBlastn = False
            gotMakeDB = False
            
            for line in lines:
                splitL = line.split("=")
                if len(splitL) == 2:
                    if splitL[0] == "blastn":
                        self.blastnExe = splitL[1].strip()
                        gotBlastn = True
                    if splitL[0] == "makeblastdb":
                        self.makeBlastDBExe = splitL[1].strip()
                        gotMakeDB = True
            
            if gotBlastn and gotMakeDB:
                self.gotBlastExes = True
        
        #Neither shutil nor cutspr.ini could provide the executables. They're neither installed nor in path. Ask the user to provide them and
        #point CutSPR directly to them. Save button writes the cutspr.ini which should then be used to find them.
        if not self.gotBlastExes:
            exes = self.openExeSelDialog("Could not find Blast Executables automatically")
            self.blastnExe = exes[0]
            self.makeBlastDBExe = exes[1]
            
    def toggleTextBlock(self, block):
        if block.isVisible():
            block.hide()
        else:
            block.show()
    
    def openAuthorsDialog(self):
        self.ad = AuthorsDialog()
        self.ad.show()
        self.ad.exec_()
        
    def openExeSelDialog(self,message=None):
        exeSelDialog = BlastExeDialog()
        if not message == None:
            exeSelDialog.openMessage.setText(message)
            
        if not self.blastnExe == None:
            exeSelDialog.blastnExeFile.setText(self.blastnExe)
        if not self.makeBlastDBExe == None:
            exeSelDialog.makeblastDBExeFile.setText(self.makeBlastDBExe)
            
        exeSelDialog.blastnExeFile.clicked.connect(lambda: exeSelDialog.fileSelectDialog(exeSelDialog.blastnExeFile))
        exeSelDialog.makeblastDBExeFile.clicked.connect(lambda: exeSelDialog.fileSelectDialog(exeSelDialog.makeblastDBExeFile))
        
        exeSelDialog.savButton.clicked.connect(lambda: exeSelDialog.writeIni(exeSelDialog.blastnExeFile.text(), exeSelDialog.makeblastDBExeFile.text()))
        
        exeSelDialog.show()
        exeSelDialog.raise_()
        exeSelDialog.activateWindow()
        exeSelDialog.exec_()
        
        return (exeSelDialog.blastnExeFile.text(),exeSelDialog.makeblastDBExeFile.text())
        
    def flipResultsList(self):
        if self.tabWidget.currentIndex() == 0:
            self.gsRNAsearch.layout().addWidget(self.resultsList)
        else:
            self.PrimerSearch.layout().addWidget(self.resultsList)
    
    def insertSeqChanged(self):
        #clean insertSeq of unwanted characters
        seqArray = self.insertSeq.toPlainText().strip().upper().split(os.linesep)
        seq = ""
        for curLine in seqArray:
            if curLine.startswith(">"):
                continue
            curLine = re.sub("[^ATCG]", "", curLine)
            seq = seq + curLine
        
        self.insert = seq
        #check for restriction sites
        if self.checkForRestrictionSites(self.insert):
            #print("Restriction Site Warning!")
            pass
        
        if len(self.insert) >= 50:
            self.insertStyleSwapper.setCurrentIndex(0)
            self.protocolSwapper.setCurrentIndex(1)
        else:
            self.insertStyleSwapper.setCurrentIndex(1)
            self.protocolSwapper.setCurrentIndex(0)
            
        self.updatePrimerList()
    
    def sgRNAPrimerSelect(self):
        self.sgRNA = str(Seq(self.resultsList.currentItem().seq, IUPAC.unambiguous_dna))
        self.updatePrimerList()
        self.cloningSitesgRNA.setEnabled(True)
        self.tabWidget.setTabEnabled(2,True)
        #check flank constraints, then enable tab 3 or alert in tab2 depending on that
        self.targetSelected = True
        if self.foundTarget:
            self.tabWidget.setTabEnabled(3,True)
        
        
    def updatePrimerList(self):
        #print("updatePrimerList")
        self.primerList.clear()
        prims = []
        
        #sgRNA primers
        sgRNAforw = Primer()
        sgRNAforw.seq = self.overhangAsgRNA + self.resultsList.currentItem().seq
        sgRNAforw.tm = mt.Tm_NN(Seq(sgRNAforw.seq))
        
        sgRNArev = Primer()
        sgRNArev.seq = self.overhangBsgRNA + str(Seq(self.resultsList.currentItem().seq, IUPAC.unambiguous_dna).reverse_complement())
        sgRNArev.tm = mt.Tm_NN(Seq(sgRNArev.seq))
        
        prims.append(sgRNAforw)
        prims.append(sgRNArev)
        
        #recombination cassette primers
        primAreg, primBreg = getRegionsForExtPrimer(self.flankA, self.flankB, self.targetFlankSizeA, self.targetFlankSizeB)
        if primAreg == None or primBreg == None:
            return
    
        try:
            primA, primB = findExternalPrimers(primAreg, primBreg,self.genomeRefFasta, self.blastnExe)
        except TypeError:
            #print("no suitable primers found")
            msg = "Could not find suitable primers at the requested flank target sizes!" + os.linesep + "Please check the input values above and increase the target flank sizes to allow for a bigger search space to find suitable primer candidates"
            self.cassetteMsgBox.setText(msg)
            self.cassetteMsgBox.show()
            self.tabWidget.setTabEnabled(3,False)
            self.tabWidget.tabBar().setTabTextColor(2, QColor(255,0,0))
            return

    
        #attach ligation overhangs to external primers
        primA.seq = self.overhangAcassette + primA.seq
        primB.seq = self.overhangBcassette + primB.seq
        
        #print(primA)
        #print(primB)
        prims.append(primA)
        for curPrim in findInternalPrimers(self.flankA, self.flankB, self.insert):
            prims.append(curPrim)
        
        prims.append(primB)
        
        count = 1
        for curPrim in prims:
            nTr = QListWidgetItem(self.primerList)
            hbox = QHBoxLayout()
            
            hbox.addWidget(QLabel(str(count)))
            seq = QLineEdit()
            #seq = QTextBrowser()
            seq.setText(curPrim.seq)
            seq.setReadOnly(True)
            seq.setFixedWidth(len(curPrim.seq)*10)
            hbox.addWidget(seq)
            
            spacer = QSpacerItem(10, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
            hbox.addItem(spacer)
            
            if curPrim.flankTm > 0.0:
                hbox.addWidget(QLabel(str(curPrim.flankTm)[:5]))
            else:
                hbox.addWidget(QLabel("   "))
            hbox.addWidget(QLabel(str(curPrim.tm)[:5]))
        
            count+=1
            
            frame = QFrame()
            frame.setLayout(hbox)
            nTr.setSizeHint(frame.minimumSizeHint())
            
            self.primerList.addItem(nTr)
            self.primerList.setItemWidget(nTr,frame)
            
    
    def checkForRestrictionSites(self, seq):
        if "GGTCTC" in seq:
            return True
        
        sfiIhit = re.search("GGCC[ATCG]{5}GGCC", seq)
        if sfiIhit:
            return True
    
        return False
            
        
    def enablePrimerSearchOld(self):
        self.tabWidget.setTabEnabled(1, True)
        self.vectorPrep.setEnabled(True)
        
        prims = []
        primAreg, primBreg = getRegionsForExtPrimer(self.flankA, self.flankB )
        primA, primB = findExternalPrimers(primAreg, primBreg,self.genomeRefFasta, self.blastnExe)
        #print(primA)
        #print(primB)
        prims.append(primA)
        for curPrim in findInternalPrimers(self.flankA, self.flankB, self.insertSeq.toPlainText()):
            prims.append(curPrim)
        
        prims.append(primB)
        
        self.primerList.clear()
        
        count = 1
        for curPrim in prims:
            nTr = QListWidgetItem(self.primerList)
            hbox = QHBoxLayout()
            
            hbox.addWidget(QLabel(str(count)))
            seq = QLineEdit()
            seq.setText(curPrim.seq)
            seq.setReadOnly(True)
            seq.setFixedWidth(len(curPrim.seq)*8)
            hbox.addWidget(seq)
            
            spacer = QSpacerItem(10, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
            hbox.addItem(spacer)
            
            hbox.addWidget(QLabel(str(curPrim.tm)[:5]))
            hbox.addWidget(QLabel(str(curPrim.flankTm)[:5]))

            
            count+=1
            
            frame = QFrame()
            frame.setLayout(hbox)
            nTr.setSizeHint(frame.minimumSizeHint())
            
            self.primerList.addItem(nTr)
            self.primerList.setItemWidget(nTr,frame)
        
        
    def openFileSelectionDialog(self):
        files = QFileDialog.getOpenFileNames(self,"Select Genome References","","Sequence Files (*.fas *.fa *.fna *.fasta *.embl *.gb *.gbk)")
        if len(files[0]) > 0:
            self.refGenFiles.clear()
            refFasta = tempfile.mkstemp()[1]
            o = open(refFasta, "w")
            for curFile in files[0]:
                self.refGenFiles.addItem(curFile)
                if curFile.endswith(".embl"):
                    for seq_record in SeqIO.parse(curFile, "embl"):
                        o.write(">" + seq_record.name + os.linesep)
                        o.write(str(seq_record.seq) + os.linesep)
                        self.sequences[seq_record.name] = str(seq_record.seq)
                
                elif curFile.endswith(".gbk") or curFile.endswith(".gb"):
                    for seq_record in SeqIO.parse(curFile, "genbank"):
                        o.write(">" + seq_record.name + os.linesep)
                        o.write(str(seq_record.seq) + os.linesep)
                        self.sequences[seq_record.name] = str(seq_record.seq)
                else:
                    f = open(curFile,"r")
                    for line in f.readlines():
                        o.write(line)
            o.flush()
            o.close()
            self.genomeRefFasta = refFasta
            
    def updateOverhangs(self, overhangA, overhangB, vecType):
        #print("change Overhangs")
        if vecType == "sgRNA":
            self.overhangAsgRNA = overhangA
            self.overhangBsgRNA = overhangB
        elif vecType == "cassette":
            self.overhangAcassette = overhangA
            self.overhangBcassette = overhangB
        
        self.updatePrimerList()
    
    def openVectorPrep(self, target):
        if hasattr(self, 'vp'):
            self.vp.close()

        if target == "sgRNA":
            self.vp = VectorPrepDialog(self.overhangAsgRNA, self.overhangBsgRNA, target)
        elif target == "cassette":
            self.vp = VectorPrepDialog(self.overhangAcassette, self.overhangBcassette, target)
        
        self.vp.tagForw.textChanged.connect(lambda: self.updateOverhangs(self.vp.tagForw.text(), self.vp.tagRev.text(),target))
        self.vp.tagRev.textChanged.connect(lambda: self.updateOverhangs(self.vp.tagForw.text(), self.vp.tagRev.text(),target))

        self.vp.show()
        self.vp.exec_()
            
    def checkFlanks(self):
        self.cassetteMsgBox.hide()
        
        if len(self.flankASize.text()) == 0 or len(self.flankBSize.text())== 0:
            return
        if int(self.flankASize.text()) == 0 or int(self.flankBSize.text()) == 0:
            return
        
        
        self.targetFlankSizeA = int(self.flankASize.text())
        self.targetFlankSizeB = int(self.flankBSize.text())
        
        if not self.circularTargetCheck.isChecked():
            self.cassetteMsgBox.setText("")
            nMsg = ""
            self.foundTarget = True
            if len(self.flankA) < self.targetFlankSizeA:  
                #print("flankA len " + str(len(self.flankA)))
                print("Flanks too small dialog!")
                nMsg = "FlankA is too short!" + os.linesep + "Maximum possible size is: " + str(len(self.flankA)) + os.linesep
                self.foundTarget = False
                
            nMsg = nMsg + os.linesep
            if len(self.flankB) < self.targetFlankSizeB:  
                #print("flankB len " + str(len(self.flankB)))
                print("Flanks too small dialog!")
                nMsg = nMsg + "FlankB is too short!" + os.linesep + "Maximum possible size is: " + str(len(self.flankB)) + os.linesep
                self.foundTarget = False

            if not self.foundTarget:
                print("Target not found Dialog!")
                self.cassetteMsgBox.setText(nMsg)
                self.cassetteMsgBox.show()
                self.tabWidget.setTabEnabled(3,False)
                
                self.tabWidget.tabBar().setTabTextColor(2, QColor(255,0,0))
                
                return
            elif self.targetSelected:
                self.tabWidget.setTabEnabled(3, True)
                self.tabWidget.tabBar().setTabTextColor(2, QColor(0,0,0))
        else:
            self.tabWidget.setTabEnabled(3, True)
            self.tabWidget.tabBar().setTabTextColor(2, QColor(0,0,0))
        self.updatePrimerList()
        
    def runCut(self):
        #self.tabWidget.setTabEnabled(1, False)
        self.resultsList.clear()
        self.cloningSitesgRNA.setEnabled(False)
        self.cassetteMsgBox.hide()
        
        self.targetFlankSizeA = int(self.flankASize.text())
        self.targetFlankSizeB = int(self.flankBSize.text())
        self.circularTarget = self.circularTargetCheck.isChecked()
        
        regex = self.getPAMregex()
        delSeq = ""
        #check if more than one > in delSeq
        delSeqRaw = self.deletionTargetSeq.toPlainText().split("\n")
        for line in delSeqRaw:
            if line.startswith(">"):
                continue
            
            delSeq = delSeq + line.replace("\n", "")
        
        delSeq = delSeq.upper()
        pamhits = {}
        
        matches = re.finditer(regex, delSeq)
        hitCount = 0
        for curhit in matches:
            hitCount+=1
            #check for bsaI and sfiI hits
            if not self.checkForRestrictionSites(curhit.group(0)):
                pamhits["PAMhit_" + str(hitCount)] = PAMhit(curhit.group(0), curhit.start(0), "+" ,curhit.end(0))
        
        delSeqRevComp = str(Seq(delSeq, IUPAC.unambiguous_dna).reverse_complement())
        matches = re.finditer(regex, delSeqRevComp)
        for curhit in matches:
            hitCount+=1
            pamhits["PAMhit_" + str(hitCount)] = PAMhit(curhit.group(0), curhit.start(0), "-" ,curhit.end(0))
            
        hitFasta = tempfile.mkstemp()[1]
        hitFastaHandle = open(hitFasta,"w")
        
        for curhit in pamhits:
            hitFastaHandle.write(">"+curhit + os.linesep)
            hitFastaHandle.write(pamhits[curhit].seq + os.linesep)
        hitFastaHandle.flush()
        hitFastaHandle.close()

        #prep blast db
        subprocess.call([self.makeBlastDBExe,"-in", self.genomeRefFasta, "-dbtype","nucl"])
        
        #blast deletion target to get Start Stop in genome for Flank cutting
        delFasta = tempfile.mkstemp()[1]
        delFastaHandle = open(delFasta, "w")
        delFastaHandle.write(">deletionTarget\n")
        delFastaHandle.write(delSeq + "\n")
        delFastaHandle.flush()
        delFastaHandle.close()
        
        #print(self.genomeRefFasta)
        subprocess.call([self.blastnExe,"-query",delFasta,"-db", self.genomeRefFasta, "-task","blastn-short", "-outfmt","5", "-out", delFasta+".xml", "-num_threads","2"])
        
        self.foundTarget = False
        
        blast_records = NCBIXML.parse(open(delFasta+".xml"))
        for blast_record in blast_records:
            queryletters = blast_record.query_letters
            for alignment in blast_record.alignments:
                ref = alignment.hit_def
                for hsp in alignment.hsps:
                    hitPerc = float(hsp.identities) / queryletters * 100
                    if hitPerc == 100.0:
                        self.foundTarget = True
                        self.flankA = self.sequences[ref][:hsp.sbjct_start-1]
                        self.flankB = self.sequences[ref][hsp.sbjct_end:]
                        
                        if not self.circularTarget:
                            self.cassetteMsgBox.setText("")
                            nMsg = ""
                            if len(self.flankA) < self.targetFlankSizeA:  
                                #print("flankA len " + str(len(self.flankA)))
                                #print("Flanks too small dialog!")
                                nMsg = "FlankA is too short!" + os.linesep + "Maximum possible size is: " + str(len(self.flankA)) + os.linesep
                                self.foundTarget = False
                                
                            nMsg = nMsg + os.linesep
                            if len(self.flankB) < self.targetFlankSizeB:  
                                #print("flankB len " + str(len(self.flankB)))
                                #print("Flanks too small dialog!")
                                nMsg = nMsg + "FlankB is too short!" + os.linesep + "Maximum possible size is: " + str(len(self.flankB)) + os.linesep
                                self.foundTarget = False
        
        if not self.foundTarget:
            print("Target not found Dialog!")
            self.cassetteMsgBox.setText(nMsg)
            self.cassetteMsgBox.show()
            self.tabWidget.setTabEnabled(3,False)
            
            self.tabWidget.tabBar().setTabTextColor(2, QColor(255,0,0))
            
            return
        elif self.targetSelected:
            self.tabWidget.setTabEnabled(3, True)
            self.tabWidget.tabBar().setTabTextColor(2, QColor(0,0,0))
        #blast our candidates
        subprocess.call([self.blastnExe,"-query",hitFasta,"-db", self.genomeRefFasta, "-task" ,"blastn-short", "-outfmt","5", "-out", hitFasta+".xml", "-num_threads","2"])
        
        blast_records = NCBIXML.parse(open(hitFasta+".xml"))
        for blast_record in blast_records:
            queryletters = blast_record.query_letters

            for alignment in blast_record.alignments:
                ref = alignment.hit_def
                for hsp in alignment.hsps:
                    hitPerc = float(hsp.identities) / queryletters * 100
                    if hitPerc == 100.0:
                        pamhits[blast_record.query].refStart = hsp.sbjct_start
                        pamhits[blast_record.query].refStop = hsp.sbjct_end
                        pamhits[blast_record.query].ref = ref
                        pamhits[blast_record.query].genRefCount+=1
                    elif blast_record.query in pamhits and pamhits[blast_record.query].secondHitPerc == 0.0:
                        pamhits[blast_record.query].secondHitPerc = hitPerc
        
        #Sort hits by their genRefCount
        orderedHits = []
        for curhit in pamhits:
            if len(orderedHits) == 0:
                orderedHits.append(curhit)
                continue
            
            targetPos = 0
            for curOrdHit in orderedHits:
                if pamhits[curhit].secondHitPerc > pamhits[curOrdHit].secondHitPerc:
                    targetPos+=1

            orderedHits.insert(targetPos, curhit)
            
        #Display results
        for curhit in orderedHits:
            if not pamhits[curhit].genRefCount == 1:
                continue #only unique sequences wanted
            nTr = QListWidgetItem(self.resultsList)
            hbox = QHBoxLayout()
            
            
            targetLength = 20
            if len(self.targetLength.text()) > 0:
                targetLength = int(self.targetLength.text())
            
            
            hitseq = QLineEdit()
            hitseq.setText(pamhits[curhit].seq[:targetLength])
            hitseq.setReadOnly(True)
            font = QFont("Monospace",9,QFont.Bold,False)   
            hitseq.setFont(font)
            hitseq.setFixedWidth(len(pamhits[curhit].seq)*8)
            hbox.addWidget(hitseq)
            
            secondHitPerc = "%.2f" % pamhits[curhit].secondHitPerc
            
            secondHitPercLabel = QLabel(secondHitPerc)
            secondHitPercLabel.setStyleSheet("QLabel { background-color: rgba(0,250,0,250); }")
            if pamhits[curhit].secondHitPerc > 50.0:
                secondHitPercLabel.setStyleSheet("QLabel { background-color: rgba(250,250,0,250); }")
            if pamhits[curhit].secondHitPerc > 70.0:
                secondHitPercLabel.setStyleSheet("QLabel { background-color: rgba(250,120,0,250); }")
            if pamhits[curhit].secondHitPerc > 90.0:
                secondHitPercLabel.setStyleSheet("QLabel { background-color: rgba(250,0,0,250); }")
            hbox.addWidget(secondHitPercLabel)
            
            nTr.__setattr__("hit", pamhits[curhit])
            nTr.__setattr__("seq", pamhits[curhit].seq[:targetLength])
            
            frame = QFrame()
            frame.setLayout(hbox)
            nTr.setSizeHint(frame.minimumSizeHint())
            
            self.resultsList.addItem(nTr)
            self.resultsList.setItemWidget(nTr,frame)
            
        
        if self.resultsList.count() == 0:
            self.resultsList.addItem("No unique hits in target")
        
    def getPAMregex(self):
        pat = self.PAMpattern.text()
        regexPat = ""
        for letter in pat.upper():
            regexPat= regexPat + regexTranslationTable[letter]
        
        targetLength = "20"
        if len(self.targetLength.text()) > 0:
            targetLength = self.targetLength.text()
            
        return "([ACTG]{" + targetLength + "})" + regexPat
        
app = QApplication(sys.argv)
screen_resolution = app.desktop().screenGeometry()
width, height = screen_resolution.width(), screen_resolution.height()

app.setStyle(QStyleFactory.create("fusion"))

mw = MainDialog(height,width)
mw.show()

#frameGeo = mw.frame.geometry()
#print(frameGeo.width())

sys.exit(app.exec_())