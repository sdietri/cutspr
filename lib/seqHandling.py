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

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os,tempfile,subprocess

from Bio.Blast import NCBIXML

class Primer:
    def __init__(self):
        self.tm = 0.0
        self.flankTm = 0.0
        self.seq = ""
        self.start = 0
        self.stop = 0
        self.strand = "+"
        self.blastRecord = None
        
def getRegionsForExtPrimer(flankA, flankB, flankSizeA, flankSizeB):
    #outsideprimers
    searchRegionA = ""
    if len(flankA) >= flankSizeA:
        searchRegionA = flankA[len(flankA)-flankSizeA:]
    else:
        #wraparound into the end of flankB
        searchRegionA = flankA
        searchRegionA = flankB[:flankSizeA-len(flankB)] + searchRegionA
        
    searchRegionB = ""
    if len(flankB) >= flankSizeB:
        searchRegionB = flankB[:flankSizeB]
    else:
        searchRegionB = flankB
        searchRegionB = searchRegionB + flankA[flankSizeB-len(flankA):]
    
    return (searchRegionA, searchRegionB)

def crawlForPrimer(seq, start, direction, targetTM):
    tm = 0.0
    startCrawl = start + 23
    extension = 0
    
    while tm < targetTM:
        primerSeq = seq[start:startCrawl+extension]
        tm = mt.Tm_NN(Seq(primerSeq))
        extension+=1
        if startCrawl+extension > len(seq):
            return
    
    primer = Primer()
    primer.tm = tm
    primer.start = start
    primer.stop = start + 23 + extension
    if direction == 1:
        primer.strand = "-"
        
    primer.seq = primerSeq

    return primer

def crawlForPrimerInt(flankA, flankB, insert=""):
    prims = []
    
    if insert=="" or len(insert) < 50:
        #flankA to B
        primer1Aportion = flankA[len(flankA)-12:]
        #print(primer1Aportion)
        tmAport = mt.Tm_NN(Seq(primer1Aportion))
        
        if tmAport > 31.0:
            while tmAport > 31.0:
                primer1Aportion = primer1Aportion[1:]
                tmAport = mt.Tm_NN(Seq(primer1Aportion))
        elif tmAport < 27.0:
            startCrawl = 12
            extension = 0
            while tmAport < 27.0:
                primer1Aportion = flankA[len(flankA)-(startCrawl+extension):]
                tmAport = mt.Tm_NN(Seq(primer1Aportion))
                extension+=1
                
        startCrawl = 23
        extension = 0
        tmBport = 0.0
        while tmBport < 55.0:
            primer1Bportion = flankB[0:startCrawl+extension]
            tmBport = mt.Tm_NN(Seq(primer1Bportion))
            extension+=1
        
        primer1 = Primer()
        primer1.tm = tmBport
        primer1.flankTm = tmAport
        primer1.strand = "+"
        primer1.seq = primer1Aportion.lower() + primer1Bportion
        
        #print("int1 " + str(primer1.tm) + " " + str(primer1.flankTm))
        
        #flank B to A
        flankArevc = str(Seq(flankA, IUPAC.unambiguous_dna).reverse_complement())
        flankBrevc = str(Seq(flankB, IUPAC.unambiguous_dna).reverse_complement())
        
        primer2Bportion = flankBrevc[len(flankBrevc)-12:]
        tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
        
        if tmBport2 > 31.0:
            while tmBport2 > 31.0:
                primer2Bportion = primer2Bportion[:len(primer2Bportion)-1]
                tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
        elif tmBport2 < 27.0:
            startCrawl = 12
            extension = 0
            while tmBport2 < 27.0:
                primer2Bportion = flankBrevc[len(flankBrevc)-(startCrawl+extension):]
                tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
                extension+=1
        
        startCrawl = 23
        extension = 0
        tmAport2 = 0.0

        
        while tmAport2 < 55.0:
            primer2Aportion = flankArevc[:startCrawl+extension]
            tmAport2 = mt.Tm_NN(Seq(primer2Aportion))
            extension+=1
        
        primer2 = Primer()
        primer2.tm = tmAport2 
        primer2.flankTm = tmBport2
        primer2.strand = "-"
        primer2.seq =primer2Bportion.lower() +  primer2Aportion
        
        #print("int2 " + str(primer2.tm) + " " + str(primer2.flankTm))
        
        #order of primers
        prims.append(primer2)
        prims.append(primer1)
        
    else: #with insert
        insert = insert.replace(os.linesep, "")
        #flankA to Insert
        primer1Aportion = flankA[len(flankA)-12:]
        tmAport = mt.Tm_NN(Seq(primer1Aportion))
        
        if tmAport > 31.0:
            while tmAport > 31.0:
                primer1Aportion = primer1Aportion[1:]
                tmAport = mt.Tm_NN(Seq(primer1Aportion))
        elif tmAport < 27.0:
            startCrawl = 12
            extension = 0
            while tmAport < 27.0:
                primer1Aportion = flankA[len(flankA)-(startCrawl+extension):]
                tmAport = mt.Tm_NN(Seq(primer1Aportion))
                extension+=1
                
        startCrawl = 23
        extension = 0
        tmBport = 0.0
        while tmBport < 55.0:
            primer1Bportion = insert[0:startCrawl+extension]
            tmBport = mt.Tm_NN(Seq(primer1Bportion))
            extension+=1
        
        primer1 = Primer()
        primer1.tm = tmBport
        primer1.flankTm = tmAport
        primer1.strand = "+"
        primer1.seq = primer1Aportion.lower().replace(" ","") + primer1Bportion.replace(" ","")
        
        #print("int1 " + str(primer1.tm) + " " + str(primer1.flankTm))
        
        #Insert to flankA
        flankArevc = str(Seq(flankA, IUPAC.unambiguous_dna).reverse_complement())
        insertrevC = str(Seq(insert, IUPAC.unambiguous_dna).reverse_complement())
        
        primer2Bportion = insertrevC[len(insertrevC)-12:]
        tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
        
        if tmBport2 > 31.0:
            while tmBport2 > 31.0:
                primer2Bportion = primer2Bportion[:len(primer2Bportion)-1]
                tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
        elif tmBport2 < 27.0:
            startCrawl = 12
            extension = 0
            while tmBport2 < 27.0:
                primer2Bportion = insertrevC[len(insertrevC)-(startCrawl+extension):]
                tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
                extension+=1
        
        startCrawl = 23
        extension = 0
        tmAport2 = 0.0

        
        while tmAport2 < 55.0:
            primer2Aportion = flankArevc[0:startCrawl+extension]
            tmAport2 = mt.Tm_NN(Seq(primer2Aportion))
            extension+=1
        
        primer2 = Primer()
        primer2.tm = tmAport2 
        primer2.flankTm = tmBport2
        primer2.strand = "-"
        primer2.seq = primer2Bportion.lower().replace(" ","") + primer2Aportion.replace(" ","")
        #print("int2 " + str(primer2.tm) + " " + str(primer2.flankTm))
        
        
        #Insert to FlankB

        primer1Aportion = insert[len(insert)-12:]
        tmAport = mt.Tm_NN(Seq(primer1Aportion))
        
        if tmAport > 31.0:
            while tmAport > 31.0:
                primer1Aportion = primer1Aportion[:len(primer1Aportion)-1]
                tmAport = mt.Tm_NN(Seq(primer1Aportion))
        elif tmAport < 27.0:
            startCrawl = 12
            extension = 0
            while tmAport < 27.0:
                primer1Aportion = insert[len(insert)-(startCrawl+extension):]
                if startCrawl+extension > len(insert):
                    #print("Insert too small")
                    return
                tmAport = mt.Tm_NN(Seq(primer1Aportion))
                extension+=1
                
        startCrawl = 23
        extension = 0
        tmBport = 0.0
        while tmBport < 55.0:
            primer1Bportion = flankB[0:startCrawl+extension]
            tmBport = mt.Tm_NN(Seq(primer1Bportion))
            extension+=1
        
        primer3 = Primer()
        primer3.tm = tmBport
        primer3.flankTm = tmAport
        primer3.strand = "+"
        primer3.seq = primer1Aportion.lower().replace(" ","") + primer1Bportion.replace(" ","")
        #print("int1 " + str(primer3.tm) + " " + str(primer3.flankTm))
        
        #FlankB to Insert
        insertrevC = str(Seq(insert, IUPAC.unambiguous_dna).reverse_complement())
        flankBrevc = str(Seq(flankB, IUPAC.unambiguous_dna).reverse_complement())
        
        primer2Bportion = flankBrevc[len(flankBrevc)-12:]
        tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
        
        if tmBport2 > 31.0:
            while tmBport2 > 31.0:
                primer2Bportion = primer2Bportion[:len(primer2Bportion)-1]
                tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
        elif tmBport2 < 27.0:
            startCrawl = 12
            extension = 0
            while tmBport2 < 27.0:
                primer2Bportion = flankBrevc[len(flankBrevc)-(startCrawl+extension):]
                tmBport2 = mt.Tm_NN(Seq(primer2Bportion))
                extension+=1
        
        startCrawl = 23
        extension = 0
        tmAport2 = 0.0

        
        while tmAport2 < 55.0:
            primer2Aportion = insertrevC[0:startCrawl+extension]
            tmAport2 = mt.Tm_NN(Seq(primer2Aportion))
            extension+=1
        
        primer4 = Primer()
        primer4.tm = tmAport2 
        primer4.flankTm = tmBport2
        primer4.strand = "-"
        primer4.seq = primer2Bportion.lower().replace(" ","") + primer2Aportion.replace(" ","")
        #print("int2 " + str(primer4.tm) + " " + str(primer4.flankTm))
        
        prims.append(primer2)
        prims.append(primer1)
        prims.append(primer4)
        prims.append(primer3)
        
    return prims
    
    
def evaluatePrimerBlastRes(blastRes):
    firstHit = True
    secondHit = None

    for blast_record in blastRes:
        queryletters = blast_record.query_letters

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hitPerc = float(hsp.identities) / queryletters * 100
                if firstHit:
                    if hitPerc == 100.0:
                        firstHit = False;
                    else:
                        print("bad start")
                        return
                else:
                    if hitPerc == 100.0:
                        print("bad primer")
                    elif secondHit == None:
                        secondHit = hsp

    if mt.Tm_NN(Seq(secondHit.query)) < 40.0:
        return True
    else:
        return False


def findInternalPrimers(flankA, flankB, insert=""):
    if len(insert) > 0:
        prims = crawlForPrimerInt(flankA, flankB, insert)
        return prims
    else:
        prims = crawlForPrimerInt(flankA, flankB)
        return prims


def findExternalPrimers(flankA, flankB, blastDB, blastnExe):
    primAcandidates = []
    
    for i in range(30,len(flankA)):
        if flankA[i] == 'G' or flankA[i] == 'C':
            primerA = crawlForPrimer(flankA, i, 0, 56)
            if primerA == None:
                return
            hitFasta = tempfile.mkstemp()[1]
            hitFastaHandle = open(hitFasta,"w")
            hitFastaHandle.write(">PrimerAext" + os.linesep)
            hitFastaHandle.write(primerA.seq + os.linesep)
            
            hitFastaHandle.flush()
            
            subprocess.call([blastnExe,"-query",hitFasta,"-db", blastDB, "-task","blastn-short", "-outfmt","5", "-out", hitFasta+".xml", "-num_threads","2"])
            #blastn = NcbiblastnCommandline(query=hitFasta, db=blastDB, task="blastn-short" ,outfmt=5,
            #                                 out=hitFasta+".xml", num_threads=2)
            #blastn.__call__()
            blastRecord = NCBIXML.parse(open(hitFasta+".xml"))
            
            if evaluatePrimerBlastRes(blastRecord):
                primerA.tm = mt.Tm_NN(Seq(primerA.seq))
                primAcandidates.append(primerA)
                break
            else:
                continue
            
            
                
        

    #print(str(primerA.tm) + " " + primerA.seq + " " + str(primerA.start) + " " + primerA.strand + " " + str(primerA.stop))
    
    revFlankB = str(Seq(flankB, IUPAC.unambiguous_dna).reverse_complement())
    primBcandidates = []

    for i in range(30,len(revFlankB)):
        if revFlankB[i] == 'G' or revFlankB[i] == 'C':

            primerB = crawlForPrimer(revFlankB, i, 1, 56)
            if primerB == None:
                return
            hitFasta = tempfile.mkstemp()[1]
            hitFastaHandle = open(hitFasta,"w")
            hitFastaHandle.write(">PrimerBext" + os.linesep)
            hitFastaHandle.write(primerB.seq + os.linesep)
            
            hitFastaHandle.flush()
            
            subprocess.call([blastnExe,"-query",hitFasta,"-db", blastDB, "-task","blastn-short", "-outfmt","5", "-out", hitFasta+".xml", "-num_threads","2"])
            #blastn = NcbiblastnCommandline(query=hitFasta, db=blastDB, task="blastn-short" ,outfmt=5,
            #                                 out=hitFasta+".xml", num_threads=2)
            #blastn.__call__()
            blastRecord = NCBIXML.parse(open(hitFasta+".xml"))
            
            if evaluatePrimerBlastRes(blastRecord):
                primerB.tm = mt.Tm_NN(Seq(primerB.seq))
                primBcandidates.append(primerB)
                break
            else:
                continue
            
    #print(str(primerB.tm) + " " + primerB.seq + " " + str(primerB.start) + " " + primerB.strand + " " + str(primerB.stop))
    
    if len(primAcandidates) == 0 or len(primBcandidates) == 0:
        return None
    
    return primAcandidates[0], primBcandidates[0]