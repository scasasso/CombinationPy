##############################################
## Class to retrieve the cross section,
## branching ratio and width for the H
## of THDM (any type)
##############################################

import os, sys
from ROOT import *
import ROOT
import math


class THDM_CSandWidth:

    def __init__(self,_initDir,_systDir):

        # Constructor
        self.initDir = _initDir
        self.systDir = _systDir
        self.subDirNameBase = "mH{0}_tanB{1}_cosBMA{2}_2HDMType{3}_{4}TeV"
        self.summFileNameBase = "SusHi_mH{0}_tanB{1}_cosBMA{2}_2HDMType{3}_{4}TeV_summary.out"
        self.loadIncludes()
        self.theSMCSW = HiggsCSandWidth("YR3",ROOT.gSystem.ExpandPathName("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/txtFiles/"))


    # Load libraries needed to compute SM XS and BRs
    def loadIncludes(self):

        ROOT.gSystem.AddIncludePath(os.path.expandvars("$CMSSW_BASE/src/"))
        ROOT.gSystem.Load("libHiggsHiggs_CS_and_Width.so")
        ROOT.gROOT.LoadMacro(os.path.expandvars("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h+"))


    # Basic trigonometry
    # Convention -1 < cos(b-a) < 1, sin(b-a) > 0
    # (tanB > 0, 0 < B < pi/2)    
    def getBeta(self,tanB):

        # python math should return always angles in [-pi/2,+pi/2]
        beta = math.atan(tanB)
        if beta < 0:
            print "[THDM_CSandWidth::getBeta]: Something weird happening, math.atan is supposed to return always angles in  [-pi/2,+pi/2]"
            exit(1)
        return beta


    def getAlpha(self,tanB,cba):

        sba = math.sqrt(1-cba*cba)
        beta = self.getBeta(tanB)
        alpha = beta - math.asin(sba)

        if not ((beta-alpha)<math.pi and (beta-alpha)>0 and beta<math.pi/2 and beta>0 and alpha>-math.pi/2 and alpha<math.pi/2):
            print "[THDM_CSandWidth::getAlpha]: trigonometry failed! beta-alpha = ",(beta-alpha),", beta = ",beta,", alpha = ",alpha
            exit(1)
            
        return alpha


    # Function to retrieve the summary file
    def getSummFile(self,mH,tanB,cba,type,sqrts):

        subDirName = self.subDirNameBase.format(str(float(mH)),str(float(tanB)),str(float(cba)),str(int(type)),str(int(sqrts)))
        fileName = self.summFileNameBase.format(str(float(mH)),str(float(tanB)),str(float(cba)),str(int(type)),str(int(sqrts)))
        fileNameFull = self.initDir + "/" + subDirName + "/" + fileName

        if not os.path.exists(fileNameFull):
            print "[THDM_CSandWidth::xs_ggH]: Input file not found!"
            print fileNameFull
            exit(1)
            
        return fileNameFull
        

    # H width
    def H_width(self,mH,tanB,cba,type,sqrts):

        xsFileName = self.getSummFile(mH,tanB,cba,type,sqrts)
        xsFile = open(xsFileName,"r")
        Hwidth = float(xsFile.readlines()[0].split()[0].strip())        
        # print ggHxs

        return Hwidth


    # ggH cross section
    def xs_ggH(self,mH,tanB,cba,type,sqrts):

        xsFileName = self.getSummFile(mH,tanB,cba,type,sqrts)
        xsFile = open(xsFileName,"r")
        ggHxs = float(xsFile.readlines()[0].split()[1].strip())        
        # print ggHxs

        return ggHxs


    # bbH cross section
    def xs_bbH(self,mH,tanB,cba,type,sqrts):

        xsFileName = self.getSummFile(mH,tanB,cba,type,sqrts)
        xsFile = open(xsFileName,"r")
        bbHxs = float(xsFile.readlines()[0].split()[2].strip())        
        # print bbHxs

        return bbHxs


    # VBF cross section
    def xs_VBF(self,mH,tanB,cba,type,sqrts):

        xsSM = self.theSMCSW.HiggsCS(2,mH,sqrts)
        VBFxs = xsSM*cba*cba

        return VBFxs


    # WH cross section
    def xs_WH(self,mH,tanB,cba,type,sqrts):

        xsSM = self.theSMCSW.HiggsCS(3,mH,sqrts)
        WHxs = xsSM*cba*cba

        return WHxs


    # ZH cross section
    def xs_ZH(self,mH,tanB,cba,type,sqrts):

        xsSM = self.theSMCSW.HiggsCS(4,mH,sqrts)
        ZHxs = xsSM*cba*cba

        return ZHxs


    # ttH cross section
    def xs_ttH(self,mH,tanB,cba,type,sqrts):

        xsSM = self.theSMCSW.HiggsCS(5,mH,sqrts)
        beta = self.getBeta(tanB)
        alpha = self.getAlpha(tanB,cba)
        ttHxs = xsSM*(math.sin(alpha)/math.sin(beta))*(math.sin(alpha)/math.sin(beta))

        return ttHxs
            


    # Branching ratio H->ZZ
    def br_HZZ(self,mH,tanB,cba,type,sqrts):

        xsFileName = self.getSummFile(mH,tanB,cba,type,sqrts)
        xsFile = open(xsFileName,"r")
        brHZZ = float(xsFile.readlines()[0].split()[3].strip())        
        # print brHZZ

        return brHZZ


    # Branching ratio H->ZZ->4l (no interference!)
    def br_H4l(self,mH,tanB,cba,type,sqrts):

        brHZZ = self.br_HZZ(mH,tanB,cba,type,sqrts)
        brZllSquared = self.theSMCSW.HiggsBR(13,mH)/self.theSMCSW.HiggsBR(11,mH)
        
        return brHZZ*brZllSquared
    
        



       

        



    


