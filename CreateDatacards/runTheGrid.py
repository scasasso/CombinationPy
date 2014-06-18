import optparse, re
import sys, os, pwd, commands
import subprocess
import shlex
import shutil, errno
import distutils.core
import glob


dcDir = "/afs/cern.ch/work/s/scasasso/private/HZZ4L_Combination_Git_forHighMass/HZZ4L_CombinationPy/CreateDatacards"
hcgDir = dcDir+"/CMSSW_6_1_1/src"
#massFile = dcDir+"/masses_HMP_above400.txt"
massFile = dcDir+"/masses_HMP_below400.txt"
#massFile = dcDir+"/masses_test.txt"
#massFile = dcDir+"/masses_HMP.txt"
scriptsDir = "/afs/cern.ch/work/s/scasasso/private/HZZ4L_Combination_Git_forHighMass/HZZ4L_CombinationPy/RunLimits/scripts"
theQueue = "1nd"

theIDBase = "cpsq{0}_brnew{1}"

# cprimes = [1.0, 0.1, 0.2, 0.3, 0.5, 0.7]
# brnews = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

cprimes = [1.0,0.1,0.2,0.3,0.5,0.7]
brnews = [0.0,0.1,0.2,0.3,0.4,0.5]

def parseOption():

    usage = ('usage: %prog runTheGrid.py [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    parser.add_option('-a', '--append', dest='append', type='string', default="",    help='String to append to the name of the datacards directory')        
    parser.add_option('-c', action='store_true', dest='create', default=False,
                      help='Create datacards: if False the program looks for cards in the directory cards_*TeV_cpsq10_brnew00_APPEND, where APPEND is the string specified with option -a')    

    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if not opt.append == "": opt.append = "_"+opt.append



def processCmd(cmd):

    args = shlex.split(cmd)
    sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()

    print err
    
    return out, err



def createDCs(aCPrime, aBRnew):


    os.chdir(dcDir)
    
    print "Creating datacards..."

    theID = theIDBase.format(str(aCPrime).replace(".",""),str(aBRnew).replace(".",""))
    # print theID

    os.chdir(dcDir)
    # Preparing the datacards
    print "   Creating 8 TeV datacards ..."
    dcCmd = "python makeDCsandWSs.py -i SM_inputs_{0}_HMP -a {0}_"+str(theID)+opt.append+" -b -j1 -c"+str(aCPrime)+" -r"+str(aBRnew)
    if aCPrime < 1.: dcCmd += " -w"
    out, err = processCmd(dcCmd.format("8TeV"))
    print "   Creating 7 TeV datacards ..."
    out, err = processCmd(dcCmd.format("7TeV"))

    oStream = open('/tmp/scasasso/stdout_'+theID+'.txt','w')
    oStream.write(out)
    eStream = open('/tmp/scasasso/stderr_'+theID+'.txt','w')
    eStream.write(err)

    oStream.close()
    eStream.close()

    os.chdir(dcDir)



# def createSVNstruct(inDir=".",outDir="/tmp/scasasso"):        
    

# def createHCGstruct(inDir="."):


def copyDCsToHCG(inDir, inSubDir, outDir):

    os.chdir(inDir)

    print "   Copying datacards to HCG directory..."
    if os.path.exists(outDir):
        print "Directory already exist: ",outDir
        exit
    os.makedirs(outDir)

    match8TeV = inSubDir.format("8","{0}")
    match7TeV = inSubDir.format("7","{0}")
    massFileBuf = open(massFile,"r")            
    for aMass in massFileBuf:
        theMassStr =  str(aMass.split()[0])
        massOutDir = outDir+"/"+theMassStr
        os.makedirs(massOutDir)
        for aFile in glob.glob(r""+match8TeV.format(theMassStr)): shutil.copy2(aFile, massOutDir)
        for aFile in glob.glob(r""+match7TeV.format(theMassStr)): shutil.copy2(aFile, massOutDir)

    os.chdir(dcDir)


def copyScripts(dest="."):
    
    shutil.copy(massFile,dest)

    for anExt in ["h","cxx","py","cc","sh","C"]:
        for aFile in glob.glob(r""+scriptsDir+"/*."+anExt): shutil.copy2(aFile,dest)
    shutil.copy2(scriptsDir+"/hadd2",dest)
    shutil.copy2(scriptsDir+"/massrange",dest)            
    distutils.dir_util.copy_tree(scriptsDir+"/cms", dest+"/cms")
    distutils.dir_util.copy_tree(scriptsDir+"/plotting",dest+"/plotting/")
    distutils.dir_util.copy_tree(scriptsDir+"/grids",dest+"/grids/")            
    
    

def combineDCsAndWSs():

    
    copyScripts()
    #Combining the datacards and the workspaces
    print "   Combining datacards and workspaces ..."
    outComb, errComb = processCmd("./makeCardsWorkspaces.sh "+massFile+" -2")
    print outComb
    print errComb    


    

def runTheLimits(theDir,aCPrime,aBRnew):

    os.chdir(theDir)

    #Running the asymptotic limits
    print "   Running the limits on lxbatch ..."    
    options = str(aCPrime)+' '+str(aBRnew)
    options += ' -E'
    # print 'DEBUG: string = python make_parallel_limits_BSM.py -M ASCLS -f '+massFile+' -q lsf -t combine -o "'+options+'"'
    outSub, errSub = processCmd('python make_parallel_limits_BSM.py -M ASCLS -f '+massFile+' -q lsf -n '+theQueue+' -t combine -o "'+options+'"')
    print outSub
    print errSub

    os.chdir(dcDir)


    

if __name__ == "__main__":

    parseOption()

    refDir = hcgDir+"/HCG_"+theIDBase.format(str(1.0).replace(".",""),str(0.0).replace(".",""))+opt.append
    
    for aCPrime in cprimes:
        for aBRnew in brnews:

            theID = theIDBase.format(str(aCPrime).replace(".",""),str(aBRnew).replace(".",""))
            print "# RUNNING ON POINT C'^2 = ",aCPrime,", BRnew = ",aBRnew,"#"
            print theID
            
            if aCPrime/(1-aBRnew) > 1.:
                print "WARNING: This point is excluded from the search -> ",aCPrime/(1-aBRnew),"> 1"
                continue
            
            hcgSubDir = hcgDir+"/HCG_"+theID+opt.append
            dcSubDir = "cards_{0}TeV_"+theID+opt.append+"_tagged/HCG/{1}/*"
            
            if aCPrime==1. and aBRnew==0.:
                if (opt.create): createDCs(aCPrime, aBRnew)
                copyDCsToHCG(dcDir,dcSubDir,hcgSubDir)
                copyScripts(hcgSubDir)                
                os.chdir(hcgSubDir)
                combineDCsAndWSs()            
            else:
                inputDir = refDir
                os.makedirs(hcgSubDir)
                os.chdir(hcgSubDir)
                copyScripts(hcgSubDir)
                massFileBuf = open(massFile,"r")            
                for aMass in massFileBuf:
                    theMassStr =  str(aMass.split()[0])
                    if theMassStr == "": break
                    os.makedirs(theMassStr)
                    os.chdir(theMassStr)
                    processCmd("ln -fs "+inputDir+"/"+theMassStr+"/comb_hzz4l.txt .")
                    processCmd("ln -fs "+inputDir+"/"+theMassStr+"/comb_hzz4l.root .")
                    os.chdir(hcgSubDir)
                
                os.chdir(dcDir)

            runTheLimits(hcgSubDir,aCPrime,aBRnew)

            

