import optparse, re
import sys, os, pwd, commands
import subprocess
import shlex
import shutil, errno
import glob


dcDir = "/afs/cern.ch/work/s/scasasso/private/HZZ4L_Combination_Git_forHighMass/HZZ4L_CombinationPy/CreateDatacards"
hcgDir = dcDir+"/CMSSW_6_1_1/src"
massFile = dcDir+"/masses_HMP_above400.txt"
#massFile = dcDir+"/masses_test.txt"
scriptsDir = "/afs/cern.ch/work/s/scasasso/private/HZZ4L_Combination_Git_forHighMass/HZZ4L_CombinationPy/RunLimits/scripts"



def processCmd(cmd):

    args = shlex.split(cmd)
    sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()

    print err
    
    return out, err



def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy2(src, dst)
        else: raise



def runTheGrid():
    
    for aCPrime in [0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
        for aBRnew in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
#     for aCPrime in [0.1]:
#         for aBRnew in [0.0]:
            
            if aCPrime/(1-aBRnew)>1.:
                print "WARNING: skipping this point because Gamma > Gamma_SM ",aCPrime," / ",aBRnew
                exit

            print "\n"
            print "Running point C'^2 = ",aCPrime,", BRnew = ",aBRnew

            theID = "cpsq"+str(aCPrime).replace(".","")+"_brnew"+str(aBRnew).replace(".","")
            print theID

            os.chdir(dcDir)
            # Preparing the datacards
            print "   Creating 8 TeV datacards ..."
            out, err = processCmd("python makeDCsandWSs.py -i SM_inputs_8TeV_HMP_onlyggH -a 8TeV_"+str(theID)+" -b -j0 -c"+str(aCPrime)+" -r"+str(aBRnew))
            print "   Creating 7 TeV datacards ..."
            out, err = processCmd("python makeDCsandWSs.py -i SM_inputs_7TeV_HMP_onlyggH -a 7TeV_"+str(theID)+" -b -j0 -c"+str(aCPrime)+" -r"+str(aBRnew))

            oStream = open('/tmp/scasasso/stdout_'+theID+'.txt','w')
            oStream.write(out)
            eStream = open('/tmp/scasasso/stderr_'+theID+'.txt','w')
            eStream.write(err)

            oStream.close()
            eStream.close()            

            # Merging the datacards
            hcgSubDir = hcgDir+"/HCG_"+theID            
            
            print "   Copying 8 TeV datacards directory ..."
            for aFile in glob.glob(r"cards_8TeV_"+theID+"/HCG/*"):
                index = len(aFile.split("/"))-1
                theMass = aFile.split("/")[index].strip()
                shutil.copytree(aFile, hcgSubDir+"/"+theMass)
            os.chdir(hcgSubDir)
            shutil.copy(massFile,".")
            
            print "   Copying 7 TeV datacards directory ..."
            massFileBuf = open(massFile,"r")
            for aMass in massFileBuf:
                theMassStr =  str(aMass.split()[0])
                for aFile in glob.glob(r""+dcDir+"/cards_7TeV_"+theID+"/HCG/"+theMassStr+"/*"):
                    shutil.copy2(aFile, hcgSubDir+"/"+theMassStr+"/")

            for anExt in ["h","cxx","py","cc","sh","C"]:
                for aFile in glob.glob(r""+scriptsDir+"/*."+anExt): shutil.copy2(aFile,".")
            shutil.copy2(scriptsDir+"/hadd2",".")
            shutil.copy2(scriptsDir+"/massrange",".")            
            shutil.copytree(scriptsDir+"/cms","./cms/")
            shutil.copytree(scriptsDir+"/plotting","./plotting/")
            shutil.copytree(scriptsDir+"/grids","./grids/")            
	
	    #Combining the datacards and the workspaces
            print "   Combining datacards and workspaces ..."
            outComb, errComb = processCmd("./makeCardsWorkspaces.sh "+massFile+" 0")
            print outComb

	    #Running the asymptotic limits
            print "   Running the limits on lxbatch ..."
            processCmd('python make_parallel_limits.py -M ASCLS -f '+massFile+' -q lsf -t combine -o "-E"')

        

if __name__ == "__main__":
    runTheGrid()
    
