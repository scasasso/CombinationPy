#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#endif

//Set Input Parameters
#include "Config.h"

//User determined Parameters
const TString destDir = "../CreateDatacards/templates2D/";
int useSqrts=0;                                                   // 0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only

TString samples[9]={"ggF","VBF","qqZZ","ggZZ","Z+X","ZH","WH","ttH","VBF Bkg"};
TString lowmassstr[2]={"Low Mass Region","High Mass Region"};
TString updownstr[5]={"JEC down","Nominal","JEC up","Alt shape","Alt shape 2"};

TFile* fggH,*fqqH,*fqqZZ,*fggZZ,*fZX,*fZH,*fWH,*fttH,*fVBFBkg;

//Global Parameters (not tested if cause issues if altered)
bool extendToHighMass = true;                                     // Include signal samples above 600 GeV
float highMzz=(extendToHighMass?1600:800);
float mBinSize=2.;

//Function Constructors 
void templateOptions(bool debug, bool findAlternatives);          // Sets debug mode (only produces unsmoothed plots) and whether alternatives are produced
void combineTemplates(bool moreAlt);                              // Combines templates into format for analysis
void buildChain(TChain* bkgMC, int sampleIndex);                  // Generates TChain from input MC
void makeTemplate(int updown,bool debug);                         // Makes template root files - NOTE: ggH, ggZZ, qqZZ, and Z+X are listed explicitly and do not rely on Config.h
TH2F* fillTemplate(int sampleIndex, bool isLowMass, int updown);  // Takes input MC to make Fisher v m4L TH2F
TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp);              // Merges low and high mass TH2F plots
TH2F* smoothtemplates(TH2F* inputdata, int sampleIndex);          // Smooths Fisher v m4L plots
TH2F* rebin(TH2F* rebinnedhist, int sampleIndex);                 // Smoothing procedure for ggH,qqH
TH2F* rebin_lowstatistics(TH2F* rebinnedhist, int sampleIndex);   // Smoothing procedure for ggZZ,qqZZ,ZH,WH
TH2F* altshapes(TH2F* originalHist, int channel, int altnum);     // Generates alternative shapes for Fisher, (channel,altnum): (0,1)=ggHAMCNLO, (0,2)=ggHMG,
                                                                  //               (1,1)=qqHd6t, (1,2)= qqHatlas, (2,1)=qqZZMG, (2,2)=mirror of qqZZMG
float altscale(float Fisher, int channel, int altnum);            // Find scale for alternative shapes given a Fisher value, channel and altnum same as altshapes()
TH2F* mirrortemplates(int sampleIndex);                           // Generates mirror alternative shapes for Fisher when no second alternative exists
bool test_bit(int mask, unsigned int iBit);                       // Used to identify correct CR events
void progressbar(int val, int tot);                               // Output progress bar
TH2F* makebkg(TH2F *SM1, TH2F *SM10, TH2F *SM25);                 // Subroutine for VBF qqZZ background, will be replaced when the full-sim come through

//---------------------------------------------------

bool test_bit( int mask, unsigned int iBit ) { return (mask >> iBit) & 1; }

void generatevbfMELATemplates() {
  bool debug=false;
  bool findAlternatives=true;
  templateOptions(debug,findAlternatives);
}
 
void templateOptions(bool debug, bool findAlternatives){
  makeTemplate(0,debug);
  if (findAlternatives){
    makeTemplate(1,debug);
    makeTemplate(-1,debug);
    makeTemplate(2,debug);
    makeTemplate(3,debug);
    combineTemplates(true);
  }
  if(!findAlternatives){
    combineTemplates(false);
  }
}

void combineTemplates(bool moreAlt){
  TString sig[9] = {"ggH","qqH","qqZZ","ZH","WH","ttH","ggZZ","Z+X","VBFBkg"};
  for (int i = 0; i < 9; i++)
    {

      TString debugname;  
      if (!moreAlt) debugname="_unnormalized";

      TFile* fishDef = new TFile(destDir+sig[i]+"_vbfMELA"+debugname+".root", "OPEN");
      TH2F* Default = (TH2F*)fishDef->Get("H_Djet");
      
      TFile* fishUp;
      if(moreAlt){
        if(i==0 || i==1 || i==2) fishUp = new TFile(destDir+sig[i]+"_vbfMELA_alt.root","OPEN");
        else if(i==3 || i==4 || i==5 || i==6 || i==8) fishUp = new TFile(destDir+sig[i]+"_vbfMELA_up.root","OPEN");
        else if(i==7) fishUp = new TFile(destDir+sig[i]+"_vbfMELA.root","OPEN");
      }
      else{
        fishUp = new TFile(destDir+sig[i]+"_vbfMELA"+debugname+".root","OPEN");
      }
      TH2F* Up = (TH2F*)fishUp->Get("H_Djet");
      
      TFile* fishDn;
      if(moreAlt){
        if(i==0 || i==1 || i==2) fishDn = new TFile(destDir+sig[i]+"_vbfMELA_alt2.root","OPEN");
        else if(i==3 || i==4 || i==5 || i==6 || i==8) fishDn = new TFile(destDir+sig[i]+"_vbfMELA_down.root","OPEN");
        else if(i==7) fishDn = new TFile(destDir+sig[i]+"_vbfMELA.root","OPEN");
      }
      else{
        fishDn = new TFile(destDir+sig[i]+"_vbfMELA"+debugname+".root","OPEN");
      }      
      TH2F* Down = (TH2F*)fishDn->Get("H_Djet");

      TFile* Out = new TFile(destDir+sig[i]+"_fisher"+debugname+".root","RECREATE");
      Out->cd();

      Default->Write("h_Fisher");
      Up->Write("h_Fisher_up");
      Down->Write("h_Fisher_dn");
  
      Out->Close();
      fishDef->Close();
      fishUp->Close();
      fishDn->Close();
    }
}

//---------------------------------------------------

void buildChain(TChain* bkgMC, int sampleIndex){
  TString filePath;
  int nPoints=0;
  int masses[100];
  if (useSqrts==1){
    if (sampleIndex==0){
      nPoints=nPoints7TeV_MINLO;
      for (int i=0;i<nPoints;i++){
      	masses[i]=masses7TeV_MINLO[i];
      }
    } else if (sampleIndex==1){
      nPoints=nVBFPoints7TeV;
      for (int i=0;i<nPoints;i++){
      	masses[i]=VBFmasses7TeV[i];
      }
    } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
      nPoints=nVHPoints7TeV;
      for (int i=0;i<nPoints;i++){
      	masses[i]=VHmasses7TeV[i];
      }
    }
    filePath = filePath7TeV;
  }
  if (useSqrts==2){
    if (sampleIndex==0){
      nPoints=nPoints8TeV_MINLO;
      for (int i=0;i<nPoints;i++){
      	masses[i]=masses8TeV_MINLO[i];
      }
    } else if (sampleIndex==1){
      nPoints=nVBFPoints8TeV;
      for (int i=0;i<nPoints;i++){
      	masses[i]=VBFmasses8TeV[i];
      }
    } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
      nPoints=nVHPoints8TeV;
      for (int i=0;i<nPoints;i++){
      	masses[i]=VHmasses8TeV[i];
      }
    }
    filePath = filePath8TeV;
  }
  if (useSqrts==0){
    if (sampleIndex==0){
      nPoints=nPoints7TeV_MINLO+nPoints8TeV_MINLO;
      for (int i=0;i<nPoints7TeV_MINLO;i++){
      	masses[i]=masses7TeV_MINLO[i];
      }
      for (int i=0;i<nPoints8TeV_MINLO;i++){
      	masses[i+nPoints7TeV_MINLO]=masses8TeV_MINLO[i];
      }
    } else if (sampleIndex==1){
      nPoints=nVBFPoints7TeV+nVBFPoints8TeV;
      for (int i=0;i<nVBFPoints7TeV;i++){
      	masses[i]=VBFmasses7TeV[i];
      }
      for (int i=0;i<nVBFPoints8TeV;i++){
      	masses[i+nVBFPoints7TeV]=VBFmasses8TeV[i];
      }
    } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
      nPoints=nVHPoints7TeV+nVHPoints8TeV;
      for (int i=0;i<nVHPoints7TeV;i++){
      	masses[i]=VHmasses7TeV[i];
      }
      for (int i=0;i<nVHPoints8TeV;i++){
      	masses[i+nVHPoints7TeV]=VHmasses8TeV[i];
      }
    }
  }
  if(sampleIndex!=2 && sampleIndex!=3 && sampleIndex!=4 && sampleIndex!=-8 && sampleIndex!=8 && sampleIndex!=9 && sampleIndex!=10){
    if(nPoints==0){
      cout<<"nPoints not set in Config.h"<<endl;
      abort();
    }
    for (int i=0; i<nPoints; i++){
      char tmp_finalInPath4mu[200],tmp_finalInPath4e[200],tmp_finalInPath2mu2e[200];
      string finalInPath4mu,finalInPath4e,finalInPath2mu2e;
      if (sampleIndex==0){
        sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_minloH%i.root",masses[i]);
        sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_minloH%i.root",masses[i]);
        sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_minloH%i.root",masses[i]);
      }else if (sampleIndex==1){
        if(masses[i]<200){
      	 sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_VBFH%i.root",masses[i]);
      	 sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_VBFH%i.root",masses[i]);
      	 sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_VBFH%i.root",masses[i]);
        } else{
         sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_powheg15VBFH%i.root",masses[i]);
         sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_powheg15VBFH%i.root",masses[i]);
         sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_powheg15VBFH%i.root",masses[i]);
        }
      }else if (sampleIndex==5){
      	sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_ZH%i.root",masses[i]);
      	sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_ZH%i.root",masses[i]);
      	sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_ZH%i.root",masses[i]);
      }else if (sampleIndex==6){
      	sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_WH%i.root",masses[i]);
      	sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_WH%i.root",masses[i]);
      	sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_WH%i.root",masses[i]);
      }else if (sampleIndex==7){
      	sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_ttH%i.root",masses[i]);
      	sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_ttH%i.root",masses[i]);
      	sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_ttH%i.root",masses[i]);
      }
      if (useSqrts!=0){
      	finalInPath4mu = filePath + tmp_finalInPath4mu;
      	finalInPath4e = filePath + tmp_finalInPath4e;
      	finalInPath2mu2e = filePath + tmp_finalInPath2mu2e;
      } else if (useSqrts==0){
      	if ((sampleIndex==0 && i<nPoints7TeV_MINLO) || (sampleIndex==1 && i<nVBFPoints7TeV) || ((sampleIndex==5 || sampleIndex==6 || sampleIndex==7) && i<nVHPoints7TeV)){
      	  finalInPath4mu = filePath7TeV + tmp_finalInPath4mu;
      	  finalInPath4e = filePath7TeV + tmp_finalInPath4e;
      	  finalInPath2mu2e = filePath7TeV + tmp_finalInPath2mu2e;
      	} else{
      	  finalInPath4mu = filePath8TeV + tmp_finalInPath4mu;
      	  finalInPath4e = filePath8TeV + tmp_finalInPath4e;
      	  finalInPath2mu2e = filePath8TeV + tmp_finalInPath2mu2e;
      	}
      }
      bkgMC->Add(finalInPath4mu.c_str());
      bkgMC->Add(finalInPath4e.c_str());
      bkgMC->Add(finalInPath2mu2e.c_str());
    }
  }
  else if (sampleIndex==2){
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo4tau.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo2e2tau.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo4tau.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo2e2tau.root");	
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo4tau.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo2e2tau.root");
    }	
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo2e2mu.root");
      //bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo4tau.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo2e2tau.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo2e2mu.root");
      //bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo4tau.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo2e2tau.root");	
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo4tau.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo2e2tau.root");
    }
  }
  else if (sampleIndex==3){
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ggZZ4l.root");
    }
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ggZZ4l.root");
    }
  }
  else if (sampleIndex==4){
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleOr_CRZLLTree.root");
    }
    //NEEDS TO BE CHANGED BEFORE PUSHING TO OTHERS
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleOr_CRZLLTree.root");
    }
  }
  else if (sampleIndex==-8){
    bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZJetsTo4L.root");
    bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZJetsTo4L.root");
    bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZJetsTo4L.root");
  }
  else if(sampleIndex==8){ //Will be updated if there is a final combination of all CJLST samples
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2muJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo4eJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo4muJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo2e2muJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4eJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4muJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo2e2muJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4eJJ_Contin.root");
    bkgMC->Add("/scratch0/hep/ianderso/CJLST/140604/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4muJJ_Contin.root");
  }   
}


//---------------------------------------------------

void makeTemplate(int updown, bool debug){

  TString jes;
  if (updown==1){
    jes="_up";
  }
  else if(updown==-1){
    jes="_down";
  }
  else if(updown==2){
    jes="_alt";
  }
  else if(updown==3){
    jes="_alt2";
  }

  TString debugname;
  if (debug) debugname="_unnormalized";

  if(updown==0 || updown==1 || updown==-1){
    fqqH = new TFile(destDir + "qqH_vbfMELA"+jes+debugname+".root","RECREATE");
    fggH = new TFile(destDir + "ggH_vbfMELA"+jes+debugname+".root","RECREATE");
    fqqZZ = new TFile(destDir + "qqZZ_vbfMELA"+jes+debugname+".root","RECREATE");
    fggZZ = new TFile(destDir + "ggZZ_vbfMELA"+jes+debugname+".root","RECREATE");
    fZX = new TFile(destDir + "Z+X_vbfMELA"+jes+debugname+".root","RECREATE");
    fZH = new TFile(destDir + "ZH_vbfMELA"+jes+debugname+".root","RECREATE");
    fWH = new TFile(destDir + "WH_vbfMELA"+jes+debugname+".root","RECREATE");
    fttH = new TFile(destDir + "ttH_vbfMELA"+jes+debugname+".root","RECREATE");
    fVBFBkg = new TFile(destDir + "VBFBkg_vbfMELA"+jes+debugname+".root","RECREATE");
  } else{
    fqqH = new TFile(destDir + "qqH_vbfMELA"+jes+debugname+".root","RECREATE");
    fggH = new TFile(destDir + "ggH_vbfMELA"+jes+debugname+".root","RECREATE");
    fqqZZ = new TFile(destDir + "qqZZ_vbfMELA"+jes+debugname+".root","RECREATE");
  }

  TH2F* low,*high,*H_Djet; 
  
  // =========================
  // ggH
  
  low = fillTemplate(0,true,updown);
  high = fillTemplate(0,false,updown);
  H_Djet = mergeTemplates(low,high);

  if (!debug && updown<2) smoothtemplates(H_Djet,0);
  if (!debug && updown==2) smoothtemplates(H_Djet,-1);
  if (!debug && updown==3) smoothtemplates(H_Djet,-2);


  fggH->cd();
  H_Djet->Write("H_Djet");
  fggH->Close();
  
  // ==========================
  // qqH
  
  low = fillTemplate(1,true,updown);
  high = fillTemplate(1,false,updown);
  H_Djet = mergeTemplates(low,high);  

  if (!debug && updown<2) smoothtemplates(H_Djet,1);
  if (!debug && updown==2) smoothtemplates(H_Djet,-3);
  //Test for p-values <- Go to smoothtemplates to change
  if (!debug && updown==3) smoothtemplates(H_Djet,-4);

  fqqH->cd();
  H_Djet->Write("H_Djet");
  fqqH->Close();
  
  
  // ==========================
  // qqZZ

  low = fillTemplate(2,true,updown);
  high = fillTemplate(2,false,updown);
  H_Djet = mergeTemplates(low,high);

  if (!debug) smoothtemplates(H_Djet,2);
  if (!debug && updown==2) altshapes(H_Djet,2,1);
  if (!debug && updown==3) altshapes(H_Djet,2,2);
  if (!debug && updown>1){
    TH1F* tempProj;
    for(int i=1; i<=H_Djet->GetNbinsX(); i++){
      tempProj = (TH1F*) H_Djet->ProjectionY("tempProj",i,i);
      float norm=tempProj->Integral();
      if (norm>0) { // Avoid introducing NaNs in the histogram
      	for(int j=1; j<=H_Djet->GetNbinsY(); j++){
      	  H_Djet->SetBinContent(i,j, H_Djet->GetBinContent(i,j)/norm);
      	}
      }
    }
  }
  //Test for p-values
  //if (!debug && updown<3) smoothtemplates(H_Djet,2);
  //if (!debug && updown==3) mirrortemplates(2);

  fqqZZ->cd();
  H_Djet->Write("H_Djet");
  fqqZZ->Close();
  
  if (updown==0 || updown==1 || updown==-1){
    // ==========================
    // ggZZ
    
    low = fillTemplate(3,true,updown);
    high = fillTemplate(3,false,updown);
    H_Djet = mergeTemplates(low,high);
    
    if (!debug) smoothtemplates(H_Djet,3);
    
    fggZZ->cd();
    H_Djet->Write("H_Djet");
    fggZZ->Close();
    
    // ==========================
    // Z+X
    
    low = fillTemplate(4,true,updown);
    high = fillTemplate(4,false,updown);
    H_Djet = mergeTemplates(low,high);
    
    if(!debug) smoothtemplates(H_Djet,4);
    
    fZX->cd();
    H_Djet->Write("H_Djet");
    fZX->Close();
    
    // ==========================
    // ZH

    low = fillTemplate(5,true,updown);
    high = fillTemplate(5,false,updown);
    H_Djet = mergeTemplates(low,high);
    
    if (!debug) smoothtemplates(H_Djet,5);
    
    fZH->cd();
    H_Djet->Write("H_Djet");
    fZH->Close();
    
    // ==========================
    // WH
    
    low = fillTemplate(6,true,updown);
    high = fillTemplate(6,false,updown);
    H_Djet = mergeTemplates(low,high);
    
    if (!debug) smoothtemplates(H_Djet,6);
    
    fWH->cd();
    H_Djet->Write("H_Djet");
    fWH->Close();
    
    // ==========================
    // ttH
    
    low = fillTemplate(7,true,updown);
    high = fillTemplate(7,false,updown);
    H_Djet = mergeTemplates(low,high);
    
    if (!debug) smoothtemplates(H_Djet,7);

    fttH->cd();
    H_Djet->Write("H_Djet");
    fttH->Close();
    
    // ==========================
    // VBF Bkg
    
    low = fillTemplate(8,true,updown);
    high = fillTemplate(8,false,updown);
    H_Djet = mergeTemplates(low,high);

    if (!debug) smoothtemplates(H_Djet,8);

    fVBFBkg->cd();
    H_Djet->Write("H_Djet");
    fVBFBkg->Close();        
  }
}

//---------------------------------------------------

TH2F* fillTemplate(int sampleIndex,bool isLowMass,int updown){
  TChain* bkgMC = new TChain("SelectedTree");
  buildChain(bkgMC,sampleIndex);

  cout << "Chain for " << samples[sampleIndex] << ", " << lowmassstr[isLowMass] << ", " << updownstr[updown+1] << " with  " << bkgMC->GetEntries() << " events." << endl;
  bkgMC->ls();

  float mass,w,phjj,pvbf,Djet,mJJ;
  Short_t njets;
  int processID;
  int CRflag;
  string channel;

  bkgMC->SetBranchAddress("ZZMass",&mass);
  bkgMC->SetBranchAddress("NJets30",&njets);
  if(updown==0 || updown==2 || updown==3){
    bkgMC->SetBranchAddress("pvbf_VAJHU_old",&pvbf);
    bkgMC->SetBranchAddress("phjj_VAJHU_old",&phjj);
  }
  else if(updown==1){
    bkgMC->SetBranchAddress("pvbf_VAJHU_old_up",&pvbf);
    bkgMC->SetBranchAddress("phjj_VAJHU_old_up",&phjj);
  }
  else if(updown==-1){
    bkgMC->SetBranchAddress("pvbf_VAJHU_old_dn",&pvbf);
    bkgMC->SetBranchAddress("phjj_VAJHU_old_dn",&phjj);
  }
  bkgMC->SetBranchAddress("MC_weight",&w);
  bkgMC->SetBranchAddress("DiJetMass",&mJJ);
  bkgMC->SetBranchAddress("genProcessId",&processID);
  if(sampleIndex==4){
    bkgMC->SetBranchAddress("CRflag",&CRflag);
  }

  TH2F* bkgHist;
  if(!isLowMass){
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-180.)/mBinSize+0.5),180,highMzz,50,0,1);
  }
  else{
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((180-100)/mBinSize+0.5),100,180,50,0,1);
  }

  //bkgHist->Sumw2();

  //float peakw=0.;

  //Fill histogram
  for(int i=0; i<bkgMC->GetEntries(); i++){
    bkgMC->GetEntry(i);
    progressbar(i,bkgMC->GetEntries());
    if((sampleIndex==4 && (test_bit(CRflag,5) || test_bit(CRflag,7) || test_bit(CRflag,9) || test_bit(CRflag,11))) || sampleIndex!=4){
      if (mass<100 || (sampleIndex!=4 && (njets<2 || phjj==-1. || pvbf==-1.)) || (sampleIndex==4 && (phjj==-1. || pvbf==-1.))) continue;
      if(sampleIndex>7 && mJJ<130.) continue;
      if(sampleIndex==0 && w<0) continue;
      //if(updown!=0) Fisher = 0.18*fabs(deta) + 0.000192*mJJ;
      Djet=pvbf/(pvbf+phjj);
      bkgHist->Fill(mass,Djet,w);
      //if(sampleIndex>7 && mass>120. && mass<130.) peakw+=w; 
    }
  }
  cout<<endl;

  //cout<<"SUM OF WEIGHTS IN PEAK: "<<peakw<<endl;

  return bkgHist;
}

//---------------------------------------------------

TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp){

  int nYbins=lowTemp->GetNbinsY();
  if (highTemp->GetNbinsY()!=nYbins) {
    cout << "ERROR: mergeTemplates: incorrect binning " << endl;
    abort();
  }

  TH2F* H_Djet = new TH2F("H_Djet","H_Djet",int((highMzz-100.)/mBinSize +0.5),100,highMzz,nYbins,0,1);

  //H_Djet->Sumw2();

  // copy lowmass into H_Djet
  for(int i=1; i<=lowTemp->GetNbinsX(); ++i){
    for(int j=1; j<=nYbins; ++j){
      H_Djet->SetBinContent(i,j, lowTemp->GetBinContent(i,j)  );
    }// end loop over Fisher
  }// end loop over m4L

  // copy high mass into H_Djet
  for(int i=1; i<=highTemp->GetNbinsX(); ++i){
    for(int j=1; j<=nYbins; ++j){
      H_Djet->SetBinContent(i+lowTemp->GetNbinsX(),j, highTemp->GetBinContent(i,j)  );
    }// end loop over Fisher
  }// end loop over m4L

  return H_Djet;
}

//---------------------------------------------------

TH2F* smoothtemplates(TH2F* inputdata, int sampleIndex){
  if(sampleIndex<3 || sampleIndex>7){
    rebin(inputdata,sampleIndex);
  }
  else if(sampleIndex>2 && sampleIndex<8){
    rebin_lowstatistics(inputdata,sampleIndex);
  }

  return inputdata;
}

//---------------------------------------------------

TH2F* rebin(TH2F* rebinnedHist, int usealt){

  int nXbins=rebinnedHist->GetNbinsX();
  int nYbins=rebinnedHist->GetNbinsY();

  double norm;
  TH1F* tempProj;

  rebinnedHist->Sumw2();

  int negativecolumns=0;
  int negativebins2=0;

  //Normalization
  for(int i=1; i<=nXbins; i++){
    tempProj = (TH1F*) rebinnedHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    for(int j=1; j<=nYbins; j++) if(rebinnedHist->GetBinContent(i,j)<0.) negativebins2++;
    if (norm!=0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
        rebinnedHist->SetBinContent(i,j, rebinnedHist->GetBinContent(i,j)/norm   );
      }
    }
    if(norm<0){
      cout<<i<<endl;
      negativecolumns++;
    }
  }

  cout<<"Initial bins, after normalization, with negative values: "<<negativebins2<<"/"<<nXbins*nYbins<<endl;
  cout<<"Initial columns, after normalization, with negative values: "<<negativecolumns<<"/"<<nXbins<<endl;

  TH2F* origHist = new TH2F (*rebinnedHist);

  origHist->Sumw2();
  //rebinnedHist->Sumw2();

  int effectiveArea=1;
  double average=0,binsUsed=0;

  //Rebin m4L
  for (int i=1; i <=nXbins; i++){
    for (int j=1; j<=nYbins; j++){
      float binMzz = rebinnedHist->GetBinCenter(i);

      if( binMzz<300 ) continue;
      if( binMzz>=300 && binMzz<350 ) effectiveArea=1;
      if( binMzz>=350 && binMzz<500 ) effectiveArea=3;
      if( binMzz>=500 && binMzz<600 ) effectiveArea=5;
      if( binMzz>=600 && binMzz<800 ) effectiveArea=7;
      if( binMzz>=800 && binMzz<1000) effectiveArea=11;
      if( binMzz>=1000 && binMzz<1200)effectiveArea=15;
      if( binMzz>=1200 && binMzz<1400) effectiveArea=25;
      if( binMzz>1400) effectiveArea=40;
      
      for(int a=-effectiveArea; a<=effectiveArea; a++){
      	if(a+i<1 || a+i>nXbins || j>nYbins || j<1) continue;
      	average+=origHist->GetBinContent(a+i,j);
      	binsUsed++;
      }
      rebinnedHist->SetBinContent(i,j,average/binsUsed);
      average=0;
      binsUsed=0;
    }
  }


  TH2F* Histstg1 = new TH2F (*rebinnedHist);

  Histstg1->Sumw2();

  //Use average of Nearest Neighbors to fill remaining zeroes
  for(int i=1; i<=nXbins;i++){
    for(int j=1; j<=nYbins;j++){
      float binvalue = rebinnedHist->GetBinContent(i,j);
      if(binvalue!=0) continue;
      for (int i2=-1;i2<=1;i2++){
      	if (i2+i<1 || i2+i>nXbins || j<1 || j>nYbins) continue;
      	average+=rebinnedHist->GetBinContent(i+i2,j);
      	binsUsed++;
      }
      for (int j2=-1;j2<=1;j2++){
      	if (i<1 || i>nXbins || j2+j<1 || j2+j>nYbins) continue;
      	average+=rebinnedHist->GetBinContent(i,j+j2);
      	binsUsed++;
      }
      rebinnedHist->SetBinContent(i,j,average/binsUsed);
      average=0;
      binsUsed=0;
    }
  }

  if (usealt==-1) altshapes(rebinnedHist,0,1);
  if (usealt==-2) altshapes(rebinnedHist,0,2);
  if (usealt==-3) altshapes(rebinnedHist,1,1);
  if (usealt==-4) altshapes(rebinnedHist,1,2);
  if (usealt==-5) altshapes(rebinnedHist,2,1);
  if (usealt==-6) altshapes(rebinnedHist,2,2);  

  int negativebins=0;

  float averagebin=0.;
  int totbins=0;
  for(int i=1; i<=nXbins; i++){
    for(int j=1; j<=nYbins; j++){
      averagebin+=rebinnedHist->GetBinContent(i,j);
      totbins++;
    }
  }
  for(int i=1;i<=nXbins;i++){
    for(int j=1; j<=nYbins; j++){
      if(rebinnedHist->GetBinContent(i,j)<0.00001*averagebin) rebinnedHist->SetBinContent(i,j,0.00001*averagebin);
      if(rebinnedHist->GetBinContent(i,j)<0.) negativebins++;
    }
  }

  cout<<"Final bins with negative values: "<<negativebins<<"/"<<nXbins*nYbins<<endl;

  //Renormalize
  for(int i=1; i<=nXbins; i++){
    tempProj = (TH1F*) rebinnedHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm!=0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
      	rebinnedHist->SetBinContent(i,j, rebinnedHist->GetBinContent(i,j)/norm   );
      }
    }
  }

  return rebinnedHist;
}

//---------------------------------------------------

TH2F* rebin_lowstatistics(TH2F* finalhist, int sampleIndex){
  int nXbins=finalhist->GetNbinsX();
  int nYbins=finalhist->GetNbinsY();

  TH2F* origHist = new TH2F (*finalhist);

  //Project Fisher v m4L to 1D Fisher plots for low, high, and full mass range
  TH1F* lowProj = (TH1F*) origHist->ProjectionY("lowProj",1,40);
  TH1F* highProj = (TH1F*) origHist->ProjectionY("highProj",41,750);
  TH1F* fullProj = (TH1F*) origHist->ProjectionY();
  double temp;


  //Plots if need to see tails
  /*if(sampleIndex!=2 && sampleIndex!=3 && sampleIndex!=5){
    TCanvas* clow = new TCanvas("clow","clow",800,800);
    clow->cd();
    lowProj->Draw();
    
    TCanvas* chigh = new TCanvas("chigh","chigh",800,800);
    chigh->cd();
    highProj->Draw();

    TCanvas* cfull = new TCanvas("cfull","cfull",800,800);
    cfull->cd();
    fullProj->Draw();
    }*/

  //qqZZ + WH + ZH
  if(sampleIndex==2 || sampleIndex==5 || sampleIndex==6 || sampleIndex==8){
    //Tail smoothing goes here if needed

    //Fill each mass point with low or high mass projections
    for(int i=1; i<=nXbins;i++){
      for(int j=1; j<=nYbins;j++){
      	float binMzz = finalhist->GetBinCenter(i);
      	if (binMzz<180){
      	  if(sampleIndex!=8) temp=lowProj->GetBinContent(j);
          else if(sampleIndex==8) temp=highProj->GetBinContent(j);
      	  finalhist->SetBinContent(i,j,temp);
      	}
      	else if (binMzz>=180){
      	  temp=highProj->GetBinContent(j);
      	  finalhist->SetBinContent(i,j,temp);
      	}
      }
    }

    //Store plots of fits, in case anything goes wrong
    if(sampleIndex==2){
      fqqZZ->cd();
    }
    else if(sampleIndex==5){
      fZH->cd();
    }
    else if(sampleIndex==6){
      fWH->cd();
    }
    else if(sampleIndex==8){
      fVBFBkg->cd();
    }
    lowProj->Write("H_low");
    highProj->Write("H_high");
    
  }
  //ggZZ + ttH + Z+X
  else if(sampleIndex==3 || sampleIndex==4 || sampleIndex==7){    
    //Fill each mass point with full projection
    for(int i=1; i<=nXbins;i++){
      for(int j=1; j<=nYbins;j++){
      	temp=fullProj->GetBinContent(j);
      	finalhist->SetBinContent(i,j,temp);
      }
    }

    //Store plots of fits, in case anything goes wrong
    if(sampleIndex==3){
      fggZZ->cd();
    }
    if(sampleIndex==4){
      fZX->cd();
    }
    if(sampleIndex==7){
      fttH->cd();
    }
    fullProj->Write("H_full");

  }

  double norm;
  TH1F* tempProj;

  //Normalize
  for(int i=1; i<=nXbins; i++){
    tempProj = (TH1F*) finalhist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
      	finalhist->SetBinContent(i,j, finalhist->GetBinContent(i,j)/norm);
      }
    }
  }

  return finalhist;
}

TH2F* altshapes(TH2F* original, int channel, int altnum){	 	 
  int nXbins=original->GetNbinsX();
  int nYbins=original->GetNbinsY();
  float binDjet,scale,altvalue;
	 	 
  for (int i=1; i <=nXbins; i++){
    for (int j=1; j<=nYbins; j++){
      binDjet = original->GetYaxis()->GetBinCenter(j);
      scale = altscale(binDjet,channel,altnum);
      altvalue=scale*(original->GetBinContent(i,j));
      original->SetBinContent(i,j,altvalue);
    }
  }
	 	 
  return original;	 	
}
	 	 
float altscale(float Djet,int channel, int altnum){
  float scale=0;

  //ggH
  if(channel==0){
    if(altnum==1){
      if (Djet<=0.1) scale=0.937288;
      if (Djet>0.1 && Djet<=0.2) scale=0.978014;
      if (Djet>0.2 && Djet<=0.3) scale=0.961835;
      if (Djet>0.3 && Djet<=0.4) scale=1.02355;
      if (Djet>0.4 && Djet<=0.5) scale=1.05243;
      if (Djet>0.5 && Djet<=0.6) scale=1.01841;
      if (Djet>0.6 && Djet<=0.7) scale=1.06476;
      if (Djet>0.7 && Djet<=0.8) scale=1.11512;
      if (Djet>0.8 && Djet<=0.9) scale=1.2112;
      if (Djet>0.9) scale=1.54671;
    }
    if(altnum==2){
      if (Djet<=0.1) scale=0.960535;
      if (Djet>0.1 && Djet<=0.2) scale=1.1733;
      if (Djet>0.2 && Djet<=0.3) scale=1.11984;
      if (Djet>0.3 && Djet<=0.4) scale=1.16218;
      if (Djet>0.4 && Djet<=0.5) scale=1.07156;
      if (Djet>0.5 && Djet<=0.6) scale=0.993379;
      if (Djet>0.6 && Djet<=0.7) scale=0.950683;
      if (Djet>0.7 && Djet<=0.8) scale=0.891218;
      if (Djet>0.8 && Djet<=0.9) scale=0.762706;
      if (Djet>0.9) scale=0.633827;
    }
  }
  //qqH
  if(channel==1){
    if(altnum==1){
      if (Djet<=0.1) scale=1.15736;
      if (Djet>0.1 && Djet<=0.2) scale=1.1752;
      if (Djet>0.2 && Djet<=0.3) scale=1.20939;
      if (Djet>0.3 && Djet<=0.4) scale=1.01484;
      if (Djet>0.4 && Djet<=0.5) scale=1.15581;
      if (Djet>0.5 && Djet<=0.6) scale=1.01961;
      if (Djet>0.6 && Djet<=0.7) scale=1.0057;
      if (Djet>0.7 && Djet<=0.8) scale=0.961047;
      if (Djet>0.8 && Djet<=0.9) scale=0.97326;
      if (Djet>0.9) scale=0.914909;
      }
    /* Atlas tuning
    if(altnum==2){
      if (Djet<=0.1) scale=1.02957;
      if (Djet>0.1 && Djet<=0.2) scale=0.816706;
      if (Djet>0.2 && Djet<=0.3) scale=1.06065;
      if (Djet>0.3 && Djet<=0.4) scale=0.917763;
      if (Djet>0.4 && Djet<=0.5) scale=1.04721;
      if (Djet>0.5 && Djet<=0.6) scale=0.930451;
      if (Djet>0.6 && Djet<=0.7) scale=0.97422;
      if (Djet>0.7 && Djet<=0.8) scale=1.06606;
      if (Djet>0.8 && Djet<=0.9) scale=1.01095;
      if (Djet>0.9) scale=0.9843;
      }*/
    //d6t tuning mirrored
    if(altnum==2){
      if (Djet<=0.1) scale=0.864039;
      if (Djet>0.1 && Djet<=0.2) scale=0.850918;
      if (Djet>0.2 && Djet<=0.3) scale=0.826864;
      if (Djet>0.3 && Djet<=0.4) scale=0.985373;
      if (Djet>0.4 && Djet<=0.5) scale=0.865197;
      if (Djet>0.5 && Djet<=0.6) scale=0.980772;
      if (Djet>0.6 && Djet<=0.7) scale=0.994331;
      if (Djet>0.7 && Djet<=0.8) scale=1.04053;
      if (Djet>0.8 && Djet<=0.9) scale=1.02748;
      if (Djet>0.9) scale=1.093;
    }   
  }
  //qqZZ
  if(channel==2){
    if(altnum==1){
      if (Djet<=0.1) scale=0.982252;
      if (Djet>0.1 && Djet<=0.2) scale=1.11411;
      if (Djet>0.2 && Djet<=0.3) scale=1.08938;
      if (Djet>0.3 && Djet<=0.4) scale=0.989144;
      if (Djet>0.4 && Djet<=0.5) scale=1.0161;
      if (Djet>0.5 && Djet<=0.6) scale=0.936218;
      if (Djet>0.6 && Djet<=0.7) scale=0.935243;
      if (Djet>0.7 && Djet<=0.8) scale=0.909441;
      if (Djet>0.8 && Djet<=0.9) scale=1.03983;
      if (Djet>0.9) scale=0.840851;
    }
    if(altnum==2){
      if (Djet<=0.1) scale=1./0.982252;
      if (Djet>0.1 && Djet<=0.2) scale=1./1.11411;
      if (Djet>0.2 && Djet<=0.3) scale=1./1.08938;
      if (Djet>0.3 && Djet<=0.4) scale=1./0.989144;
      if (Djet>0.4 && Djet<=0.5) scale=1./1.0161;
      if (Djet>0.5 && Djet<=0.6) scale=1./0.936218;
      if (Djet>0.6 && Djet<=0.7) scale=1./0.935243;
      if (Djet>0.7 && Djet<=0.8) scale=1./0.909441;
      if (Djet>0.8 && Djet<=0.9) scale=1./1.03983;
      if (Djet>0.9) scale=1./0.840851;
    }
  }
	 	 
  if(scale==0) cout<<"ERROR: Djet Template has values <0 or >1."<<endl;
  
  return scale;
  
}
	 	 
TH2F* mirrortemplates(int sampleIndex){
  TFile* alt1,*orig;
  if (sampleIndex==1){
    alt1 = new TFile(destDir + "qqH_vbfMELA_alt.root","OPEN");
    orig = new TFile(destDir + "qqH_vbfMELA.root","OPEN");
  }
  if (sampleIndex==2){
    //cout<<"here?"<<endl;
    alt1 = new TFile(destDir + "qqZZ_vbfMELA_alt.root","OPEN");
    orig = new TFile(destDir + "qqZZ_vbfMELA.root","OPEN");
  }
  TH2F* altfisher = (TH2F*)alt1->Get("H_Djet");
  TH2F* origfisher = (TH2F*)orig->Get("H_Djet");
  TH2F* alt2fisher = new TH2F("H_Djet","H_Djet",int((highMzz-100.)/mBinSize+0.5),100,highMzz,50,0.,1.);

  gStyle->SetPalette(1);

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  altfisher->Draw("colz");

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  origfisher->Draw("colz");

  for(int i=1;i<=origfisher->GetNbinsX();i++){
    for(int j=1;j<=origfisher->GetNbinsY();j++){
      float origval=origfisher->GetBinContent(i,j);
      float alt1val=altfisher->GetBinContent(i,j);
      float scale=origval/alt1val;
      float alt2val=scale*origval;
      alt2fisher->SetBinContent(i,j,alt2val);
    }
  }
  
  TCanvas* c3 = new TCanvas("c3","c3",800,800);
  c3->cd();
  alt2fisher->Draw("colz");

  //Renormalize
  TH1F* tempProj;
  double norm;
  for(int i=1; i<=alt2fisher->GetNbinsX(); i++){
    tempProj = (TH1F*) alt2fisher->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=alt2fisher->GetNbinsY(); j++){
      	alt2fisher->SetBinContent(i,j,alt2fisher->GetBinContent(i,j)/norm);
      }
    }
  }

  TCanvas* c4 = new TCanvas("c4","c4",800,800);
  c4->cd();
  alt2fisher->Draw("colz");
	 	 
  return alt2fisher;
}

void progressbar(int val, int tot){
  int percent=floor(0.01*tot);

  if(val%percent==0 && val!=tot){
    cout<<"[ "<<setw(3)<<val/percent<<"% |";
    for(int k=1;k<val/percent;k++) cout<<"=";
    if(val%percent!=100) cout<<">";
    for(int k=val/percent+1;k<100;k++) cout<<" ";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
  else if(val==tot){
    cout<<"[ 100% |";
    for(int k=0;k<100;k++) cout<<"=";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');        
  }
}

TH2F* makebkg(TH2F *SM1, TH2F *SM10, TH2F *SM25){
  float denomscale = 60.-24.*sqrt(10.);

  TH2F* madebkg = (TH2F*) SM1->Clone("madebkg");
  madebkg->Scale(50.-25.*sqrt(10.));
  madebkg->Add(SM10,20);
  madebkg->Add(SM25,-10.+sqrt(10.));
  //madebkg->Scale(-5.+sqrt(10.));
  //madebkg->Add(SM10,4);
  //madebkg->Add(SM25,-1.+sqrt(10.));
  madebkg->Scale(1/denomscale);
  //madebkg->Scale(1./madebkg->Integral());

  return madebkg;
}