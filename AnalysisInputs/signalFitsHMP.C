/* 
 * Compute shape parameters parametrization as a function of the invariant mass  for signals and write them in a card fragment.
 * Committed for the High Mass paper
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b signalFitsHMP.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "TStyle.h"

/*
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TSystem.h"

#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooPlot.h"
*/

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

// New flags and options
TString appendStr("_C0.2_fitted");
Double_t cprimeVal = 0.2;

bool fixResPars = false;
bool copyToWeb = false;
const char* webDir = "/afs/cern.ch/user/s/scasasso/www/H4l/HighMass/SignalShapesReco/%iTeV/";

//Parameters to choose which samples to examine
//HCP option is useggH=true, useVBF,useVH = false and usedijet,usenondijet=true
bool debug = false;
bool useggH = true;
bool useVBF = false;
bool useVH = false;
bool usedijet = true;
bool usenondijet = true;
bool setParSeed = true;

TFile* massfits;

using namespace RooFit;
using namespace std;

//Declaration
void signalFitsHMP(int channel, int sqrts, ofstream* interpCode=0);
float WidthValue(float mHStarWidth);
float highparameter(TF1 *generatedfit, int gamma);

void signalFitsHMP(){
  gSystem->Exec("mkdir -p sigFigs7TeV_HMP");
  gSystem->Exec("mkdir -p sigFigs8TeV_HMP");

  gSystem->Load("../CreateDatacards/CMSSW_6_1_1/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so");

  if (!usedijet && !usenondijet){
    cout << "Neither dijet tagging category was chosen. Please choose one or both." << endl;
    abort();
  }

  if (!useggH && !useVBF && !useVH){
    cout << "Please choose at least one production method: ggH, VBF, VH."<<endl;
    abort();
  }
  
  ofstream txtInterp;
  string outInterpCodeName = "interpolation.txt";
  txtInterp.open(outInterpCodeName.c_str(),fstream::out);

  signalFitsHMP(1,8,&txtInterp);
  signalFitsHMP(2,8,&txtInterp);
  signalFitsHMP(3,8,&txtInterp);
  signalFitsHMP(1,7,&txtInterp);
  signalFitsHMP(2,7,&txtInterp);
  signalFitsHMP(3,7,&txtInterp);

  txtInterp.close();

  return;
}

//The actual job
void signalFitsHMP(int channel, int sqrts, ofstream* interpCode){

  TString cprimeTS;
  if (cprimeVal==0.2) cprimeTS = "C0.2";
  else cprimeTS = "C1.0";

  TString ssqrts;
  if (sqrts==7) ssqrts = "7TeV";
  else if (sqrts==8) ssqrts = "8TeV";
  else {
    std::cout << "Choose 7 or 8 TeV as center of mass energy" << std::endl;
    exit(1);
  }

  (*interpCode) << "      if (sqrts==" << sqrts << "){" << endl;
  (*interpCode) << "	    if (channel==" << channel << "){" << endl;
  (*interpCode) << "  	      // Fitted on the C" << cprimeVal << endl;

  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;
  Bool_t is8TeV = true;
  if (sqrts==7) is8TeV=false;


  //Pick the correct mass points and paths
  TString filePath;
  int nPointsggH,nPointsVBF,nPointsVH;
  double* massesggH,*massesVBF,*massesVH;

  if (sqrts==7) {
    nPointsggH = nPoints7TeV_p15_HM;
    massesggH  = mHVal7TeV_p15_HM;
    nPointsVBF = nVBFPoints7TeV_HM;
    massesVBF  = mHVBFVal7TeV_HM;
    nPointsVH = nVHPoints7TeV;
    massesVH  = mHVHVal7TeV;
    filePath = filePath7TeV;  
  } else if (sqrts==8) {
    nPointsggH = nPoints8TeV_p15_HM;
    massesggH  = mHVal8TeV_p15_HM;
    nPointsVBF = nVBFPoints8TeV_HM;
    massesVBF  = mHVBFVal8TeV_HM;
    nPointsVH = nVHPoints8TeV;
    massesVH  = mHVHVal8TeV;
    filePath = filePath8TeV;
  }
  else abort();

  int nPoints;
  double masses[200];
  double masses_err[200];
  for (int i=0;i<200;i++){
    masses[i]=-1;
    masses_err[i]=0;
  }
  if (useggH){
    for (int i=0; i<nPointsggH; i++){
      masses[i]=massesggH[i];
    }
    nPoints=nPointsggH;
  }else if(!useggH && useVBF){
    for (int i=0; i<nPointsVBF; i++){
      masses[i]=massesVBF[i];
    }
    nPoints=nPointsVBF;
  }else if(!useggH && !useVBF && useVH){
    for (int i=0; i<nPointsVH; i++){
      masses[i]=massesVH[i];
    }
      nPoints=nPointsVH;
  }


  bool flag=0;
  //If any mass points are in VBF/VH but not in ggH, this will include them
  if (useggH && useVBF){
    for (int i=0; i<nPointsVBF; i++){
      flag=1;
      for (int j=0; j<nPoints; j++){
	if (massesVBF[i]==masses[j]) flag=0;
      }
      if (flag==1){
	masses[nPoints]=massesVBF[i];
	nPoints++;
      }
    }
  }
  if ((useggH || useVBF) && useVH){
    for (int i=0; i<nPointsVH; i++){
      flag=1;
      for (int j=0; j<nPoints; j++){
	if (massesVH[i]==masses[j]) flag=0;
      }
      if (flag==1){
	masses[nPoints]=massesVH[i];
	nPoints++;
      }
    }
  }

  filePath.Append(schannel=="2e2mu"?"2mu2e":schannel);

  //Prepare to store all the shape parameters for the mass points	
  const int arraySize=200;
  assert(arraySize >= nPoints);

  Double_t a_meanCB[arraySize];    Double_t a_meanCB_err[arraySize]; 
  Double_t a_sigmaCB[arraySize];   Double_t a_sigmaCB_err[arraySize];
  Double_t a_alphaCB_1[arraySize]; Double_t a_alphaCB_1_err[arraySize];
  Double_t a_nCB_1[arraySize];     Double_t a_nCB_1_err[arraySize];
  Double_t a_Gamma[arraySize];     Double_t a_Gamma_err[arraySize];
  Double_t a_alphaCB_2[arraySize]; Double_t a_alphaCB_2_err[arraySize];
  Double_t a_nCB_2[arraySize];     Double_t a_nCB_2_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];

  char outfile[192];
  sprintf(outfile,"sigFigs%iTeV_HMP",sqrts);
                
  //Loop over the mass points
  for (int i = 0; i < nPoints; i++){

    bool flagVBF==0;
    bool flagVH==0;

    if (debug && masses[i]!=126.) continue;
		
    //Open input file with shapes and retrieve the tree
    char tmp_finalInPathggH[200],tmp_finalInPathVBF[200],tmp_finalInPathZH[200],tmp_finalInPathWH[200],tmp_finalInPathttH[200];
    sprintf(tmp_finalInPathggH,"/HZZ4lTree_powheg15H%i.root",masses[i]);
    sprintf(tmp_finalInPathVBF,"/HZZ4lTree_powheg15VBFH%i.root",masses[i]);
    sprintf(tmp_finalInPathZH,"/HZZ4lTree_ZH%i.root",masses[i]);
    sprintf(tmp_finalInPathWH,"/HZZ4lTree_WH%i.root",masses[i]);
    sprintf(tmp_finalInPathttH,"/HZZ4lTree_ttH%i.root",masses[i]);
    TString finalInPathggH = filePath + tmp_finalInPathggH;
    TString finalInPathVBF = filePath + tmp_finalInPathVBF;
    TString finalInPathZH = filePath + tmp_finalInPathZH;
    TString finalInPathWH = filePath + tmp_finalInPathWH;
    TString finalInPathttH = filePath + tmp_finalInPathttH;

    TChain *f = new TChain("SelectedTree");
    if (useggH){
      if (i<nPointsggH) f->Add(finalInPathggH);
      else{
	cout<<"No ggH sample at this mass point."<<endl;
      }
    }
    else if (!useggH && useVBF){
      if (i<nPointsVBF) f->Add(finalInPathVBF);
      else{
	cout<<"No qqH sample at this mass point."<<endl;
      }
    }
    else if (!useggH && !useVBF && useVH){
      if (i<nPointsVH){
	f->Add(finalInPathZH);
	f->Add(finalInPathWH);
	f->Add(finalInPathttH);
      }
      else{
	cout<<"No VH samples at this mass point."<<endl;
      }
    }    
    if (useggH && useVBF){
      for (int j=0;j<nPointsVBF;j++){
	if (massesVBF[j]=masses[i]) flagVBF=1;
      }
      if (flagVBF==1) f->Add(finalInPathVBF);
      if (flagVBF==0) cout<<"No qqH sample at this mass point."<<endl;
    }
    if ((useggH || useVBF) && useVH){
      for (int j=0;j<nPointsVH;j++){
	if (massesVH[j]=masses[i]) flagVH=1;
      }
      if (flagVH==1){
	f->Add(finalInPathZH);
	f->Add(finalInPathWH);
	f->Add(finalInPathttH);
      }
      if (flagVH==0) cout<<"No VH samples at this mass point."<<endl;
    }

    double valueWidth = WidthValue(masses[i]);
    Double_t gammaVal = valueWidth;
    double windowVal = max(valueWidth,1.);
    double lowside = 100.;
    if(masses[i] > 300) lowside = 200.;
    double low_M = max( (masses[i] - 15.*windowVal), lowside) ;
    double high_M = min( (masses[i] + 10.*windowVal), 1400.);


    // FIXME: as soon as fit ranges can be defined in a continuous way, replace this twofold ranges definition
    if(masses[i] > 399.){

      if (cprimeVal>0.9){
	// For the C'^2 = 1.0 case
	if (useggH){
	  low_M = max( (masses[i] - 2.*windowVal), 250.) ;
	  high_M = min( (masses[i] + 2.*windowVal), 1450.);
	  // high_M = min( (masses[i] + 2.*windowVal), 1300.);
	}
	else if (!useggH && useVBF){
	  low_M = max( (masses[i] - 2.*windowVal), 250.) ;
	  high_M = min( (masses[i] + 2.*windowVal), 1500.);
	}
      }
      else{
	// For the C'^2 = 0.2 case
	if (masses[i]<401.){
	  low_M = 360;
	  high_M = 430; 
	}
	else if (masses[i]<451.){
	  low_M = 400;
	  high_M = 490; 	
	}
	else if (masses[i]<501.){
	  low_M = 440;
	  high_M = 540; 	
	}
	else if (masses[i]<551.){
	  low_M = 485;
	  high_M = 600; 	
	}
	else if (masses[i]<601.){
	  low_M = 520.;
	  high_M = 660.; 	
	}
	else if (masses[i]<651.){
	  low_M = 570.;
	  high_M = 710.; 	
	}
	else if (masses[i]<701.){
	  low_M = 610.;
	  high_M = 765.; 	
	}
	else if (masses[i]<801.){
	  low_M = 690.;
	  high_M = 870.; 	
	}
	else if (masses[i]<901.){
	  low_M = 740.;
	  high_M = 980.; 	
	}
	else if (masses[i]<1001.){
	  low_M = 820.;
	  high_M = 1100.; 	
	}
	else {
	  low_M = 300;
	  high_M = 1600; 	
	}
      }
      
    }

    double low_M_plot = low_M;
    double high_M_plot = high_M;
    

    cout << "lowM = " << low_M << ", highM = " << high_M << endl;

    //Set the observable and get the RooDataSomething
    RooRealVar ZZMass("ZZMass","ZZMass",low_M_plot,high_M_plot);
    ZZMass.setRange("plotRange",low_M_plot,high_M_plot);
    ZZMass.setRange("fitRange",low_M,high_M);
    RooRealVar MC_weight("MC_weight","MC_weight",0.,10.);
    RooRealVar NJets("NJets","NJets",0.,100.);
    RooRealVar genProcessId("genProcessId","genProcessId",0.,150.);

    ZZMass.setBins(50);

    RooDataSet* set;

    if (!flagVH){
      if (usedijet && !usenondijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets), "NJets>1", "MC_weight");
      }
      else if (usenondijet && !usedijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets), "NJets<2", "MC_weight");
      } else{
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight), "", "MC_weight");
      }
    }
    if (flagVH){
      if (usedijet && !usenondijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets,genProcessId), "NJets>1 && (genProcessId==24 || genProcessId==26 || genProcessId==121 || genProcessId==122 || genProcessId== 10011 )", "MC_weight");
      }
      else if (usenondijet && !usedijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets,genProcessId), "NJets<2 && (genProcessId==24 || genProcessId==26 || genProcessId==121 || genProcessId==122 || genProcessId== 10011)", "MC_weight");
      } else{
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets,genProcessId), "(genProcessId==24 || genProcessId==26 || genProcessId==121 || genProcessId==122 || genProcessId== 10011)", "MC_weight");
      }
    }

    RooDataHist *hist = (RooDataHist*)set->binnedClone("datahist","datahist"); 

    //Theoretical signal model  
    RooRealVar MHStar("MHStar","MHStar",masses[i]);
    RooRealVar Gamma_TOT("Gamma_TOT","Gamma_TOT",gammaVal,0.,700.);
    Gamma_TOT.setConstant(true);
    RooRealVar one("one","one",1.0);
    one.setConstant(kTRUE);

    RooRealVar CSquared("CSquared","C'^{2}",cprimeVal);
    RooRealVar BRnew("BRnew","BR_{new}",0.);
    RooRealVar IntStr("IntStr","#IntStr",1.);

    RooAbsPdf* SignalTheor;
    if (useggH) {
      SignalTheor = new RooBWHighMassGGH("model","model",ZZMass,MHStar,CSquared,BRnew,IntStr,is8TeV);
    }
    else if (useVBF) {
      SignalTheor = new RooCPSHighMassVBF("model","model",ZZMass,MHStar,CSquared,BRnew,IntStr,is8TeV);
    }
    else {
      std::cout<<"You have to set at least 'useggH' or 'useVBF' to true"<<std::endl;
      exit(1);
    }

    //Experimental resolution
    RooRealVar meanCB("meanCB","meanCB",0.,-1.0,1.0);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.,5.);
    RooRealVar sigmaCB_high("sigmaCB_high","sigmaCB_high",5,3.,8.);
    RooRealVar alphaCB_1("alphaCB_1","alphaCB_1",1.,0.4,2.);
    RooRealVar nCB_1("nCB_1","nCB_1",5.,0.,12.);
    nCB_1.setConstant(kTRUE);
    RooRealVar alphaCB_2("alphaCB_2","alphaCB_2",1.,0.4,2.);
    RooRealVar nCB_2("nCB_2","nCB_2",20.,0.,12.);
    nCB_2.setConstant(kTRUE);

    //Initialize to decent values
    float m = masses[i];


    // Fixing the values to the one extracted from the C'^2 = 0.2 case (extracted from ggH case)
    if (useggH && setParSeed){ 
      if (sqrts==8){
	if (channel==1){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(4.53024+-0.00426436*m+2.8606e-06*m*m+1.18084e-08*m*m*m+6.76142e-12*m*m*m*m+-1.65141e-14*m*m*m*m*m);
	  meanCB.setVal(0.880437+-0.00175511*m+-5.3235e-07*m*m+1.54134e-09*m*m*m+1.45748e-12*m*m*m*m+-1.01434e-15*m*m*m*m*m);
	  alphaCB_1.setVal(1.21602+-0.000110275*m+-6.71223e-07*m*m+-2.29648e-10*m*m*m+3.78332e-13*m*m*m*m+5.47917e-16*m*m*m*m*m);
	  alphaCB_2.setVal(2.86583+-0.00436325*m+-1.34961e-06*m*m+4.90553e-09*m*m*m+4.51724e-12*m*m*m*m+-6.01361e-15*m*m*m*m*m);
        }
      }
      if (sqrts==8){
	if (channel==2){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(10.9323+-0.0375246*m+5.03703e-05*m*m+3.11412e-08*m*m*m+-8.38849e-11*m*m*m*m+3.44008e-14*m*m*m*m*m);
	  meanCB.setVal(-12.3553+0.0674625*m+-9.97082e-05*m*m+-1.84552e-08*m*m*m+1.35133e-10*m*m*m*m+-7.20926e-14*m*m*m*m*m);
	  alphaCB_1.setVal(8.39501+-0.0369883*m+4.75209e-05*m*m+2.54207e-08*m*m*m+-8.60408e-11*m*m*m*m+4.26562e-14*m*m*m*m*m);
	  alphaCB_2.setVal(2.08321+-0.00079577*m+-7.3537e-06*m*m+1.32439e-08*m*m*m+-6.05768e-12*m*m*m*m+-4.60935e-16*m*m*m*m*m);
        }
      }
      if (sqrts==8){
	if (channel==3){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(-16.93+0.112952*m+-0.000152425*m*m+-4.66275e-08*m*m*m+2.18637e-10*m*m*m*m+-1.08984e-13*m*m*m*m*m);
	  meanCB.setVal(-19.4116+0.102812*m+-0.000152008*m*m+-3.40849e-08*m*m*m+2.27168e-10*m*m*m*m+-1.25426e-13*m*m*m*m*m);
	  alphaCB_1.setVal(1.67181+0.000115226*m+-3.93008e-06*m*m+-9.42681e-10*m*m*m+8.39824e-12*m*m*m*m+-4.77952e-15*m*m*m*m*m);
	  alphaCB_2.setVal(2.94976+-0.00277342*m+-4.74973e-06*m*m+1.96315e-09*m*m*m+1.10653e-11*m*m*m*m+-8.02456e-15*m*m*m*m*m);
        }
      }
      if (sqrts==7){
	if (channel==1){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(31.6905+-0.1627*m+0.000268251*m*m+4.893e-08*m*m*m+-4.0735e-10*m*m*m*m+2.29134e-13*m*m*m*m*m);
	  meanCB.setVal(-16.4877+0.0891798*m+-0.000137443*m*m+-2.34985e-08*m*m*m+2.067e-10*m*m*m*m+-1.18404e-13*m*m*m*m*m);
	  alphaCB_1.setVal(5.81384+-0.0242094*m+3.30583e-05*m*m+1.18931e-08*m*m*m+-5.59469e-11*m*m*m*m+2.99477e-14*m*m*m*m*m);
	  alphaCB_2.setVal(6.27062+-0.0252813*m+3.26626e-05*m*m+1.32288e-08*m*m*m+-5.45376e-11*m*m*m*m+2.81544e-14*m*m*m*m*m);
        }
      }
      if (sqrts==7){
	if (channel==2){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(-12.7371+0.084158*m+-0.000114061*m*m+-2.97707e-08*m*m*m+1.6838e-10*m*m*m*m+-9.08215e-14*m*m*m*m*m);
	  meanCB.setVal(18.1302+-0.0980912*m+0.000143725*m*m+3.71438e-08*m*m*m+-2.13802e-10*m*m*m*m+1.1303e-13*m*m*m*m*m);
	  alphaCB_1.setVal(3.02238+-0.00968557*m+1.28765e-05*m*m+3.51457e-09*m*m*m+-1.77337e-11*m*m*m*m+9.03178e-15*m*m*m*m*m);
	  alphaCB_2.setVal(3.31674+-0.0100594*m+1.23872e-05*m*m+3.5756e-09*m*m*m+-1.72705e-11*m*m*m*m+9.00383e-15*m*m*m*m*m);
        }
      }
      if (sqrts==7){
	if (channel==3){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(7.85707+-0.0329166*m+7.185e-05*m*m+4.24507e-09*m*m*m+-1.09577e-10*m*m*m*m+6.50213e-14*m*m*m*m*m);
	  meanCB.setVal(1.26046+-0.00388926*m+-3.98554e-07*m*m+4.76425e-09*m*m*m+3.86687e-12*m*m*m*m+-5.05511e-15*m*m*m*m*m);
	  alphaCB_1.setVal(0.521945+0.00160607*m+3.12609e-07*m*m+-3.14525e-09*m*m*m+-2.06383e-12*m*m*m*m+3.4318e-15*m*m*m*m*m);
	  alphaCB_2.setVal(1.85608+-0.00160577*m+-4.84323e-07*m*m+3.25779e-10*m*m*m+1.04032e-12*m*m*m*m+-5.80937e-16*m*m*m*m*m);
	}
      }
    }

    // Fixing the values to the one extracted from the C'^2 = 0.2 case (extracted from VBF case)
    if (!useggH && useVBF && setParSeed){ 
      if (sqrts==8){
	    if (channel==1){
  	      // Fitted on the C0.2
              sigmaCB_high.setVal(12.5956+-0.0470057*m+6.83473e-05*m*m+2.51305e-08*m*m*m+-1.06192e-10*m*m*m*m+5.28732e-14*m*m*m*m*m);
              meanCB.setVal(-3.78949+0.010069*m+2.0648e-06*m*m+-1.19493e-08*m*m*m+-1.04701e-11*m*m*m*m+1.43896e-14*m*m*m*m*m);
              alphaCB_1.setVal(-0.331026+0.00794445*m+-1.38703e-05*m*m+2.74076e-09*m*m*m+1.2858e-11*m*m*m*m+-8.7033e-15*m*m*m*m*m);
              alphaCB_2.setVal(-0.948574+0.0103706*m+-1.39029e-05*m*m+-1.13422e-09*m*m*m+9.85803e-12*m*m*m*m+-3.33141e-15*m*m*m*m*m);
        }
      }
      if (sqrts==8){
	    if (channel==2){
  	      // Fitted on the C0.2
              sigmaCB_high.setVal(-10.0184+0.0697895*m+-9.60631e-05*m*m+-2.45528e-08*m*m*m+1.53632e-10*m*m*m*m+-8.91256e-14*m*m*m*m*m);
              meanCB.setVal(-2.1427+0.0052267*m+1.29187e-06*m*m+-5.58547e-09*m*m*m+-4.95629e-12*m*m*m*m+6.37905e-15*m*m*m*m*m);
              alphaCB_1.setVal(-3.18761+0.0244051*m+-3.81919e-05*m*m+-1.22216e-08*m*m*m+6.95187e-11*m*m*m*m+-3.99139e-14*m*m*m*m*m);
              alphaCB_2.setVal(-3.11972+0.0246487*m+-3.82539e-05*m*m+-1.30719e-08*m*m*m+6.87467e-11*m*m*m*m+-3.85156e-14*m*m*m*m*m);
        }
      }
      if (sqrts==8){
	    if (channel==3){
  	      // Fitted on the C0.2
              sigmaCB_high.setVal(20.4808+-0.0935922*m+0.000149709*m*m+3.27243e-08*m*m*m+-2.16736e-10*m*m*m*m+1.13581e-13*m*m*m*m*m);
              meanCB.setVal(0.135166+-0.000940303*m+-2.9491e-09*m*m+1.79791e-09*m*m*m+1.70383e-12*m*m*m*m+-1.96563e-15*m*m*m*m*m);
              alphaCB_1.setVal(3.26293+-0.0113857*m+1.59068e-05*m*m+3.53184e-09*m*m*m+-2.27015e-11*m*m*m*m+1.18839e-14*m*m*m*m*m);
              alphaCB_2.setVal(7.20668+-0.0307051*m+4.1712e-05*m*m+1.16758e-08*m*m*m+-6.05437e-11*m*m*m*m+3.10839e-14*m*m*m*m*m);
        }
      }
      if (sqrts==7){
	    if (channel==1){
  	      // Fitted on the C0.2
              sigmaCB_high.setVal(29.7087+-0.14648*m+0.000232374*m*m+5.19291e-08*m*m*m+-3.46633e-10*m*m*m*m+1.85597e-13*m*m*m*m*m);
              meanCB.setVal(-2.00322+0.00524247*m+8.5533e-07*m*m+-6.15464e-09*m*m*m+-4.92105e-12*m*m*m*m+7.54849e-15*m*m*m*m*m);
              alphaCB_1.setVal(2.19357+-0.00549671*m+6.63663e-06*m*m+4.56663e-09*m*m*m+-1.58807e-11*m*m*m*m+8.514e-15*m*m*m*m*m);
              alphaCB_2.setVal(7.65399+-0.0323573*m+4.23417e-05*m*m+1.59591e-08*m*m*m+-6.72994e-11*m*m*m*m+3.41317e-14*m*m*m*m*m);
        }
      }
      if (sqrts==7){
	    if (channel==2){
  	      // Fitted on the C0.2
              sigmaCB_high.setVal(14.2734+-0.0639824*m+0.000112014*m*m+1.84252e-08*m*m*m+-1.8111e-10*m*m*m*m+1.08221e-13*m*m*m*m*m);
              meanCB.setVal(-1.20102+0.00270735*m+1.05305e-06*m*m+-2.86439e-09*m*m*m+-3.06745e-12*m*m*m*m+3.48754e-15*m*m*m*m*m);
              alphaCB_1.setVal(0.334235+0.00567805*m+-1.29738e-05*m*m+5.36553e-09*m*m*m+1.11176e-11*m*m*m*m+-8.79236e-15*m*m*m*m*m);
              alphaCB_2.setVal(0.67474+0.00505365*m+-1.33184e-05*m*m+5.83472e-09*m*m*m+1.17403e-11*m*m*m*m+-9.24834e-15*m*m*m*m*m);
        }
      }
      if (sqrts==7){
	    if (channel==3){
  	      // Fitted on the C0.2
              sigmaCB_high.setVal(2.0798+-0.00293226*m+3.09867e-05*m*m+-1.77069e-08*m*m*m+-3.10432e-11*m*m*m*m+2.46706e-14*m*m*m*m*m);
              meanCB.setVal(-9.8286+0.048186*m+-6.31769e-05*m*m+-2.39942e-08*m*m*m+9.55702e-11*m*m*m*m+-4.61232e-14*m*m*m*m*m);
              alphaCB_1.setVal(4.31033+-0.0168256*m+2.65514e-05*m*m+-5.65718e-09*m*m*m+-1.92652e-11*m*m*m*m+1.1468e-14*m*m*m*m*m);
              alphaCB_2.setVal(4.39198+-0.0166661*m+2.61775e-05*m*m+-6.16957e-09*m*m*m+-1.92985e-11*m*m*m*m+1.20878e-14*m*m*m*m*m);
        }
      }
    }


    if (fixResPars){
      sigmaCB_high.setConstant(kTRUE);
      meanCB.setConstant(kTRUE);	     
      alphaCB_1.setConstant(kTRUE);   
      alphaCB_2.setConstant(kTRUE);   
      nCB_1.setConstant(kTRUE);   
      nCB_2.setConstant(kTRUE);         
    }


    RooDoubleCB massRes("massRes","Double Crystal Ball",ZZMass,meanCB,sigmaCB,alphaCB_1,nCB_1,alphaCB_2,nCB_2);
    RooDoubleCB massResH("massResH","DCB Highmass",ZZMass,meanCB,sigmaCB_high,alphaCB_1,nCB_1,alphaCB_2,nCB_2);

    //Convolute theoretical shape and resolution
    RooFFTConvPdf *sigPDF;
    if(masses[i] < 399.) {
      sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,*SignalTheor,massRes);
    }
    else{
      sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,*SignalTheor,massResH);
    }
    sigPDF->setBufferFraction(0.2);

    RooPlot *xplot = ZZMass.frame();
    RooPlot *xplot_z = ZZMass.frame();
    TCanvas *canv = new TCanvas("canv","canv",1200,800);
    TCanvas *canv_z = new TCanvas("canv_z","canv_z",1200,800);

    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%i_%iTeV_",masses[i],sqrts);
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel;
    TString rootTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel;
    if (useggH){
      plotFileTitle+="_ggH";
      rootTitle+="_ggH";
    }
    if (useVBF){
      plotFileTitle+="_VBF";
      rootTitle+="_VBF";
    }
    if (useVH){
      plotFileTitle+="_VH";
      rootTitle+="_VH";
    }
    if (usedijet && !usenondijet){
      plotFileTitle += "_1";
      rootTitle += "_1";
    }
    if (usenondijet && !usedijet){
      plotFileTitle += "_0";
      rootTitle += "_0";
    }

    double mass,mean,sigma,a1,n1,a2,n2,gamma;
    TH1F* parameters;

    //Fit the shape
    massfits = new TFile(rootTitle + ".root","RECREATE");
    
    RooFitResult *fitRes = sigPDF->fitTo(*hist,Save(1), SumW2Error(kTRUE), Range("fitRange"));
    // RooFitResult *fitRes = sigPDF->chi2FitTo(*hist,Save(1), SumW2Error(kTRUE), Range("fitRange"), Strategy(2));

//     a_fitEDM[i] = fitRes->edm();
//     a_fitCovQual[i] = fitRes->covQual();
//     a_fitStatus[i] = fitRes->status();

    mass = masses[i];
    mean = meanCB.getVal();
    if (mass > 399.) sigma = sigmaCB_high.getVal();
    else sigma = sigmaCB.getVal();
    a1=alphaCB_1.getVal();
    n1=nCB_1.getVal();
    a2=alphaCB_2.getVal();
    n2=nCB_2.getVal();
    gamma=Gamma_TOT.getVal();

    a_meanCB[i]  = mean;
    a_sigmaCB[i] = sigma;
    a_alphaCB_1[i]  = a1;
    a_nCB_1[i]     = n1;
    a_alphaCB_2[i]  = a2;
    a_nCB_2[i]     = n2;
    a_Gamma[i] = gamma;
    
    a_meanCB_err[i]  = meanCB.getError();
    if (masses[i] > 399.) a_sigmaCB_err[i] = sigmaCB_high.getError();
    else a_sigmaCB_err[i] = sigmaCB.getError();
    a_alphaCB_1_err[i]  = alphaCB_1.getError();
    a_nCB_1_err[i]     = nCB_1.getError();
    a_alphaCB_2_err[i]  = alphaCB_2.getError();
    a_nCB_2_err[i]     = 0;
    if(masses[i] > 399.) a_Gamma_err[i] = Gamma_TOT.getError();
    else a_Gamma_err[i] = 0.;

    //Plot in the figures directory
    hist->plotOn(xplot);
    sigPDF->plotOn(xplot,Range("plotRange"));
    canv->cd();
    xplot->Draw();

    hist->plotOn(xplot_z);
    sigPDF->plotOn(xplot_z,Range("fitRange"));
    canv_z->cd();
    xplot_z->Draw();

    TString plotFileTitleTS(plotFileTitle.c_str());
    TString plotgif = plotFileTitleTS + ""+appendStr+".png";
    TString plotpdf = plotFileTitleTS + ""+appendStr+".pdf";
    char baseWWW[192];
    sprintf(baseWWW,webDir,sqrts);
    TString TS_baseWWW(baseWWW);
    TS_baseWWW+=cprimeTS;

    canv->SaveAs(plotgif);
    canv->SaveAs(plotpdf);
    gSystem->Exec("cp "+plotgif+" "+TS_baseWWW);
    plotgif.ReplaceAll(".png","_z.png");
    plotpdf.ReplaceAll(".pdf","_z.pdf");
    // gSystem->Exec("cp "+plotgif+" "+TS_baseWWW);
    canv_z->SaveAs(plotgif);
    canv_z->SaveAs(plotpdf);

    massfits->cd();
    set->Write("MassData");
    parameters = new TH1F("","",8,0,8);
    parameters->Fill(0,mass);
    parameters->SetBinError(1,0);
    parameters->Fill(1,mean);
    parameters->SetBinError(2,a_meanCB_err[i]);
    parameters->Fill(2,sigma);
    parameters->SetBinError(3,a_sigmaCB_err[i]);
    parameters->Fill(3,a1);
    parameters->SetBinError(4,a_alphaCB_1_err[i]);
    parameters->Fill(4,n1);
    parameters->SetBinError(5,a_nCB_1_err[i]);
    parameters->Fill(5,a2);
    parameters->SetBinError(6,a_alphaCB_2_err[i]);
    parameters->Fill(6,n2);
    parameters->SetBinError(7,a_nCB_2_err[i]);
    parameters->Fill(7,gamma);
    parameters->SetBinError(8,a_Gamma_err[i]);
    parameters->Write("Parameters");
    massfits->Close();
    
  }


  if (debug) return;

  TGraph* gr_meanCB  = new TGraph(nPoints, masses, a_meanCB);
  TGraph* gr_sigmaCB = new TGraph(nPoints, masses, a_sigmaCB);
  TGraph* gr_alphaCB_1 = new TGraph(nPoints, masses, a_alphaCB_1);
  TGraph* gr_nCB_1     = new TGraph(nPoints, masses, a_nCB_1);
  TGraph* gr_alphaCB_2 = new TGraph(nPoints, masses, a_alphaCB_2);
  TGraph* gr_nCB_2     = new TGraph(nPoints, masses, a_nCB_2);
  TGraph* gr_Gamma   = new TGraph(nPoints, masses, a_Gamma);

//   TGraphErrors* gr_meanCB  = new TGraphErrors(nPoints, masses, a_meanCB, masses_err, a_meanCB_err);
//   TGraphErrors* gr_sigmaCB = new TGraphErrors(nPoints, masses, a_sigmaCB, masses_err, a_sigmaCB_err);
//   TGraphErrors* gr_alphaCB_1 = new TGraphErrors(nPoints, masses, a_alphaCB_1, masses_err, a_alphaCB_1_err);
//   TGraphErrors* gr_nCB_1     = new TGraphErrors(nPoints, masses, a_nCB_1, masses_err, a_nCB_1_err);
//   TGraphErrors* gr_alphaCB_2 = new TGraphErrors(nPoints, masses, a_alphaCB_2, masses_err, a_alphaCB_2_err);
//   TGraphErrors* gr_nCB_2     = new TGraphErrors(nPoints, masses, a_nCB_2, masses_err, a_nCB_2_err);
//   TGraphErrors* gr_Gamma   = new TGraphErrors(nPoints, masses, a_Gamma, masses_err, a_Gamma_err);

  gr_meanCB ->SetMarkerStyle(20);
  gr_sigmaCB->SetMarkerStyle(20);
  gr_alphaCB_1->SetMarkerStyle(20);
  gr_nCB_1->SetMarkerStyle(20);
  gr_alphaCB_2->SetMarkerStyle(20);
  gr_nCB_2->SetMarkerStyle(20);
  gr_Gamma->SetMarkerStyle(20);  
  
  // TF1 *paramfit = new TF1("paramfit","(x<400)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>=400)*(([0]+([1]-[6])*400+([2]-[7])*pow(400,2)+([3]-[8])*pow(400,3)+[4]*pow(400,4)+[5]*pow(400,5))+[6]*x+[7]*x*x+[8]*x*x*x)",115,1000);
  // TF1 *paramfit = new TF1("paramfit","(([0]+([1]-[6])*400+([2]-[7])*pow(400,2)+([3]-[8])*pow(400,3)+[4]*pow(400,4)+[5]*pow(400,5))+[6]*x+[7]*x*x+[8]*x*x*x)",400,1000);
  TF1 *paramfit = new TF1("paramfit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x",400,1000);
  // TF1 *linearfit = new TF1("linearfit","[0]+[1]*x",400,1000);
  TF1 *gammafit = new TF1("gammafit","(x<400)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>=400)*(([0]+([6]-[2])*pow(400,2)+2*([7]-[3])*pow(400,3)-3*[4]*pow(400,4)-4*[5]*pow(400,5))+([1]+2*([2]-[6])*400+3*([3]-[7])*pow(400,2)+4*[4]*pow(400,3)+5*[5]*pow(400,4))*x+[6]*x*x+[7]*x*x*x)",115,1000);

  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the mean of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_meanCB->Fit("paramfit","MER");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the sigma of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_sigmaCB->Fit("paramfit","MER");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the alpha1 (R) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_alphaCB_1->Fit("paramfit","MER");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the N1 (R) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_nCB_1->Fit("paramfit","MER");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the alpha2 (L) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_alphaCB_2->Fit("paramfit","MER");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the N2 (L) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_nCB_2->Fit("pol0");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the gamma"<<endl;
  cout<<"#################################"<<endl;
  gr_Gamma->Fit("gammafit");

  TF1 *fit_meanCB  = gr_meanCB->GetListOfFunctions()->First();
  TF1 *fit_sigmaCB = gr_sigmaCB->GetListOfFunctions()->First();
  TF1 *fit_alphaCB_1 = gr_alphaCB_1->GetListOfFunctions()->First();
  TF1 *fit_nCB_1     = gr_nCB_1->GetListOfFunctions()->First();
  TF1 *fit_alphaCB_2 = gr_alphaCB_2->GetListOfFunctions()->First();
  TF1 *fit_nCB_2     = gr_nCB_2->GetListOfFunctions()->First();
  TF1 *fit_Gamma   = gr_Gamma->GetListOfFunctions()->First();

  
  (*interpCode) << "              sigmaCB_high.setVal("<<fit_sigmaCB->GetParameter(0)<<"+"<<fit_sigmaCB->GetParameter(1)<<"*m+"<<fit_sigmaCB->GetParameter(2)<<"*m*m+"<<fit_sigmaCB->GetParameter(3)<<"*m*m*m+"<<fit_sigmaCB->GetParameter(4)<<"*m*m*m*m+"<<fit_sigmaCB->GetParameter(5)<<"*m*m*m*m*m);"<<endl;
  (*interpCode) << "              meanCB.setVal("<<fit_meanCB->GetParameter(0)<<"+"<<fit_meanCB->GetParameter(1)<<"*m+"<<fit_meanCB->GetParameter(2)<<"*m*m+"<<fit_meanCB->GetParameter(3)<<"*m*m*m+"<<fit_meanCB->GetParameter(4)<<"*m*m*m*m+"<<fit_meanCB->GetParameter(5)<<"*m*m*m*m*m)"<<";"<<endl;
  (*interpCode) << "              alphaCB_1.setVal("<<fit_alphaCB_1->GetParameter(0)<<"+"<<fit_alphaCB_1->GetParameter(1)<<"*m+"<<fit_alphaCB_1->GetParameter(2)<<"*m*m+"<<fit_alphaCB_1->GetParameter(3)<<"*m*m*m+"<<fit_alphaCB_1->GetParameter(4)<<"*m*m*m*m+"<<fit_alphaCB_1->GetParameter(5)<<"*m*m*m*m*m)"<<";"<<endl;
  (*interpCode) << "              alphaCB_2.setVal("<<fit_alphaCB_2->GetParameter(0)<<"+"<<fit_alphaCB_2->GetParameter(1)<<"*m+"<<fit_alphaCB_2->GetParameter(2)<<"*m*m+"<<fit_alphaCB_2->GetParameter(3)<<"*m*m*m+"<<fit_alphaCB_2->GetParameter(4)<<"*m*m*m*m+"<<fit_alphaCB_2->GetParameter(5)<<"*m*m*m*m*m)"<<";"<<endl;

  (*interpCode) << "        }" << endl;
  (*interpCode) << "      }" << endl;

  gr_meanCB->SetTitle("Mean value of the DCB function");
  gr_sigmaCB->SetTitle("Sigma of the DCB function");
  gr_alphaCB_1->SetTitle("Alpha parameter of the R leg of DCB function");
  gr_nCB_1->SetTitle("n parameter of the R leg of DCB function");
  gr_alphaCB_2->SetTitle("Alpha parameter of the L leg of DCB function");
  gr_nCB_2->SetTitle("n parameter of the L leg of DCB function");
  gr_Gamma->SetTitle("#Gamma of the BW function");

//   gr_meanCB->SetTitle("");
//   gr_sigmaCB->SetTitle("");
//   gr_alphaCB_1->SetTitle("");
//   gr_nCB_1->SetTitle("");
//   gr_alphaCB_2->SetTitle("");
//   gr_nCB_2->SetTitle("");
//   gr_Gamma->SetTitle("");

//   TCanvas *canv2 = new TCanvas("canv2","canv2",1600,800);
//   canv2->Divide(4,2);

//   canv2->cd(1); gr_meanCB->Draw("A*");  fit_meanCB->Draw("SAME");
//   canv2->cd(2); gr_sigmaCB->Draw("A*"); fit_sigmaCB->Draw("SAME");
//   canv2->cd(3); gr_alphaCB_1->Draw("A*"); fit_alphaCB_1->Draw("SAME");
//   canv2->cd(4); gr_nCB_1->Draw("A*");     fit_nCB_1->Draw("SAME");
//   canv2->cd(5); gr_alphaCB_2->Draw("A*"); fit_alphaCB_2->Draw("SAME");
//   canv2->cd(6); gr_nCB_2->Draw("A*");     fit_nCB_2->Draw("SAME");
//   canv2->cd(7); gr_Gamma->Draw("A*");   fit_Gamma->Draw("SAME");

  TCanvas *canv2 = new TCanvas("canv2_"+TString(schannel)+"_"+ssqrts,"canv2_"+TString(schannel)+"_"+ssqrts,1600,400);
  canv2->Divide(4);

  canv2->cd(1); gr_meanCB->Draw("AP");  fit_meanCB->Draw("SAME");
  canv2->cd(2); gr_sigmaCB->Draw("AP"); fit_sigmaCB->Draw("SAME");
  canv2->cd(3); gr_alphaCB_1->Draw("AP"); fit_alphaCB_1->Draw("SAME");
  canv2->cd(4); gr_alphaCB_2->Draw("AP"); fit_alphaCB_2->Draw("SAME");


  gr_meanCB->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_sigmaCB->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_alphaCB_1->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_alphaCB_2->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_nCB_1->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_nCB_2->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_Gamma->GetXaxis()->SetTitle("m_{H} (GeV)");

  gr_meanCB->GetYaxis()->SetRangeUser(-3,3);
  gr_sigmaCB->GetYaxis()->SetRangeUser(2,12);
  gr_alphaCB_1->GetYaxis()->SetRangeUser(0,3);
  gr_alphaCB_2->GetYaxis()->SetRangeUser(0,3);

  string tmp_paramPlotFileTitle;
  tmp_paramPlotFileTitle.insert(0,outfile);
  tmp_paramPlotFileTitle += "/fitParam_";
  char tmp2_paramPlotFileTitle[200];
  sprintf(tmp2_paramPlotFileTitle,"%iTeV_",sqrts);
  string paramPlotFileTitle = tmp_paramPlotFileTitle + tmp2_paramPlotFileTitle + schannel;
  if (useggH) paramPlotFileTitle+="_ggH";
  if (useVBF) paramPlotFileTitle+="_VBF";
  if (useVH) paramPlotFileTitle+="_VH";
  string paramgif=paramPlotFileTitle;
  string parampdf=paramPlotFileTitle;
  TString paramgifTS(paramPlotFileTitle.c_str());
  TString parampdfTS(paramPlotFileTitle.c_str());
  if (usedijet && !usenondijet){
    paramgifTS+="_0"+appendStr+".png";
    parampdfTS+="_0"+appendStr+".pdf";
  }else if (usenondijet && !usedijet){
    paramgifTS+="_1"+appendStr+".png";
    parampdfTS+="_1"+appendStr+".pdf";
  }else if (usenondijet && usedijet){
    paramgifTS+="deriv5"+appendStr+".png";
    parampdfTS+="deriv5"+appendStr+".pdf";
  }
  canv2->SaveAs(paramgifTS);
  canv2->SaveAs(parampdfTS);
  if (copyToWeb) gSystem->Exec("cp "+paramgifTS+" "+TS_baseWWW);

  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"%iTeV_",sqrts);
  string prependName = "CardFragments/signalFunctionsHMP_";
  if (useggH && !useVBF) prependName += "ggH_";
  else if (!useggH && useVBF) prependName += "VBF_";
  else prependName += "";
  string appendName = ".txt";
  string outCardName =  prependName + tmp_outCardName + schannel + appendName;

  float highn1,higha1,higha2,highmean,highsigma,highgamma1,highgamma2;
  highn1=highparameter(fit_nCB_1,0);
  higha1=highparameter(fit_alphaCB_1,0);
  higha2=highparameter(fit_alphaCB_2,0);
  highmean=highparameter(fit_meanCB,0);
  highsigma=highparameter(fit_sigmaCB,0);  
  highgamma1=highparameter(fit_Gamma,1);
  highgamma2=highparameter(fit_Gamma,2);

  TString tag;
  if (useggH && !useVBF) tag = "HighMasssignalShapeggH";
  else if (!useggH && useVBF) tag = "HighMasssignalShapeVBF";
  else tag = "HighMasssignalShape";

  ofstream ofsCard;
  if (usedijet && usenondijet){
    ofsCard.open(outCardName.c_str(),fstream::out);
    ofsCard << "## signal functions --- no spaces! ##" << endl;
    ofsCard << tag << " n_CB 5"<< endl;
    ofsCard << tag << " alpha_CB " << fit_alphaCB_1->GetParameter(0) << "+(" << fit_alphaCB_1->GetParameter(1) << "*@0)+(" << fit_alphaCB_1->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB_1->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB_1->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_alphaCB_1->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << tag << " n2_CB 20" << endl;
    ofsCard << tag << " alpha2_CB " << fit_alphaCB_2->GetParameter(0) << "+(" << fit_alphaCB_2->GetParameter(1) << "*@0)+(" << fit_alphaCB_2->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB_2->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB_2->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_alphaCB_2->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << tag << " mean_CB " << fit_meanCB->GetParameter(0) << "+(" << fit_meanCB->GetParameter(1) << "*@0)+(" << fit_meanCB->GetParameter(2) << "*@0*@0)+(" << fit_meanCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_meanCB->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_meanCB->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << tag << " sigma_CB " << fit_sigmaCB->GetParameter(0) << "+(" << fit_sigmaCB->GetParameter(1) << "*@0)+(" << fit_sigmaCB->GetParameter(2) << "*@0*@0)+(" << fit_sigmaCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_sigmaCB->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_sigmaCB->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << endl;

  }


  return;
}

float WidthValue(float mHStarWidth){
  ostringstream MassString;
  MassString << mHStarWidth;
  
  ifstream widthFile("widthvalues.txt");
  
  string line;

  bool FindedMass = false;

  float Gamma_ggCal, Gamma_ZZCal, Gamma_TOTCal;
  
  while (getline(widthFile,line)) {
    if( line == "" || line[0] == '#' ) continue;
    
    if(line[0]== MassString.str()[0] && line[1]== MassString.str()[1] && line[2]== MassString.str()[2]){
      
      stringstream stringline;
      stringline << line;    
      string masschar;
      stringline>>masschar>>Gamma_ggCal>>Gamma_ZZCal>>Gamma_TOTCal;
      
      FindedMass = true;
    }
  }
  if(!FindedMass) abort();

  return Gamma_TOTCal;
}

float highparameter(TF1 *generatedfit, int gamma){
  float highzero;

  if (gamma==0){
    highzero=generatedfit->GetParameter(0);
    highzero+=(generatedfit->GetParameter(1)-generatedfit->GetParameter(6))*400;
    highzero+=(generatedfit->GetParameter(2)-generatedfit->GetParameter(7))*pow(400,2);
    highzero+=(generatedfit->GetParameter(3)-generatedfit->GetParameter(8))*pow(400,3);
    highzero+=generatedfit->GetParameter(4)*pow(400,4);
    highzero+=generatedfit->GetParameter(5)*pow(400,5);
  } else if (gamma==1){
    highzero=generatedfit->GetParameter(0);
    highzero+=(generatedfit->GetParameter(6)-generatedfit->GetParameter(2))*pow(400,2);
    highzero+=2*(generatedfit->GetParameter(7)-generatedfit->GetParameter(3))*pow(400,3);
    highzero-=3*(generatedfit->GetParameter(4))*pow(400,4);
    highzero-=4*(generatedfit->GetParameter(5))*pow(400,5);
  } else if (gamma==2){
    highzero=generatedfit->GetParameter(1);
    highzero+=2*(generatedfit->GetParameter(2)-generatedfit->GetParameter(6))*400;
    highzero+=3*(generatedfit->GetParameter(3)-generatedfit->GetParameter(7))*pow(400,2);
    highzero+=4*generatedfit->GetParameter(4)*pow(400,3);
    highzero+=5*generatedfit->GetParameter(5)*pow(400,4);
  }
  return highzero;
}
