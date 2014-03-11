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
TString appendStr("_C1.0_fitted");
Double_t cprimeVal = 1.0;
bool fixResPars = false;

//Parameters to choose which samples to examine
//HCP option is useggH=true, useVBF,useVH = false and usedijet,usenondijet=true
bool debug = false;
bool useggH = true;
bool useVBF = false;
bool useVH = false;
bool usedijet = true;
bool usenondijet = true;

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

  (*interpCode) << "    if (fixResPars){ " << endl;
  (*interpCode) << "      if (sqrts==" << sqrts << "){" << endl;
  (*interpCode) << "	    if (channel==" << channel << "){" << endl;
  (*interpCode) << "  	      // Fitted on the C" << cprimeVal << endl;

  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;

  char inParFile[192];
  sprintf(inParFile,"ParamsGrid_SigInt_%iTeV.root",sqrts);
  TString TS_inParFile(inParFile);

  TFile* fPar = TFile::Open(TS_inParFile,"READ");
  TH2D* hR = (TH2D*)fPar->Get("h_r");
  TH2D* hAlpha = (TH2D*)fPar->Get("h_Alpha");
  TH2D* hBeta = (TH2D*)fPar->Get("h_Beta");
  TH2D* hGamma = (TH2D*)fPar->Get("h_Gamma");
  TH2D* hDelta = (TH2D*)fPar->Get("h_Delta");


  //Pick the correct mass points and paths
  TString filePath;
  int nPointsggH,nPointsVBF,nPointsVH;
  double* massesggH,*massesVBF,*massesVH;

  if (sqrts==7) {
    nPointsggH = nPoints7TeV_p15_HM;
    massesggH  = mHVal7TeV_p15_HM;
    nPointsVBF = nVBFPoints7TeV;
    massesVBF  = mHVBFVal7TeV;
    nPointsVH = nVHPoints7TeV;
    massesVH  = mHVHVal7TeV;
    filePath = filePath7TeV;  
  } else if (sqrts==8) {
    nPointsggH = nPoints8TeV_p15_HM;
    massesggH  = mHVal8TeV_p15_HM;
    nPointsVBF = nVBFPoints8TeV;
    massesVBF  = mHVBFVal8TeV;
    nPointsVH = nVHPoints8TeV;
    massesVH  = mHVHVal8TeV;
    filePath = filePath8TeV;
  }
  else abort();

  int nPoints;
  double masses[200];
  for (int i=0;i<200;i++){
    masses[i]=-1;
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
    sprintf(tmp_finalInPathVBF,"/HZZ4lTree_VBFH%i.root",masses[i]);
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
    Double_t gammaVal = hGamma->GetBinContent(hBeta->FindBin(masses[i],cprimeVal));
    double windowVal = max(valueWidth,1.);
    double lowside = 100.;
    if(masses[i] > 300) lowside = 200.;
    double low_M = max( (masses[i] - 15.*windowVal), lowside) ;
    double high_M = min( (masses[i] + 10.*windowVal), 1400.);

    // FIXME: as soon as fit ranges can be defined in a continuous way, replace this twofold ranges definition
    if(masses[i] > 399.){

      if (cprimeVal>0.9){
	// For the C'^2 = 1.0 case
	low_M = max( (masses[i] - 2.*windowVal), 250.) ;
	// high_M = min( (masses[i] + 2.*windowVal), 1600.);
	high_M = min( (masses[i] + 2.*windowVal), 1300.);
	if (fabs(masses[i]-800)<5) high_M = 1150;
      }
      else{
	// For the C'^2 = 0.2 case
	if (masses[i]<401.){
	  low_M = 360;
	  high_M = 430; 
	}
	else if (masses[i]<501.){
	  low_M = 440;
	  high_M = 540; 	
	}
	else if (masses[i]<601.){
	  low_M = 520.;
	  high_M = 660.; 	
	}
	else if (masses[i]<701.){
	  low_M = 600.;
	  high_M = 760.; 	
	}
	else if (masses[i]<801.){
	  low_M = 680.;
	  high_M = 860.; 	
	}
	else if (masses[i]<901.){
	  low_M = 730.;
	  high_M = 970.; 	
	}
	else if (masses[i]<1001.){
	  low_M = 800.;
	  high_M = 1100.; 	
	}
	else {
	  low_M = 300;
	  high_M = 1600; 	
	}
      }
      
    }

    cout << "lowM = " << low_M << ", highM = " << high_M << endl;

    //Set the observable and get the RooDataSomething
    RooRealVar ZZMass("ZZMass","ZZMass",low_M,high_M);
    ZZMass.setRange("fitRange",low_M,high_M);
    ZZMass.setRange("plotRange",max( (masses[i] - 2.*windowVal), 250.),min( (masses[i] + 2.*windowVal), 1300.));
    RooRealVar MC_weight("MC_weight","MC_weight",0.,10.);
    RooRealVar NJets("NJets","NJets",0.,100.);
    RooRealVar genProcessId("genProcessId","genProcessId",0.,150.);

    if(channel == 2) ZZMass.setBins(50);
    if(channel == 3) ZZMass.setBins(50);

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


    Double_t rVal = hR->GetBinContent(hR->FindBin(masses[i],cprimeVal));
    Double_t alphaVal = hAlpha->GetBinContent(hAlpha->FindBin(masses[i],cprimeVal));
    Double_t betaVal = hBeta->GetBinContent(hBeta->FindBin(masses[i],cprimeVal));
    Double_t deltaVal = hDelta->GetBinContent(hDelta->FindBin(masses[i],cprimeVal));

    RooRealVar k("k","k",0.25);
    RooRealVar delta("delta","delta",deltaVal);
    RooRealVar CSquared("CSquared","C'^{2}",cprimeVal);
    RooRealVar BRnew("BRnew","BR_{new}",0.);
    RooRealVar alpha("alpha","#alpha",alphaVal);
    RooRealVar beta("beta","#beta",betaVal);
    RooRealVar r("r","r",rVal,rVal-0.2*rVal,rVal+0.2*rVal);
    r.setConstant(kTRUE);

    RooSigPlusInt SignalTheor("model","model",ZZMass,MHStar,delta,Gamma_TOT,k,CSquared,BRnew,alpha,beta,r);

    //Experimental resolution
    RooRealVar meanCB("meanCB","meanCB",0.,-0.8,0.8);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.,5.);
    RooRealVar sigmaCB_high("sigmaCB_high","sigmaCB_high",6.5,3.,10.);
    RooRealVar alphaCB_1("alphaCB_1","alphaCB_1",1.,0.4,2.);
    RooRealVar nCB_1("nCB_1","nCB_1",5.,0.,12.);
    nCB_1.setConstant(kTRUE);
    RooRealVar alphaCB_2("alphaCB_2","alphaCB_2",1.,0.4,2.);
    RooRealVar nCB_2("nCB_2","nCB_2",20.,0.,12.);
    nCB_2.setConstant(kTRUE);

    //Initialize to decent values
    float m = masses[i];


    // Fixing the values to the one extracted from the C'^2 = 0.2 case
    if (fixResPars){
      if (sqrts==7){
	if (channel==1){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(((0.786102+(0.00196526-0.0437839)*400+(4.91313e-06--5.26213e-05)*pow(400,2)+(1.22829e-08-2.114e-08)*pow(400,3)+3.07071e-11*pow(400,4)+7.67677e-14*pow(400,5))+0.0437839*m+-5.26213e-05*m*m+2.114e-08*m*m*m)); sigmaCB_high.setConstant(kTRUE);
	  meanCB.setVal(((-0.0295214+(-7.38034e-05-0.0292789)*400+(-1.84508e-07--4.44546e-05)*pow(400,2)+(-4.61271e-10-2.15402e-08)*pow(400,3)+-1.15318e-12*pow(400,4)+-2.88308e-15*pow(400,5))+0.0292789*m+-4.44546e-05*m*m+2.15402e-08*m*m*m)); meanCB.setConstant(kTRUE);
	  alphaCB_1.setVal(((0.183363+(0.000458408-0.00590286)*400+(1.14601e-06--1.06111e-05)*pow(400,2)+(2.86506e-09-5.61392e-09)*pow(400,3)+7.16261e-12*pow(400,4)+1.79065e-14*pow(400,5))+0.00590286*m+-1.06111e-05*m*m+5.61392e-09*m*m*m)); alphaCB_1.setConstant(kTRUE);
	  alphaCB_2.setVal(((0.200009+(0.000500024--0.00136617)*400+(1.25005e-06--3.27474e-07)*pow(400,2)+(3.12516e-09-4.49393e-10)*pow(400,3)+7.81288e-12*pow(400,4)+1.95321e-14*pow(400,5))+-0.00136617*m+-3.27474e-07*m*m+4.49393e-10*m*m*m)); alphaCB_2.setConstant(kTRUE);
	}
	else if (channel==2){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(((0.816228+(0.00204057--0.137186)*400+(5.10143e-06-0.000227819)*pow(400,2)+(1.27536e-08--1.16673e-07)*pow(400,3)+3.18839e-11*pow(400,4)+7.97096e-14*pow(400,5))+-0.137186*m+0.000227819*m*m+-1.16673e-07*m*m*m)); sigmaCB_high.setConstant(kTRUE);
	  meanCB.setVal(((-0.0135269+(-3.38163e-05--0.0135245)*400+(-8.45421e-08-2.33866e-05)*pow(400,2)+(-2.11354e-10--1.14332e-08)*pow(400,3)+-5.28385e-13*pow(400,4)+-1.32114e-15*pow(400,5))+-0.0135245*m+2.33866e-05*m*m+-1.14332e-08*m*m*m)); meanCB.setConstant(kTRUE);
	  alphaCB_1.setVal(((0.197481+(0.000493703--0.0162022)*400+(1.23426e-06-2.3071e-05)*pow(400,2)+(3.08564e-09--1.00702e-08)*pow(400,3)+7.7141e-12*pow(400,4)+1.92851e-14*pow(400,5))+-0.0162022*m+2.3071e-05*m*m+-1.00702e-08*m*m*m)); alphaCB_1.setConstant(kTRUE);
	  alphaCB_2.setVal(((0.256949+(0.000642374--0.140568)*400+(1.60593e-06-0.000213493)*pow(400,2)+(4.01483e-09--1.01699e-07)*pow(400,3)+1.00371e-11*pow(400,4)+2.50925e-14*pow(400,5))+-0.140568*m+0.000213493*m*m+-1.01699e-07*m*m*m)); alphaCB_2.setConstant(kTRUE);
	}
	else{
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(((0.785874+(0.00196469-0.0442967)*400+(4.91171e-06--5.34087e-05)*pow(400,2)+(1.22793e-08-2.15273e-08)*pow(400,3)+3.06982e-11*pow(400,4)+7.67454e-14*pow(400,5))+0.0442967*m+-5.34087e-05*m*m+2.15273e-08*m*m*m)); sigmaCB_high.setConstant(kTRUE);
 	  meanCB.setVal(((-0.0297066+(-7.42664e-05-0.0301311)*400+(-1.85666e-07--4.58667e-05)*pow(400,2)+(-4.64165e-10-2.22915e-08)*pow(400,3)+-1.16041e-12*pow(400,4)+-2.90105e-15*pow(400,5))+0.0301311*m+-4.58667e-05*m*m+2.22915e-08*m*m*m)); meanCB.setConstant(kTRUE);
	  alphaCB_1.setVal(((0.183317+(0.000458293-0.00616808)*400+(1.14573e-06--1.10595e-05)*pow(400,2)+(2.86432e-09-5.8565e-09)*pow(400,3)+7.16079e-12*pow(400,4)+1.7902e-14*pow(400,5))+0.00616808*m+-1.10595e-05*m*m+5.8565e-09*m*m*m)); alphaCB_1.setConstant(kTRUE);
	  alphaCB_2.setVal(((0.200336+(0.000500842--0.001427)*400+(1.25211e-06--3.27917e-07)*pow(400,2)+(3.13025e-09-5.02815e-10)*pow(400,3)+7.82563e-12*pow(400,4)+1.95641e-14*pow(400,5))+-0.001427*m+-3.27917e-07*m*m+5.02815e-10*m*m*m)); alphaCB_2.setConstant(kTRUE);
	}
      }// 7TeV
      
      else{
	if (channel==1){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(((1.12482+(0.00281204-0.0241159)*400+(7.0301e-06--3.83572e-05)*pow(400,2)+(1.75753e-08-1.93221e-08)*pow(400,3)+4.39381e-11*pow(400,4)+1.09845e-13*pow(400,5))+0.0241159*m+-3.83572e-05*m*m+1.93221e-08*m*m*m)); sigmaCB_high.setConstant(kTRUE);
	  meanCB.setVal(((-0.14709+(-0.000367724--0.00392939)*400+(-9.19309e-07-1.67077e-05)*pow(400,2)+(-2.29826e-09--1.13231e-08)*pow(400,3)+-5.74565e-12*pow(400,4)+-1.43642e-14*pow(400,5))+-0.00392939*m+1.67077e-05*m*m+-1.13231e-08*m*m*m)); meanCB.setConstant(kTRUE);
	  alphaCB_1.setVal(((0.247082+(0.000617706--0.00859185)*400+(1.54427e-06-5.71764e-06)*pow(400,2)+(3.86068e-09--5.80474e-10)*pow(400,3)+9.65169e-12*pow(400,4)+2.41292e-14*pow(400,5))+-0.00859185*m+5.71764e-06*m*m+-5.80474e-10*m*m*m)); alphaCB_1.setConstant(kTRUE);
	  alphaCB_2.setVal(((0.257664+(0.00064416--0.0149572)*400+(1.6104e-06-1.56168e-05)*pow(400,2)+(4.02601e-09--5.65521e-09)*pow(400,3)+1.0065e-11*pow(400,4)+2.51625e-14*pow(400,5))+-0.0149572*m+1.56168e-05*m*m+-5.65521e-09*m*m*m)); alphaCB_2.setConstant(kTRUE);
	}
	else if (channel==2){
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(((0.745765+(0.00186441--0.0156989)*400+(4.66103e-06-4.34532e-05)*pow(400,2)+(1.16526e-08--2.89894e-08)*pow(400,3)+2.91315e-11*pow(400,4)+7.28285e-14*pow(400,5))+-0.0156989*m+4.34532e-05*m*m+-2.89894e-08*m*m*m)); sigmaCB_high.setConstant(kTRUE);
	  meanCB.setVal(((-0.0386066+(-9.65166e-05-0.0457496)*400+(-2.41291e-07--6.74899e-05)*pow(400,2)+(-6.03229e-10-3.22642e-08)*pow(400,3)+-1.50807e-12*pow(400,4)+-3.77035e-15*pow(400,5))+0.0457496*m+-6.74899e-05*m*m+3.22642e-08*m*m*m)); meanCB.setConstant(kTRUE);
	  alphaCB_1.setVal(((0.190398+(0.000475994-0.00910674)*400+(1.18999e-06--1.60194e-05)*pow(400,2)+(2.97496e-09-8.7402e-09)*pow(400,3)+7.43741e-12*pow(400,4)+1.85934e-14*pow(400,5))+0.00910674*m+-1.60194e-05*m*m+8.7402e-09*m*m*m)); alphaCB_1.setConstant(kTRUE);
	  alphaCB_2.setVal(((0.196422+(0.000491054--0.000297765)*400+(1.22763e-06--2.33901e-07)*pow(400,2)+(3.06908e-09--2.75437e-10)*pow(400,3)+7.67272e-12*pow(400,4)+1.91816e-14*pow(400,5))+-0.000297765*m+-2.33901e-07*m*m+-2.75437e-10*m*m*m)); alphaCB_2.setConstant(kTRUE);
	}
	else {
	  // Fitted on the C0.2
	  sigmaCB_high.setVal(((0.99157+(0.00247892-0.0693237)*400+(6.19731e-06--9.30082e-05)*pow(400,2)+(1.54933e-08-3.92305e-08)*pow(400,3)+3.87332e-11*pow(400,4)+9.6833e-14*pow(400,5))+0.0693237*m+-9.30082e-05*m*m+3.92305e-08*m*m*m)); sigmaCB_high.setConstant(kTRUE);
	  meanCB.setVal(((-0.0271076+(-6.7769e-05-0.0010161)*400+(-1.69422e-07-7.82239e-08)*pow(400,2)+(-4.23556e-10--6.91491e-10)*pow(400,3)+-1.05889e-12*pow(400,4)+-2.64723e-15*pow(400,5))+0.0010161*m+7.82239e-08*m*m+-6.91491e-10*m*m*m)); meanCB.setConstant(kTRUE);
	  alphaCB_1.setVal(((0.215028+(0.00053757-0.00923058)*400+(1.34393e-06--1.7239e-05)*pow(400,2)+(3.35981e-09-9.14737e-09)*pow(400,3)+8.39954e-12*pow(400,4)+2.09988e-14*pow(400,5))+0.00923058*m+-1.7239e-05*m*m+9.14737e-09*m*m*m)); alphaCB_1.setConstant(kTRUE);
	  alphaCB_2.setVal(((0.266434+(0.000666085--0.0110694)*400+(1.66521e-06-1.08711e-05)*pow(400,2)+(4.16304e-09--3.8998e-09)*pow(400,3)+1.04076e-11*pow(400,4)+2.60189e-14*pow(400,5))+-0.0110694*m+1.08711e-05*m*m+-3.8998e-09*m*m*m)); alphaCB_2.setConstant(kTRUE);
	}
      } //8TeV

    }


    RooDoubleCB massRes("massRes","Double Crystal Ball",ZZMass,meanCB,sigmaCB,alphaCB_1,nCB_1,alphaCB_2,nCB_2);
    RooDoubleCB massResH("massResH","DCB Highmass",ZZMass,meanCB,sigmaCB_high,alphaCB_1,nCB_1,alphaCB_2,nCB_2);

    //Convolute theoretical shape and resolution
    RooFFTConvPdf *sigPDF;
    if(masses[i] < 399.) sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,SignalTheor,massRes);
    else sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,SignalTheor,massResH);
    sigPDF->setBufferFraction(0.2);

    RooPlot *xplot = ZZMass.frame();
    TCanvas *canv = new TCanvas("canv","canv",1200,800);

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
    sigPDF->plotOn(xplot);
    // sigPDF->paramOn(xplot);
    canv->cd();
    xplot->Draw();

    TString plotFileTitleTS(plotFileTitle.c_str());
    TString plotgif = plotFileTitleTS + ""+appendStr+".png";
    canv->SaveAs(plotgif);
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

  // TF1 *paramfit = new TF1("paramfit","(x<400)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>=400)*(([0]+([1]-[6])*400+([2]-[7])*pow(400,2)+([3]-[8])*pow(400,3)+[4]*pow(400,4)+[5]*pow(400,5))+[6]*x+[7]*x*x+[8]*x*x*x)",115,1000);
  // TF1 *paramfit = new TF1("paramfit","(([0]+([1]-[6])*400+([2]-[7])*pow(400,2)+([3]-[8])*pow(400,3)+[4]*pow(400,4)+[5]*pow(400,5))+[6]*x+[7]*x*x+[8]*x*x*x)",400,1000);
  TF1 *paramfit = new TF1("paramfit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x",400,1000);
  // TF1 *linearfit = new TF1("linearfit","[0]+[1]*x",400,1000);
  TF1 *gammafit = new TF1("gammafit","(x<400)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>=400)*(([0]+([6]-[2])*pow(400,2)+2*([7]-[3])*pow(400,3)-3*[4]*pow(400,4)-4*[5]*pow(400,5))+([1]+2*([2]-[6])*400+3*([3]-[7])*pow(400,2)+4*[4]*pow(400,3)+5*[5]*pow(400,4))*x+[6]*x*x+[7]*x*x*x)",115,1000);

  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the mean of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_meanCB->Fit("paramfit");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the sigma of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_sigmaCB->Fit("paramfit");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the alpha1 (R) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_alphaCB_1->Fit("paramfit");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the N1 (R) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_nCB_1->Fit("paramfit");
  cout<<"#################################"<<endl;
  cout<<"Fitting the trend of the alpha2 (L) of the CB"<<endl;
  cout<<"#################################"<<endl;
  gr_alphaCB_2->Fit("paramfit");
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

  
  (*interpCode) << "              sigmaCB_high.setVal("<<fit_sigmaCB->GetParameter(0)<<"+"<<fit_sigmaCB->GetParameter(1)<<"*m+"<<fit_sigmaCB->GetParameter(2)<<"*m*m+"<<fit_sigmaCB->GetParameter(3)<<"*m*m*m)"<<fit_sigmaCB->GetParameter(4)<<"*m*m*m*m)"<<fit_sigmaCB->GetParameter(5)<<"*m*m*m*m*m)); sigmaCB_high.setConstant(kTRUE);"<<endl;
  (*interpCode) << "              meanCB.setVal("<<fit_meanCB->GetParameter(0)<<"+"<<fit_meanCB->GetParameter(1)<<"*m+"<<fit_meanCB->GetParameter(2)<<"*m*m+"<<fit_meanCB->GetParameter(3)<<"*m*m*m)"<<fit_meanCB->GetParameter(4)<<"*m*m*m*m)"<<fit_meanCB->GetParameter(5)<<"*m*m*m*m*m)"<<"); meanCB.setConstant(kTRUE);"<<endl;
  (*interpCode) << "              alphaCB_1.setVal("<<fit_alphaCB_1->GetParameter(0)<<"+"<<fit_alphaCB_1->GetParameter(1)<<"*m+"<<fit_alphaCB_1->GetParameter(2)<<"*m*m+"<<fit_alphaCB_1->GetParameter(3)<<"*m*m*m)"<<fit_alphaCB_1->GetParameter(4)<<"*m*m*m*m)"<<fit_alphaCB_1->GetParameter(5)<<"*m*m*m*m*m)"<<"); alphaCB_1.setConstant(kTRUE);"<<endl;
  (*interpCode) << "              alphaCB_2.setVal("<<fit_alphaCB_2->GetParameter(0)<<"+"<<fit_alphaCB_2->GetParameter(1)<<"*m+"<<fit_alphaCB_2->GetParameter(2)<<"*m*m+"<<fit_alphaCB_2->GetParameter(3)<<"*m*m*m)"<<fit_alphaCB_2->GetParameter(4)<<"*m*m*m*m)"<<fit_alphaCB_2->GetParameter(5)<<"*m*m*m*m*m)"<<"); alphaCB_2.setConstant(kTRUE);"<<endl;

  (*interpCode) << "        }" << endl;
  (*interpCode) << "      }" << endl;
  (*interpCode) << "    }" << endl;

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

  TCanvas *canv2 = new TCanvas("canv2","canv2",1600,800);
  canv2->Divide(4,2);

  canv2->cd(1); gr_meanCB->Draw("A*");  fit_meanCB->Draw("SAME");
  canv2->cd(2); gr_sigmaCB->Draw("A*"); fit_sigmaCB->Draw("SAME");
  canv2->cd(3); gr_alphaCB_1->Draw("A*"); fit_alphaCB_1->Draw("SAME");
  canv2->cd(4); gr_nCB_1->Draw("A*");     fit_nCB_1->Draw("SAME");
  canv2->cd(5); gr_alphaCB_2->Draw("A*"); fit_alphaCB_2->Draw("SAME");
  canv2->cd(6); gr_nCB_2->Draw("A*");     fit_nCB_2->Draw("SAME");
  canv2->cd(7); gr_Gamma->Draw("A*");   fit_Gamma->Draw("SAME");
  gr_meanCB->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_sigmaCB->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_alphaCB_1->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_alphaCB_2->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_nCB_1->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_nCB_2->GetXaxis()->SetTitle("m_{H} (GeV)");
  gr_Gamma->GetXaxis()->SetTitle("m_{H} (GeV)");

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
  TString paramgifTS(paramPlotFileTitle.c_str());
  if (usedijet && !usenondijet){
    paramgifTS+="_0"+appendStr+".png";
  }else if (usenondijet && !usedijet){
    paramgifTS+="_1"+appendStr+".png";
  }else if (usenondijet && usedijet){
    paramgifTS+="deriv5"+appendStr+".png";
  }
  canv2->SaveAs(paramgifTS);

  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"%iTeV_",sqrts);
  string prependName = "CardFragments/signalFunctions_";
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

  ofstream ofsCard;
  if (usedijet && usenondijet){
    ofsCard.open(outCardName.c_str(),fstream::out);
    ofsCard << "## signal functions --- no spaces! ##" << endl;
    // ofsCard << "HighMasssignalShape n_CB " << fit_nCB_1->GetParameter(0) << "+(" << fit_nCB_1->GetParameter(1) << "*@0)+(" << fit_nCB_1->GetParameter(2) << "*@0*@0)+(" << fit_nCB_1->GetParameter(3) << "*@0*@0*@0)+(" << fit_nCB_1->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_nCB_1->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape n_CB 5"<< endl;
    ofsCard << "HighMasssignalShape alpha_CB " << fit_alphaCB_1->GetParameter(0) << "+(" << fit_alphaCB_1->GetParameter(1) << "*@0)+(" << fit_alphaCB_1->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB_1->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB_1->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_alphaCB_1->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape n2_CB 20" << endl;
    ofsCard << "HighMasssignalShape alpha2_CB " << fit_alphaCB_2->GetParameter(0) << "+(" << fit_alphaCB_2->GetParameter(1) << "*@0)+(" << fit_alphaCB_2->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB_2->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB_2->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_alphaCB_2->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape mean_CB " << fit_meanCB->GetParameter(0) << "+(" << fit_meanCB->GetParameter(1) << "*@0)+(" << fit_meanCB->GetParameter(2) << "*@0*@0)+(" << fit_meanCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_meanCB->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_meanCB->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape sigma_CB " << fit_sigmaCB->GetParameter(0) << "+(" << fit_sigmaCB->GetParameter(1) << "*@0)+(" << fit_sigmaCB->GetParameter(2) << "*@0*@0)+(" << fit_sigmaCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_sigmaCB->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_sigmaCB->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
//     ofsCard << "HighMasssignalShape n_CB " << highn1 << "+(" << fit_nCB_1->GetParameter(6) << "*@0)+(" << fit_nCB_1->GetParameter(7) << "*@0*@0)+(" << fit_nCB_1->GetParameter(8) << "*@0*@0*@0)" <<endl;
//     ofsCard << "HighMasssignalShape alpha_CB " << higha1 << "+(" << fit_alphaCB_1->GetParameter(6) << "*@0)+(" << fit_alphaCB_1->GetParameter(7) << "*@0*@0)+(" << fit_alphaCB_1->GetParameter(8) << "*@0*@0*@0)" << endl;
//     ofsCard << "HighMasssignalShape n2_CB " << fit_nCB_2->GetParameter(0) << endl;
//     ofsCard << "HighMasssignalShape alpha2_CB " << higha2 << "+(" << fit_alphaCB_2->GetParameter(6) << "*@0)+(" << fit_alphaCB_2->GetParameter(7) << "*@0*@0)+(" << fit_alphaCB_2->GetParameter(8) << "*@0*@0*@0)" << endl;
//     ofsCard << "HighMasssignalShape mean_CB " << highmean << "+(" << fit_meanCB->GetParameter(6) << "*@0)+(" << fit_meanCB->GetParameter(7) << "*@0*@0)+(" << fit_meanCB->GetParameter(8) << "*@0*@0*@0)" << endl;
//     ofsCard << "HighMasssignalShape sigma_CB " << highsigma << "+(" << fit_sigmaCB->GetParameter(6) << "*@0)+(" << fit_sigmaCB->GetParameter(7) << "*@0*@0)+(" << fit_sigmaCB->GetParameter(8) << "*@0*@0*@0)" << endl;
//     ofsCard << "HighMasssignalShape gamma_BW " << highgamma1 << "+(" << highgamma2 << "*@0)+(" << fit_Gamma->GetParameter(6) << "*@0*@0)+(" << fit_Gamma->GetParameter(7) << "*@0*@0*@0)" << endl;
    ofsCard << endl;

  }

  fPar->Close();

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
