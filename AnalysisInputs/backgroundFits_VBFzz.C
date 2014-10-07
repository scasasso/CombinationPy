/* 
 * Fit qqZZ background shapes and write parameters in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b backgroundFits_VBFzz.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

/*
  #ifndef __CINT__
  #include "RooGlobalFunc.h"
  #endif
  #include "RooRealVar.h"
  #include "RooDataSet.h"
  #include "RooGaussian.h"
  #include "RooConstVar.h"
  #include "RooChebychev.h"
  #include "RooAddPdf.h"
  #include "RooWorkspace.h"
  #include "RooPlot.h"
  #include "TCanvas.h"
  #include "TAxis.h"
  #include "TFile.h"
  #include "TH1.h"
*/


#include <iostream>
#include <iomanip>

using namespace RooFit ;


Double_t minM = 0.;
Double_t maxM = 1000.;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

void backgroundFits_VBFzz(int channel, int sqrts, int VBFtag, ofstream* interpCode=0);

// Run all final states and sqrts in one go
void backgroundFits_VBFzz() {

  gSystem->Exec("mkdir -p bkgFigs7TeV");
  gSystem->Exec("mkdir -p bkgFigs8TeV");

  ofstream txtInterp;
  string outInterpCodeName = "interpolation.txt";
  txtInterp.open(outInterpCodeName.c_str(),fstream::out);


//   backgroundFits_VBFzz(1,7,0,&txtInterp);
//   backgroundFits_VBFzz(2,7,0,&txtInterp);
  backgroundFits_VBFzz(3,7,0,&txtInterp);

//   backgroundFits_VBFzz(1,7,1,&txtInterp);
//   backgroundFits_VBFzz(2,7,1,&txtInterp);
  backgroundFits_VBFzz(3,7,1,&txtInterp);


//   backgroundFits_VBFzz(1,8,0,&txtInterp);
//   backgroundFits_VBFzz(2,8,0,&txtInterp);
//   backgroundFits_VBFzz(3,8,0,&txtInterp);

//   backgroundFits_VBFzz(1,8,1,&txtInterp);
//   backgroundFits_VBFzz(2,8,1,&txtInterp);
//   backgroundFits_VBFzz(3,8,1,&txtInterp);

}

// The actual job
void backgroundFits_VBFzz(int channel, int sqrts, int VBFtag, ofstream* interpCode)
{
  // if(sqrts==7 && channel==3) return;
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << schannel << endl;

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "schannel = " << schannel << "  sqrts = " << sqrts << " VBFtag = " << VBFtag << endl;

  TString outfile;
  if(VBFtag<2) outfile = "CardFragments/VBFzzBackgroundFit_" + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  if(VBFtag==2) outfile = "CardFragments/VBFzzBackgroundFit_" + ssqrts + "_" + schannel + ".txt";
  ofstream of(outfile,ios_base::out);
  of << "### background functions ###" << endl;

  (*interpCode) << "  if (sqrts==" << sqrts << "){" << endl;
  (*interpCode) << "    if (channel==" << channel << "){" << endl;
  (*interpCode) << "      if (VBFtag==" << VBFtag << "){" << endl;


  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  gROOT->ProcessLine(".L ../CreateDatacards/include/tdrstyle.cc");
  setTDRStyle(false);
  gStyle->SetPadLeftMargin(0.16);

  TString filepath;
  if (sqrts==7) {
      filepath = "root://lxcms00//data3/2014/HZZ_stat/140604/PRODFSR/";
      // filepath = filePath7TeV;
  } else if (sqrts==8) {
    filepath = "root://lxcms00//data3/2014/HZZ_stat/140604/PRODFSR_8TeV/";
    // filepath = filePath8TeV;
  }

  TChain* tree = new TChain("SelectedTree");

  cout << filepath << endl;
  tree->Add( filepath + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_ZZTo*JJ_Contin.root");

  RooRealVar* MC_weight = new RooRealVar("MC_weight","MC_weight",0.,2.) ; 
  RooRealVar* ZZMass = new RooRealVar("ZZMass","ZZMass",minM,maxM);
  RooRealVar* NJets30 = new RooRealVar("NJets30","NJets30",0.,100.);
  RooArgSet ntupleVarSet(*ZZMass,*NJets30,*MC_weight);
  RooDataSet *set = new RooDataSet("set","set",ntupleVarSet,WeightVar("MC_weight"));

  ZZMass->setRange("fullrange",100.,1000.);
  ZZMass->setRange("fitrange",100.,600.);
  ZZMass->setRange("zoomrange",100.,200.);    


  Float_t myMC,myMass;
  Short_t myNJets;
  int nentries = tree->GetEntries();

  tree->SetBranchAddress("ZZMass",&myMass);
  tree->SetBranchAddress("MC_weight",&myMC);
  tree->SetBranchAddress("NJets30",&myNJets);

  // cout << "nentries = " << nentries << endl;

  for(int i =0;i<nentries;i++) {
    tree->GetEntry(i);
    if(VBFtag==1 && myNJets<2)continue;
    if(VBFtag==0 && myNJets>1)continue;

    ntupleVarSet.setRealValue("ZZMass",myMass);
    ntupleVarSet.setRealValue("MC_weight",myMC);
    ntupleVarSet.setRealValue("NJets30",(double)myNJets);

    set->add(ntupleVarSet, myMC);
  }

  double totalweight = 0.;
  double totalweight_z = 0.;
  for (int i=0 ; i<set->numEntries() ; i++) { 
    //set->get(i) ; 
    RooArgSet* row = set->get(i) ;
    //row->Print("v");
    totalweight += set->weight();
    if (row->getRealValue("ZZMass") < 200) totalweight_z += set->weight();
  } 
  cout << "nEntries: " << set->numEntries() << ", totalweight: " << totalweight << ", totalweight_z: " << totalweight_z << endl;

  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
	


  //// ---------------------------------------
  //Background
  RooRealVar CMS_VBFzzbkg_a0("CMS_VBFzzbkg_a0","CMS_VBFzzbkg_a0",200.,10.,500.);
  RooRealVar CMS_VBFzzbkg_a1("CMS_VBFzzbkg_a1","CMS_VBFzzbkg_a1",200.,10.,500.);
  RooRealVar CMS_VBFzzbkg_a2("CMS_VBFzzbkg_a2","CMS_VBFzzbkg_a2",160.,10.,420.);
  RooRealVar CMS_VBFzzbkg_a3("CMS_VBFzzbkg_a3","CMS_VBFzzbkg_a3",0.15,0.,2.);
  RooRealVar CMS_VBFzzbkg_a4("CMS_VBFzzbkg_a4","CMS_VBFzzbkg_a4",185.,50.,320.);
  RooRealVar CMS_VBFzzbkg_a5("CMS_VBFzzbkg_a5","CMS_VBFzzbkg_a5",16.,2.,50.);
  RooRealVar CMS_VBFzzbkg_a6("CMS_VBFzzbkg_a6","CMS_VBFzzbkg_a6",34.,50.,150.);
  RooRealVar CMS_VBFzzbkg_a7("CMS_VBFzzbkg_a7","CMS_VBFzzbkg_a7",0.1,0.,2.);
  RooRealVar CMS_VBFzzbkg_a8("CMS_VBFzzbkg_a8","CMS_VBFzzbkg_a8",58.,10.,200.);
  RooRealVar CMS_VBFzzbkg_a9("CMS_VBFzzbkg_a9","CMS_VBFzzbkg_a9",0.5,0.,2.);
  RooRealVar CMS_VBFzzbkg_a10("CMS_VBFzzbkg_a10","CMS_VBFzzbkg_a10",170.,10.,320.);
  RooRealVar CMS_VBFzzbkg_a11("CMS_VBFzzbkg_a11","CMS_VBFzzbkg_a11",-15.,-0.,-150.);
  RooRealVar CMS_VBFzzbkg_a12("CMS_VBFzzbkg_a12","CMS_VBFzzbkg_a12",160.,10.,520.);
  RooRealVar CMS_VBFzzbkg_a13("CMS_VBFzzbkg_a13","CMS_VBFzzbkg_a13",0.2,0.,2.);


  if (VBFtag==0){
    if (channel == 1){
      ///* 4mu
      CMS_VBFzzbkg_a0.setVal(103.854);
      CMS_VBFzzbkg_a1.setVal(10.0718);
      CMS_VBFzzbkg_a2.setVal(117.551);
      CMS_VBFzzbkg_a3.setVal(0.0450287);
      CMS_VBFzzbkg_a4.setVal(185.262);
      CMS_VBFzzbkg_a5.setVal(7.99428);
      CMS_VBFzzbkg_a6.setVal(39.7813);
      CMS_VBFzzbkg_a7.setVal(0.0986891);
      CMS_VBFzzbkg_a8.setVal(49.1325);
      CMS_VBFzzbkg_a9.setVal(0.0389984);
      CMS_VBFzzbkg_a10.setVal(98.6645);
      CMS_VBFzzbkg_a11.setVal(-7.02043);
      CMS_VBFzzbkg_a12.setVal(5694.66);
      CMS_VBFzzbkg_a13.setVal(0.0774525);
      //*/
    }
    else if (channel == 2){
      ///* 4e
      CMS_VBFzzbkg_a0.setVal(111.165);
      CMS_VBFzzbkg_a1.setVal(19.8178);
      CMS_VBFzzbkg_a2.setVal(120.89);
      CMS_VBFzzbkg_a3.setVal(0.0546639);
      CMS_VBFzzbkg_a4.setVal(184.878);
      CMS_VBFzzbkg_a5.setVal(11.7041);
      CMS_VBFzzbkg_a6.setVal(33.2659);
      CMS_VBFzzbkg_a7.setVal(0.140858);
      CMS_VBFzzbkg_a8.setVal(56.1226);
      CMS_VBFzzbkg_a9.setVal(0.0957699);
      CMS_VBFzzbkg_a10.setVal(98.3662);
      CMS_VBFzzbkg_a11.setVal(-6.98701);
      CMS_VBFzzbkg_a12.setVal(10.0536);
      CMS_VBFzzbkg_a13.setVal(0.110576);
      //*/
    }
    else if (channel == 3){
      ///* 2e2mu
      CMS_VBFzzbkg_a0.setVal(110.293);
      CMS_VBFzzbkg_a1.setVal(11.8334);
      CMS_VBFzzbkg_a2.setVal(116.91);
      CMS_VBFzzbkg_a3.setVal(0.0433151);
      CMS_VBFzzbkg_a4.setVal(185.817);
      CMS_VBFzzbkg_a5.setVal(10.5945);
      CMS_VBFzzbkg_a6.setVal(29.6208);
      CMS_VBFzzbkg_a7.setVal(0.0826);
      CMS_VBFzzbkg_a8.setVal(53.1346);
      CMS_VBFzzbkg_a9.setVal(0.0882081);
      CMS_VBFzzbkg_a10.setVal(85.3776);
      CMS_VBFzzbkg_a11.setVal(-13.3836);
      CMS_VBFzzbkg_a12.setVal(7587.95);
      CMS_VBFzzbkg_a13.setVal(0.325621);
      //*/
    }
    else {
      cout << "disaster" << endl;
    }
  }


  if (VBFtag==1){
    if (channel == 1){
      ///* 4mu
      CMS_VBFzzbkg_a0.setVal(187.427);
      CMS_VBFzzbkg_a1.setVal(13.4226);
      CMS_VBFzzbkg_a2.setVal(99.4485);
      CMS_VBFzzbkg_a3.setVal(1.144);
      CMS_VBFzzbkg_a4.setVal(131.256);
      CMS_VBFzzbkg_a5.setVal(22.0388);
      CMS_VBFzzbkg_a6.setVal(53.6544);
      CMS_VBFzzbkg_a7.setVal(0.118926);
      CMS_VBFzzbkg_a8.setVal(55.4692);
      CMS_VBFzzbkg_a9.setVal(0.0693987);
      CMS_VBFzzbkg_a10.setVal(100.004);
      CMS_VBFzzbkg_a11.setVal(0);
      CMS_VBFzzbkg_a12.setVal(519.709);
      CMS_VBFzzbkg_a13.setVal(0.0824696);
      //*/
    }
    else if (channel == 2){
      ///* 4e
      CMS_VBFzzbkg_a0.setVal(111.165);
      CMS_VBFzzbkg_a1.setVal(19.8178);
      CMS_VBFzzbkg_a2.setVal(120.89);
      CMS_VBFzzbkg_a3.setVal(0.0546639);
      CMS_VBFzzbkg_a4.setVal(184.878);
      CMS_VBFzzbkg_a5.setVal(11.7041);
      CMS_VBFzzbkg_a6.setVal(33.2659);
      CMS_VBFzzbkg_a7.setVal(0.140858);
      CMS_VBFzzbkg_a8.setVal(56.1226);
      CMS_VBFzzbkg_a9.setVal(0.0957699);
      CMS_VBFzzbkg_a10.setVal(98.3662);
      CMS_VBFzzbkg_a11.setVal(-6.98701);
      CMS_VBFzzbkg_a12.setVal(10.0536);
      CMS_VBFzzbkg_a13.setVal(0.110576);
      //*/
    }
    else if (channel == 3){
      ///* 2e2mu
      CMS_VBFzzbkg_a0.setVal(110.293);
      CMS_VBFzzbkg_a1.setVal(11.8334);
      CMS_VBFzzbkg_a2.setVal(116.91);
      CMS_VBFzzbkg_a3.setVal(0.0433151);
      CMS_VBFzzbkg_a4.setVal(185.817);
      CMS_VBFzzbkg_a5.setVal(10.5945);
      CMS_VBFzzbkg_a6.setVal(29.6208);
      CMS_VBFzzbkg_a7.setVal(0.0826);
      CMS_VBFzzbkg_a8.setVal(53.1346);
      CMS_VBFzzbkg_a9.setVal(0.0882081);
      CMS_VBFzzbkg_a10.setVal(85.3776);
      CMS_VBFzzbkg_a11.setVal(-13.3836);
      CMS_VBFzzbkg_a12.setVal(7587.95);
      CMS_VBFzzbkg_a13.setVal(0.325621);
      //*/
    }
    else {
      cout << "disaster" << endl;
    }
  }


    
  RooqqZZPdf_v2* bkg_VBFzz = new RooqqZZPdf_v2("bkg_VBFzz","bkg_VBFzz",*ZZMass,
					      CMS_VBFzzbkg_a0,CMS_VBFzzbkg_a1,CMS_VBFzzbkg_a2,CMS_VBFzzbkg_a3,CMS_VBFzzbkg_a4,
					      CMS_VBFzzbkg_a5,CMS_VBFzzbkg_a6,CMS_VBFzzbkg_a7,CMS_VBFzzbkg_a8,
					      CMS_VBFzzbkg_a9,CMS_VBFzzbkg_a10,CMS_VBFzzbkg_a11,CMS_VBFzzbkg_a12,CMS_VBFzzbkg_a13);


 
  RooFitResult *r1 = bkg_VBFzz->fitTo( *set, Save(kTRUE), SumW2Error(kTRUE), Range("fitrange") );//, Save(kTRUE), SumW2Error(kTRUE)) ;


  cout << endl;
  cout << "------- Parameters for " << schannel << " sqrts=" << sqrts << endl;
  cout << "  a0_bkgd = " << CMS_VBFzzbkg_a0.getVal() << endl;
  cout << "  a1_bkgd = " << CMS_VBFzzbkg_a1.getVal() << endl;
  cout << "  a2_bkgd = " << CMS_VBFzzbkg_a2.getVal() << endl;
  cout << "  a3_bkgd = " << CMS_VBFzzbkg_a3.getVal() << endl;
  cout << "  a4_bkgd = " << CMS_VBFzzbkg_a4.getVal() << endl;
  cout << "  a5_bkgd = " << CMS_VBFzzbkg_a5.getVal() << endl;
  cout << "  a6_bkgd = " << CMS_VBFzzbkg_a6.getVal() << endl;
  cout << "  a7_bkgd = " << CMS_VBFzzbkg_a7.getVal() << endl;
  cout << "  a8_bkgd = " << CMS_VBFzzbkg_a8.getVal() << endl;
  cout << "  a9_bkgd = " << CMS_VBFzzbkg_a9.getVal() << endl;
  cout << "  a10_bkgd = " << CMS_VBFzzbkg_a10.getVal() << endl;
  cout << "  a11_bkgd = " << CMS_VBFzzbkg_a11.getVal() << endl;
  cout << "  a12_bkgd = " << CMS_VBFzzbkg_a12.getVal() << endl;
  cout << "  a13_bkgd = " << CMS_VBFzzbkg_a13.getVal() << endl;
  cout << "}" << endl;
  cout << "---------------------------" << endl;

  of << "VBFZZshape a0_bkgd   " << CMS_VBFzzbkg_a0.getVal() << endl;
  of << "VBFZZshape a1_bkgd   " << CMS_VBFzzbkg_a1.getVal() << endl;
  of << "VBFZZshape a2_bkgd   " << CMS_VBFzzbkg_a2.getVal() << endl;
  of << "VBFZZshape a3_bkgd   " << CMS_VBFzzbkg_a3.getVal() << endl;
  of << "VBFZZshape a4_bkgd   " << CMS_VBFzzbkg_a4.getVal() << endl;
  of << "VBFZZshape a5_bkgd   " << CMS_VBFzzbkg_a5.getVal() << endl;
  of << "VBFZZshape a6_bkgd   " << CMS_VBFzzbkg_a6.getVal() << endl;
  of << "VBFZZshape a7_bkgd   " << CMS_VBFzzbkg_a7.getVal() << endl;
  of << "VBFZZshape a8_bkgd   " << CMS_VBFzzbkg_a8.getVal() << endl;
  of << "VBFZZshape a9_bkgd   " << CMS_VBFzzbkg_a9.getVal() << endl;
  of << "VBFZZshape a10_bkgd  " << CMS_VBFzzbkg_a10.getVal() << endl;
  of << "VBFZZshape a11_bkgd  " << CMS_VBFzzbkg_a11.getVal() << endl;
  of << "VBFZZshape a12_bkgd  " << CMS_VBFzzbkg_a12.getVal() << endl;
  of << "VBFZZshape a13_bkgd  " << CMS_VBFzzbkg_a13.getVal() << endl;
  of << endl << endl;
  of.close();

  cout << endl << "Output written to: " << outfile << endl;


  (*interpCode) << endl << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a0.setVal(" << CMS_VBFzzbkg_a0.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a1.setVal(" << CMS_VBFzzbkg_a1.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a2.setVal(" << CMS_VBFzzbkg_a2.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a3.setVal(" << CMS_VBFzzbkg_a3.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a4.setVal(" << CMS_VBFzzbkg_a4.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a5.setVal(" << CMS_VBFzzbkg_a5.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a6.setVal(" << CMS_VBFzzbkg_a6.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a7.setVal(" << CMS_VBFzzbkg_a7.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a8.setVal(" << CMS_VBFzzbkg_a8.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a9.setVal(" << CMS_VBFzzbkg_a9.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a10.setVal(" << CMS_VBFzzbkg_a10.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a11.setVal(" << CMS_VBFzzbkg_a11.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a12.setVal(" << CMS_VBFzzbkg_a12.getVal() << ");" << endl;
  (*interpCode) << "      CMS_VBFzzbkg_a13.setVal(" << CMS_VBFzzbkg_a13.getVal() << ");" << endl;

  (*interpCode) << "      }" << endl;
  (*interpCode) << "    }" << endl;
  (*interpCode) << "  }" << endl;

  
    
  // Plot m4l and
  // RooPlot* frameM4l = ZZMass->frame(Title("M4L"),Range(100,600),Bins(250)) ;
  RooPlot* frameM4l = ZZMass->frame(Title("M4L"),Range("fitrange"),Bins(250)) ;
  set->plotOn(frameM4l, MarkerStyle(20));
  
//   //set->plotOn(frameM4l) ;
//   RooPlot* frameM4lz = ZZMass->frame(Title("M4L"),Range(100,200),Bins(100)) ;
//   set->plotOn(frameM4lz, MarkerStyle(20)) ;


  int iLineColor = 1;
  string lab = "blah";
  if (channel == 1) { iLineColor = 2; lab = "4#mu"; }
  if (channel == 3) { iLineColor = 4; lab = "2e2#mu"; }
  if (channel == 2) { iLineColor = 6; lab = "4e"; }

  bkg_VBFzz->plotOn(frameM4l,LineColor(iLineColor),NormRange("fitrange")) ;
//   bkg_VBFzz->plotOn(frameM4lz,LineColor(iLineColor),NormRange("zoomrange")) ;
  // bkg_VBFzz->plotOn(frameM4l,LineColor(iLineColor),NormRange("fullrange")) ;
    
// //second shape to compare with (if previous comparison code unceommented)
//   //bkg_VBFzz_bkgd->plotOn(frameM4l,LineColor(1),NormRange("largerange")) ;
//   //bkg_VBFzz_bkgd->plotOn(frameM4lz,LineColor(1),NormRange("zoomrange")) ;
    
  
// //   double normalizationBackground_qqzz = bkg_VBFzz->createIntegral( RooArgSet(*ZZMass), Range("fullrange") )->getVal();
// //   cout << "Norm all = " << normalizationBackground_qqzz << endl;
    
  frameM4l->GetXaxis()->SetTitle("m_{4l} [GeV]");
  frameM4l->GetYaxis()->SetTitle("a.u.");
// //   frameM4lz->GetXaxis()->SetTitle("m_{4l} [GeV]");
// //   frameM4lz->GetYaxis()->SetTitle("a.u.");

  char lname[192];
  sprintf(lname,"VBF #rightarrow ZZ #rightarrow %s", lab.c_str() );
  char lname2[192];
  sprintf(lname2,"Shape Model, %s", lab.c_str() );
  // dummy!
  TF1* dummyF = new TF1("dummyF","1",0.,1.);
  TH1F* dummyH = new TH1F("dummyH","",1, 0.,1.);
  dummyF->SetLineColor( iLineColor );
  dummyF->SetLineWidth( 2 );

  dummyH->SetLineColor( kBlue );
  TLegend * box2 = new TLegend(0.4,0.70,0.80,0.90);
  box2->SetFillColor(0);
  box2->SetBorderSize(0);
  box2->AddEntry(dummyH,"Simulation (POWHEG+Pythia)  ","pe");
  box2->AddEntry(dummyH,lname,"");
  box2->AddEntry(dummyH,"","");
  box2->AddEntry(dummyF,lname2,"l");
    
  TPaveText *pt = new TPaveText(0.15,0.955,0.4,0.99,"NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->AddText("CMS Preliminary 2012");
  TPaveText *pt2 = new TPaveText(0.84,0.955,0.99,0.99,"NDC");
  pt2->SetFillColor(0);
  pt2->SetBorderSize(0);
  TString entag;entag.Form("#sqrt{s} = %d TeV",sqrts);
  pt2->AddText(entag.Data());

  TCanvas *c = new TCanvas("c","c",800,600);
  c->cd();
  frameM4l->Draw();
  // frameM4l->GetYaxis()->SetRangeUser(0,0.4);
  // if(channel == 3)frameM4l->GetYaxis()->SetRangeUser(0,0.7);
  box2->Draw();
  pt->Draw();
  pt2->Draw();
  TString outputPath = "bkgFigs";
  outputPath = outputPath+ (long) sqrts + "TeV/";
  TString outputName;
  if(VBFtag<2) outputName =  outputPath + "bkgVBFzz_" + schannel + "_" + Form("%d",int(VBFtag));
  if(VBFtag==2) outputName =  outputPath + "bkgVBFzz_" + schannel;
  c->SaveAs(outputName + ".eps");
  c->SaveAs(outputName + ".png");
  c->SaveAs(outputName + ".root");
    
//   TCanvas *c2 = new TCanvas("c2","c2",1000,500);
//   c2->Divide(2,1);
//   c2->cd(1);
//   frameM4l->Draw();
//   box2->Draw("same");
//   c2->cd(2);
// //   frameM4lz->Draw();
// //   box2->Draw("same");
  
//   if (VBFtag<2) outputName = outputPath + "bkgVBFzz_" + schannel + "_z" + "_" + Form("%d",int(VBFtag));
//   if (VBFtag==2) outputName = outputPath + "bkgVBFzz_" + schannel + "_z";
//   c2->SaveAs(outputName + ".eps");
//   c2->SaveAs(outputName + ".png");

//   /* TO make the ratio btw 2 shapes, if needed for compairson
//   TCanvas *c3 = new TCanvas("c3","c3",1000,500);
//    if(sqrts==7)
//     sprintf(outputName, "bkgFigs7TeV/bkgVBFzz_%s_ratio.eps",schannel.c_str());
//   else if(sqrts==8)
//     sprintf(outputName, "bkgFigs8TeV/bkgVBFzz_%s_ratio.eps",schannel.c_str());

//    const int nPoints = 501.;
//   double masses[nPoints] ;
//   int j=0;
//   for (int i=100; i<601; i++){
//     masses[j] = i;
//     j++;
//   }
//   cout<<j<<endl;
//   double effDiff[nPoints];
//   for (int i = 0; i < nPoints; i++){
//     ZZMass->setVal(masses[i]);
//     double eval = (bkg_VBFzz_bkgd->getVal(otherASet)-bkg_VBFzz->getVal(myASet))/(bkg_VBFzz->getVal(myASet));
//     //cout<<bkg_VBFzz_bkgd->getVal(otherASet)<<" "<<bkg_VBFzz->getVal(myASet)<<" "<<eval<<endl;
//     effDiff[i]=eval;
//   }
//   TGraph* grEffDiff = new TGraph( nPoints, masses, effDiff );
//   grEffDiff->SetMarkerStyle(20);
//   grEffDiff->Draw("AL");

//   //c3->SaveAs(outputName);
//   */

//   if (VBFtag<2) outputName = outputPath + "bkgVBFzz_" + schannel + "_z" + "_" + Form("%d",int(VBFtag)) + ".root";
//   if (VBFtag==2) outputName = outputPath + "bkgVBFzz_" + schannel + "_z" + ".root";
//   TFile* outF = new TFile(outputName,"RECREATE");
//   outF->cd();
//   c2->Write();
//   frameM4l->Write();
// //   frameM4lz->Write();	
//   outF->Close();


  delete c;
//   delete c2;
}

