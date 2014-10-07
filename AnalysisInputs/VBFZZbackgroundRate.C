/* 
 * Compute ZZ background rates and write them in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b VBFZZbackgroundRate.C+
 * This runs on 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLegend.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include "TTree.h"
#include "TText.h"
#include "TROOT.h"
#include "TStyle.h"

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//VBFtag convention: 0->0/1 Jets; 1->2+ jets; 2->No tagging

Double_t minM = 100.;
Double_t maxM = 1600.;

using namespace std;
void compute(TString name, TString name_Phantom, int sqrts, double lumi, int VBFtag);
double sumWeights(TString name, double lumi, int selVBF, Double_t minM=0., Double_t maxM=2000.);


// Run both sqrts in one go
void VBFZZbackgroundRate() {

//   compute(filePath7TeV, filePath7TeV_Phantom, 7, lumi7TeV, 1);
//   compute(filePath8TeV, filePath8TeV_Phantom, 8, lumi8TeV, 1);

//   compute(filePath7TeV, filePath7TeV_Phantom, 7, lumi7TeV, 0);
//   compute(filePath8TeV, filePath8TeV_Phantom, 8, lumi8TeV, 0);

  compute(filePath7TeV, filePath7TeV_Phantom, 7, lumi7TeV, 2);
  compute(filePath8TeV, filePath8TeV_Phantom, 8, lumi8TeV, 2);

}


// The actual job
void compute(TString filePath, TString filePath_Phantom, int sqrts, double lumi, int VBFtag){

  double VBFZZ_h126Win[3];
  double VBFHZZ_h126Win[3];
  double VBFZZ[3];

  double VBFHZZ_yield[3];
  if (sqrts==7){    
    VBFHZZ_yield[0] = 9.2458836E-02;
    VBFHZZ_yield[1] = 5.1755897E-02;
    VBFHZZ_yield[2] = 1.2861921E-01;
  }
  else if (sqrts==8){
    VBFHZZ_yield[0] = 4.6798807E-01;
    VBFHZZ_yield[1] = 2.4788553E-01;
    VBFHZZ_yield[2] = 6.1781689E-01;
  }
  else {
    cout << "Invalid sqrts, choose 7 or 8: " << sqrts << endl;
    exit(1);
  }

  
  cout<<"VBFZZ 4mu "<< "sqrts = " << sqrts <<  " VBFtag = " << VBFtag << endl;
  double VBFZZ_4mu_ev4mu = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo4muJJ_Contin.root",lumi,VBFtag,minM,maxM);
  double VBFZZ_4mu_ev4e = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo4eJJ_Contin.root",lumi,VBFtag,minM,maxM);
  // double VBFZZ_4mu_ev4e = 0;
  double VBFZZ_4mu_ev2e2mu = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo2e2muJJ_Contin.root",lumi,VBFtag,minM,maxM);
  // double VBFZZ_4mu_ev2e2mu = 0;
  double VBFZZ_4mu_ev4mu_h126Win = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo4muJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  double VBFZZ_4mu_ev4e_h126Win = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo4eJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  // double VBFZZ_4mu_ev4e_h126Win = 0;
  double VBFZZ_4mu_ev2e2mu_h126Win = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo2e2muJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  // double VBFZZ_4mu_ev2e2mu_h126Win = 0;
  double VBFHZZ_4mu_ev4mu_h126Win = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo4muJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  double VBFHZZ_4mu_ev4e_h126Win = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo4eJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  double VBFHZZ_4mu_ev2e2mu_h126Win = sumWeights(filePath_Phantom + "/4mu/HZZ4lTree_ZZTo2e2muJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  VBFZZ_h126Win[0] = VBFZZ_4mu_ev4mu_h126Win+VBFZZ_4mu_ev4e_h126Win+VBFZZ_4mu_ev2e2mu_h126Win;
  VBFHZZ_h126Win[0] = VBFHZZ_4mu_ev4mu_h126Win+VBFHZZ_4mu_ev4e_h126Win+VBFHZZ_4mu_ev2e2mu_h126Win;
  // VBFZZ[0] = (VBFZZ_4mu_ev4mu+VBFZZ_4mu_ev4e+VBFZZ_4mu_ev2e2mu)*(VBFHZZ_yield[0]/(VBFHZZ_h126Win[0]-VBFZZ_h126Win[0]));    
  VBFZZ[0] = VBFZZ_4mu_ev4mu+VBFZZ_4mu_ev4e+VBFZZ_4mu_ev2e2mu;    
  
  cout<<VBFZZ_4mu_ev4mu<<" "<<VBFZZ_4mu_ev4e<<" "<<VBFZZ_4mu_ev2e2mu
      <<" -> ("<<VBFZZ_4mu_ev4mu<<" + "
      <<VBFZZ_4mu_ev4e+VBFZZ_4mu_ev2e2mu<<")*("<<VBFHZZ_yield[0]/(VBFHZZ_h126Win[0]-VBFZZ_h126Win[0])<<")"
      <<" -> "<<VBFZZ[0]<<endl;

  cout<<"VBFZZ 4e "<< "sqrts = " << sqrts <<  " VBFtag = " << VBFtag << endl;
  double VBFZZ_4e_ev4mu = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo4muJJ_Contin.root",lumi,VBFtag,minM,maxM);
  double VBFZZ_4e_ev4e = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo4eJJ_Contin.root",lumi,VBFtag,minM,maxM);
  // double VBFZZ_4e_ev4e = 0;
  double VBFZZ_4e_ev2e2mu = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo2e2muJJ_Contin.root",lumi,VBFtag,minM,maxM);
  // double VBFZZ_4e_ev2e2mu = 0;
  double VBFZZ_4e_ev4mu_h126Win = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo4muJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  double VBFZZ_4e_ev4e_h126Win = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo4eJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  // double VBFZZ_4e_ev4e_h126Win = 0;
  double VBFZZ_4e_ev2e2mu_h126Win = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo2e2muJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  // double VBFZZ_4e_ev2e2mu_h126Win = 0;
  double VBFHZZ_4e_ev4mu_h126Win = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo4muJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  double VBFHZZ_4e_ev4e_h126Win = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo4eJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  double VBFHZZ_4e_ev2e2mu_h126Win = sumWeights(filePath_Phantom + "/4e/HZZ4lTree_ZZTo2e2muJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  VBFZZ_h126Win[1] = VBFZZ_4e_ev4mu_h126Win+VBFZZ_4e_ev4e_h126Win+VBFZZ_4e_ev2e2mu_h126Win;
  VBFHZZ_h126Win[1] = VBFHZZ_4e_ev4mu_h126Win+VBFHZZ_4e_ev4e_h126Win+VBFHZZ_4e_ev2e2mu_h126Win;
  // VBFZZ[1] = (VBFZZ_4e_ev4mu+VBFZZ_4e_ev4e+VBFZZ_4e_ev2e2mu)*(VBFHZZ_yield[1]/(VBFHZZ_h126Win[1]-VBFZZ_h126Win[1]));    
  VBFZZ[1] = VBFZZ_4e_ev4mu+VBFZZ_4e_ev4e+VBFZZ_4e_ev2e2mu;    
  
  cout<<VBFZZ_4e_ev4mu<<" "<<VBFZZ_4e_ev4e<<" "<<VBFZZ_4e_ev2e2mu
      <<" -> ("<<VBFZZ_4e_ev4e<<" + "
      <<VBFZZ_4e_ev4mu+VBFZZ_4e_ev2e2mu<<")*("<<VBFHZZ_yield[1]/(VBFHZZ_h126Win[1]-VBFZZ_h126Win[1])<<")"
      <<" -> "<<VBFZZ[1]<<endl;


  cout<<"VBFZZ 2e2mu "<< "sqrts = " << sqrts <<  " VBFtag = " << VBFtag << endl;
  double VBFZZ_2e2mu_ev4mu = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo4muJJ_Contin.root",lumi,VBFtag,minM,maxM);
  double VBFZZ_2e2mu_ev4e = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo4eJJ_Contin.root",lumi,VBFtag,minM,maxM);
  // double VBFZZ_2e2mu_ev4e = 0;
  double VBFZZ_2e2mu_ev2e2mu = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo2e2muJJ_Contin.root",lumi,VBFtag,minM,maxM);
  // double VBFZZ_2e2mu_ev2e2mu = 0;
  double VBFZZ_2e2mu_ev4mu_h126Win = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo4muJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  double VBFZZ_2e2mu_ev4e_h126Win = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo4eJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  // double VBFZZ_2e2mu_ev4e_h126Win = 0;
  double VBFZZ_2e2mu_ev2e2mu_h126Win = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo2e2muJJ_Contin.root",lumi,VBFtag,105.6,140.6);
  // double VBFZZ_2e2mu_ev2e2mu_h126Win = 0;
  double VBFHZZ_2e2mu_ev4mu_h126Win = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo4muJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  double VBFHZZ_2e2mu_ev4e_h126Win = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo4eJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  double VBFHZZ_2e2mu_ev2e2mu_h126Win = sumWeights(filePath_Phantom + "/2mu2e/HZZ4lTree_ZZTo2e2muJJ_SMHContinInterf_H125.6.root",lumi,VBFtag,105.6,140.6);
  VBFZZ_h126Win[2] = VBFZZ_2e2mu_ev4mu_h126Win+VBFZZ_2e2mu_ev4e_h126Win+VBFZZ_2e2mu_ev2e2mu_h126Win;
  VBFHZZ_h126Win[2] = VBFHZZ_2e2mu_ev4mu_h126Win+VBFHZZ_2e2mu_ev4e_h126Win+VBFHZZ_2e2mu_ev2e2mu_h126Win;
  // VBFZZ[2] = (VBFZZ_2e2mu_ev4mu+VBFZZ_2e2mu_ev4e+VBFZZ_2e2mu_ev2e2mu)*(VBFHZZ_yield[2]/(VBFHZZ_h126Win[2]-VBFZZ_h126Win[2]));    
  VBFZZ[2] = VBFZZ_2e2mu_ev4mu+VBFZZ_2e2mu_ev4e+VBFZZ_2e2mu_ev2e2mu;    
  
  cout<<VBFZZ_2e2mu_ev4mu<<" "<<VBFZZ_2e2mu_ev4e<<" "<<VBFZZ_2e2mu_ev2e2mu
      <<" -> ("<<VBFZZ_2e2mu_ev2e2mu<<" + "
      <<VBFZZ_2e2mu_ev4mu+VBFZZ_2e2mu_ev4e<<")*("<<VBFHZZ_yield[2]/(VBFHZZ_h126Win[2]-VBFZZ_h126Win[2])<<")"
      <<" -> "<<VBFZZ[2]<<endl;


  TString schannel[3] = {"4mu","4e","2e2mu"};
  TString ssqrts = (long) sqrts + TString("TeV");
  for (int i=0; i<3; ++i) {
    TString outfile;
    if (VBFtag<2) outfile = "CardFragments/VBFZZRates_" + ssqrts + "_" + schannel[i] + "_" + Form("%d",int(VBFtag)) + ".txt";
    if (VBFtag==2) outfile = "CardFragments/VBFZZRates_" + ssqrts + "_" + schannel[i] + ".txt"; 
    ofstream of(outfile,ios_base::out);
    of << "## rates --- format = chan N lumi ##" << endl
       << "## if lumi is blank, lumi for cards used ##" << endl;
    of << "rate VBFZZ  " << VBFZZ[i] << endl;
    of.close();
    cout << "Output written to: " << outfile << endl;
  }

}


double sumWeights(TString name, double lumi, int selVBF, Double_t minMass, Double_t maxMass){
  TChain* tree = new TChain("SelectedTree");
  tree->Add((TString)(name));

  float MC_weight;
  Short_t NJets;
  double totEvents=0;
  float m4l;

  tree->SetBranchAddress("MC_weight",&MC_weight);  
  tree->SetBranchAddress("NJets30", &NJets);
  tree->SetBranchAddress("ZZMass", &m4l);

  for(int iEvt=0; iEvt<tree->GetEntries(); iEvt++){
    tree->GetEntry(iEvt);
    // cout << "DEBUG sumWeights: " << MC_weight << " / " << NJets << " / " << m4l << endl; 
    if( ( (selVBF == 1 && NJets > 1) || (selVBF == 0 && NJets < 2) || (selVBF==2) ) && ( m4l > minMass && m4l < maxMass ))
      totEvents=totEvents+MC_weight;
  }

  return totEvents*lumi;
}
