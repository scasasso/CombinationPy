//----------> SET INPUT VARIABLES HERE

// Input trees

// ICHEP inputs
// TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/240612/PRODFSR";
// TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/240612/PRODFSR_8TeV";

// With LD fixes (phi1, Z1/Z2)
// TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812/PRODFSR";
// TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812/PRODFSR_8TeV";

// With LD fixes (phi1, Z1/Z2), but use rescaled 7-TeV gZZ for 8 TeV
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled/PRODFSR_8TeV/";

// With high mass reweights
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_HighMassReweight/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_HighMassReweight/PRODFSR_8TeV/";

// Experimental
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled_HighMassReweight/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled_HighMassReweight/PRODFSR_8TeV/";

//171012, HCP unblinding
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/171012/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/171012/PRODFSR_8TeV/";

//181012, HCP unblinding
// TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/191012/PRODFSR/";
// TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/191012/PRODFSR_8TeV/";

//Final set, with ICHEP data/MC SFs
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012_ichep/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012_ichep/PRODFSR_8TeV/";

//Final set, HCP
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012/PRODFSR_8TeV/";

//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/";

//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130120/v2/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130120/v2/PRODFSR_8TeV/";

// Moriond, final set
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130205/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130205/PRODFSR_8TeV/";

// Legacy paper
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702b/PRODFSR/";	 
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702b/PRODFSR_8TeV/";

//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702b/PRODFSR/"; //FIXME: statistical tree not yet ready, take those from previous production
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130715/PRODFSR_8TeV/";

/* TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130720d/PRODFSR/"; */
/* TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130720d/PRODFSR_8TeV/"; */

// High mass paper reweighted samples (needed to extract resolution function or ggH)
/* TString filePath7TeV = "root://lxcms02///data/Higgs/rootuplesOut/140311_HMPaper/C0.2/PRODFSR/";	  */
/* TString filePath8TeV = "root://lxcms02///data/Higgs/rootuplesOut/140311_HMPaper/C0.2/PRODFSR_8TeV/";	  */

/* TString filePath7TeV = "root://lxcms02///data/Higgs/rootuplesOut/140311_HMPaper/C1.0/PRODFSR/";	  */
/* TString filePath8TeV = "root://lxcms02///data/Higgs/rootuplesOut/140311_HMPaper/C1.0/PRODFSR_8TeV/";	  */


// High mass paper reweighted samples (for signal shapes)
/* TString filePath7TeV = "/afs/cern.ch/work/s/scasasso/private/H4l_HighMass_Git/CMSSW_5_3_9/src/ZZAnalysis/AnalysisStep/test/Macros/140515_LOPlusInt/C02/PRODFSR/"; */
/* TString filePath8TeV = "/afs/cern.ch/work/s/scasasso/private/H4l_HighMass_Git/CMSSW_5_3_9/src/ZZAnalysis/AnalysisStep/test/Macros/140515_LOPlusInt/C02/PRODFSR_8TeV/"; */

/* TString filePath7TeV = "/afs/cern.ch/work/s/scasasso/private/H4l_HighMass_Git/CMSSW_5_3_9/src/ZZAnalysis/AnalysisStep/test/Macros/140515_LOPlusInt/C10/PRODFSR/"; */
/* TString filePath8TeV = "/afs/cern.ch/work/s/scasasso/private/H4l_HighMass_Git/CMSSW_5_3_9/src/ZZAnalysis/AnalysisStep/test/Macros/140515_LOPlusInt/C10/PRODFSR_8TeV/"; */

//TString filePath7TeV = "root://lxcms02///data/Higgs/rootuplesOut/140311_HMPaper/C1.0/PRODFSR/";	 
//TString filePath8TeV = "root://lxcms02///data/Higgs/rootuplesOut/140311_HMPaper/C1.0/PRODFSR_8TeV/";	

//Full iteration? Check later
//TString filePath7TeV = "root://lxcms00//data3/2014/HZZ_stat/140519/PRODFSR/";
//TString filePath8TeV = "root://lxcms00//data3/2014/HZZ_stat/140519/PRODFSR_8TeV/";

//Local copy
TString filePath7TeV = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR/";
TString filePath8TeV = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV/"; 


// Luminosity, as float and as string to be used in file names, etc.
double lumi7TeV = 5.051;
double lumi8TeV = 19.712;
TString lumistr7TeV = "5.051";
TString lumistr8TeV = "19.712";


// Location of output root files containing data events
TString DataRootFilePath = "../CreateDatacards/CMSdata/"; 

//--------------------
// The number and values of mass points for which you have signal trees, for 7 and 8 TeV

//--- ggH

// Old powheg ggH samples
const int nPoints7TeV = 34;
int masses7TeV[nPoints7TeV]   = {120,124,125,126,130,140,150,160,170,180,190,200,210,220,250,275,300,325,350,400,425,450,475,525,550,575,600,650,700,750,800,900,950,1000};
double mHVal7TeV[nPoints7TeV] = {120,124,125,126,130,140,150,160,170,180,190,200,210,220,250,275,300,325,350,400,425,450,475,525,550,575,600,650,700,750,800,900,950,1000};

const int nPoints8TeV = 49;
int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};
double mHVal8TeV[nPoints8TeV] = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};

// Full list of new powheg samples: you must use "powheg15_jhuGenV3" samples for m<=200, "powheg15" ones above. 
const int nPoints7TeV_p15 = 34;
int masses7TeV_p15[nPoints7TeV_p15]   = {115,120,122,124,125,126,128,130,135,140,145,150,160,170,175,180,185,190,200,225,250,275,300,350,400,450,500,550,600,650,700,800,900,1000};
double mHVal7TeV_p15[nPoints7TeV_p15] = {115,120,122,124,125,126,128,130,135,140,145,150,160,170,175,180,185,190,200,225,250,275,300,350,400,450,500,550,600,650,700,800,900,1000};

const int nPoints8TeV_p15 = 34;
int masses8TeV_p15[nPoints8TeV_p15]   = {115,120,122,124,125,126,128,130,135,140,145,150,160,170,175,180,185,190,200,225,250,275,300,350,400,450,500,550,600,650,700,800,900,1000};
double mHVal8TeV_p15[nPoints8TeV_p15] = {115,120,122,124,125,126,128,130,135,140,145,150,160,170,175,180,185,190,200,225,250,275,300,350,400,450,500,550,600,650,700,800,900,1000};

// List of new powheg samples for high mass (mH>400 GeV) studies
const int nPoints7TeV_p15_HM = 10;
int masses7TeV_p15_HM[nPoints7TeV_p15_HM]   = {400,450,500,550,600,650,700,800,900,1000};
double mHVal7TeV_p15_HM[nPoints7TeV_p15_HM] = {400,450,500,550,600,650,700,800,900,1000};

const int nPoints8TeV_p15_HM = 10;
int masses8TeV_p15_HM[nPoints8TeV_p15_HM]   = {400,450,500,550,600,650,700,800,900,1000};
double mHVal8TeV_p15_HM[nPoints8TeV_p15_HM] = {400,450,500,550,600,650,700,800,900,1000};

//MINLO samples
const int nPoints7TeV_MINLO = 15;
int masses7TeV_MINLO[nPoints7TeV_MINLO]   = {125,126,190,200,250,300,350,400,450,500,600,750,900,950,1000};
double mHVal7TeV_MINLO[nPoints7TeV_MINLO] = {125,126,190,200,250,300,350,400,450,500,600,750,900,950,1000};

const int nPoints8TeV_MINLO = 35;
int masses8TeV_MINLO[nPoints8TeV_MINLO]   = {90,95,100,105,110,115,120,124,125,126,130,135,140,145,150,155,170,180,190,200,250,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};
double mHVal8TeV_MINLO[nPoints8TeV_MINLO] = {90,95,100,105,110,115,120,124,125,126,130,135,140,145,150,155,170,180,190,200,250,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};


//List of "powheg15" ggH samples alone (available for high mass only +125,126); just for debugging purposes
// const int nPoints7TeV_p15 = 12;
// int masses7TeV_p15[nPoints7TeV_p15]   = {125,126,400,450,500,550,600,650,700,800,900,1000};
// double mHVal7TeV_p15[nPoints7TeV_p15] = {125,126,400,450,500,550,600,650,700,800,900,1000};

// const int nPoints8TeV_p15 = 17;
// int masses8TeV_p15[nPoints8TeV_p15]   = {125,126,225,250,275,300,350,400,450,500,550,600,650,700,800,900,1000};
// double mHVal8TeV_p15[nPoints8TeV_p15] = {125,126,225,250,275,300,350,400,450,500,550,600,650,700,800,900,1000};


//--- VBF
const int nVBFPoints7TeV_HM = 10; 
int VBFmasses7TeV_HM[nVBFPoints7TeV_HM]   = {400,450,500,550,600,650,700,800,900,1000};
double mHVBFVal7TeV_HM[nVBFPoints7TeV_HM] = {400,450,500,550,600,650,700,800,900,1000};

const int nVBFPoints8TeV_HM = 10;
int VBFmasses8TeV_HM[nVBFPoints8TeV_HM]   = {400,450,500,550,600,650,700,800,900,1000};
double mHVBFVal8TeV_HM[nVBFPoints8TeV_HM] = {400,450,500,550,600,650,700,800,900,1000};


const int nVBFPoints7TeV = 25; //FIXME: 950 GeV sample is not there since high mass weights are not available
int VBFmasses7TeV[nVBFPoints7TeV]   = {115,120,125,130,140,150,160,170,180,190,200,225,275,300,350,400,450,500,550,600,650,700,800,900,1000};
double mHVBFVal7TeV[nVBFPoints7TeV] = {115,120,125,130,140,150,160,170,180,190,200,225,275,300,350,400,450,500,550,600,650,700,800,900,1000};

const int nVBFPoints8TeV = 42;
int VBFmasses8TeV[nVBFPoints8TeV]   = {116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,225,250,275,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};
double mHVBFVal8TeV[nVBFPoints8TeV] = {116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,225,250,275,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};

//--- WH/ZH/ttH
const int nVHPoints7TeV = 11;
int VHmasses7TeV[nVHPoints7TeV]   = {110,115,120,125,126,130,140,150,160,180,200}; // 
double mHVHVal7TeV[nVHPoints7TeV] = {110,115,120,125,126,130,140,150,160,180,200}; // 

const int nVHPoints8TeV = 11;
int VHmasses8TeV[nVHPoints8TeV]   = {110,115,120,125,126,130,140,150,160,180,200};
double mHVHVal8TeV[nVHPoints8TeV] = {110,115,120,125,126,130,140,150,160,180,200};

double filter_eff_WH_8TeV  = 0.010275; // (47888+98897)/(4761904+9523809)
double filter_eff_ZH_8TeV  = 0.028616; // (51246+58135)/(1785714+2036645)
double filter_eff_ttH_8TeV = 0.029688;   // (30306+49210)/(1013513+1664865)

