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
TString filePath7TeV = "root://lxcms00//data3/2014/HZZ_stat/140328/PRODFSR/";
TString filePath8TeV = "root://lxcms00//data3/2014/HZZ_stat/140211/PRODFSR_8TeV/";

//Local copy
//TString filePath7TeV = "/scratch0/hep/ianderso/CJLST/140328/PRODFSR/";
//TString filePath8TeV = "/scratch0/hep/ianderso/CJLST/140211/PRODFSR_8TeV/"; 


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
<<<<<<< HEAD
const int nVBFPoints7TeV_HM = 10; 
int VBFmasses7TeV_HM[nVBFPoints7TeV_HM]   = {400,450,500,550,600,650,700,800,900,1000};
double mHVBFVal7TeV_HM[nVBFPoints7TeV_HM] = {400,450,500,550,600,650,700,800,900,1000};

const int nVBFPoints8TeV_HM = 10;
int VBFmasses8TeV_HM[nVBFPoints8TeV_HM]   = {400,450,500,550,600,650,700,800,900,1000};
double mHVBFVal8TeV_HM[nVBFPoints8TeV_HM] = {400,450,500,550,600,650,700,800,900,1000};


const int nVBFPoints7TeV = 32; //FIXME: 950 GeV sample is not there since high mass weights are not available
int VBFmasses7TeV[nVBFPoints7TeV]   = {115,120,125,130,140,150,160,170,180,190,200,210,220,230,250,275,300,325,350,375,400,425,450,475,500,575,600,650,700,800,900,/*950,*/1000};
double mHVBFVal7TeV[nVBFPoints7TeV] = {115,120,125,130,140,150,160,170,180,190,200,210,220,230,250,275,300,325,350,375,400,425,450,475,500,575,600,650,700,800,900,/*950,*/1000};
=======
const int nVBFPoints7TeV = 25; //FIXME: 950 GeV sample is not there since high mass weights are not available
int VBFmasses7TeV[nVBFPoints7TeV]   = {115,120,125,130,140,150,160,170,180,190,200,225,275,300,350,400,450,500,550,600,650,700,800,900,/*950,*/1000};
double mHVBFVal7TeV[nVBFPoints7TeV] = {115,120,125,130,140,150,160,170,180,190,200,225,275,300,350,400,450,500,550,600,650,700,800,900,/*950,*/1000};
>>>>>>> b6c3a520c9fed687e1091566997144e9ab068969

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

/*
//qqZZ
const int nqqZZPoints = 331;
int qqZZmasses[nqqZZPoints] = {100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202
,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320,3
22,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440
,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558,5
60,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,656,658,660,662,664,666,668,670,672,674,676,678
,680,682,684,686,688,690,692,694,696,698,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000};
double mHqqZZVal[nqqZZPoints] = {100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,2
02,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320
,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,4
40,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558
,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,656,658,660,662,664,666,668,670,672,674,676,6
78,680,682,684,686,688,690,692,694,696,698,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000};

//ggZZ
const int nggZZPoints = 91;
int ggZZmasses[nggZZPoints] = {100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610
,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000};
double mHggZZVal[nggZZPoints] = {100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,6
10,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000};*/


//Values from Roberto, for sum of 7+8 TeV
double filter_eff_WH_8TeV  = 0.010275; // (47888+98897)/(4761904+9523809)
double filter_eff_ZH_8TeV  = 0.028616; // (51246+58135)/(1785714+2036645)
double filter_eff_ttH_8TeV = 0.029688;   // (30306+49210)/(1013513+1664865)
