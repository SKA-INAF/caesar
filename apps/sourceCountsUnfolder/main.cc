// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************

#include "SpectrumUtils.h"
#include "ForwardFolder.h"

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>
#include <MathUtils.h>
#include <StatsUtils.h>


//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace Caesar;

void Usage(char* exeName)
{
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input=[INPUT_FILE] \t Input ROOT file with source count histogram stored to be unfolded"<<endl;	
	cout<<"-H, --histoname=[INPUT_FILE] \t Name of histogram with spectrum data to be unfolded (default=SourceCounts)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (ROOT format) where to store outputs (default=output.root)"<<endl;
	cout<<"-g, --gamma1=[GAMMA1] \t 1st spectral index of spectrum model (default=1)"<<endl;
	cout<<"-G, --gamma2=[GAMMA2] \t 2nd spectral index of spectrum model (default=2.7)"<<endl;
	cout<<"-l, --gamma3=[GAMMA3] \t 3rd spectral index of spectrum model (default=3.5)"<<endl;
	cout<<"-c, --lgS-break=[LGS_BREAK] \t Spectrum break par (default=-2.88)"<<endl;
	cout<<"-C, --lgS-cutoff=[LGS_CUTOFF] \t Spectrum cutoff par (default=-0.67)"<<endl;
	cout<<"-E, --effmodel=[EFFICIENCY_MODEL] \t Completeness model (const,sigmoid)"<<endl;
	cout<<"-e, --effpars=[EFFICIENCY_PARS] \t Completeness model pars"<<endl;
	cout<<"-B, --biasmodel=[BIAS_MODEL] \t Bias model (const,exp,stepexp)"<<endl;
	cout<<"-b, --biaspars=[BIAS_PARS] \t Flux bias model pars"<<endl;
	cout<<"-R, --resomodel=[RESO_MODEL] \t Reso model (const,exp)"<<endl;
	cout<<"-r, --resopars=[RESO_PARS] \t Flux reso model pars"<<endl;
	cout<<"-S, --smodel=[SPECTRUM_MODEL] \t Spectrum model used for response matrix creation (flat,powerlaw,brokenpowerlaw,pol3)"<<endl;
	cout<<"-s, --smodel-fit=[SPECTRUM_MODEL] \t Spectrum model used as fit model in unfolding (flat,powerlaw,brokenpowerlaw,pol3)"<<endl;
	cout<<"-d, --spolpars=[SPECTRUM_POL_PARS] \t Flux spectrum polynomial model pars"<<endl;
	cout<<"-j, --useResoBiasCorrection \t If enabled, include resolution bias correction in response matrix"<<endl;
	cout<<"-u, --uncertainties \t Compute stats & syst uncertainties in unfolded flux"<<endl;
	cout<<"-t, --struebins=[STRUE_BINS] \t True flux bins (set equal to rec bins if not provided)"<<endl;
	cout<<"-w, --useErrorsInChi2 \t If enabled, include bin errors in chi2 definition, otherwise set all errors to 1"<<endl;
	cout<<"-L, --likelihoodFit \t If enabled, minimize -log likelihood"<<endl;
	cout<<"-m, --useFitRange \t If enabled, limit fit to specified flux density range (default=false)"<<endl;
	cout<<"-f, --lgSfit-min=[LGSFIT_MIN] \t Minimum flux density range used in fitting (default=-4)"<<endl;
	cout<<"-F, --lgSfit-max=[LGSFIT_MAX] \t Maximum flux density range used in fitting (default=1.5)"<<endl;
	cout<<"-n, --freq=[FREQ] \t Data frequency in GHz used to convert source count literature fit models (default=1.4 GHz)"<<endl;	
	cout<<"-N, --sindex=[SPECTRAL_INDEX] \t Data frequency spectral index used to convert source count literature fit models (default=-0.9)"<<endl;
	cout<<"-A, --area=[AREA] \t Normalization area in deg squared (default=37.7)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;

}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "histoname", required_argument, 0, 'H' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", required_argument, 0, 'o' },
	{ "gamma1", required_argument, 0, 'g' },
	{ "gamma2", required_argument, 0, 'G' },
	{ "gamma3", required_argument, 0, 'l' },
	{ "lgS-break", required_argument, 0, 'c' },
	{ "lgS-cutoff", required_argument, 0, 'C' },
	{ "effpars", required_argument, 0, 'e' },	
	{ "biaspars", required_argument, 0, 'b' },
	{ "resopars", required_argument, 0, 'r' },
	{ "effmodel", required_argument, 0, 'E' },	
	{ "biasmodel", required_argument, 0, 'B' },
	{ "resomodel", required_argument, 0, 'R' },
	{ "smodel", required_argument, 0, 'S' },
	{ "smodel-fit", required_argument, 0, 's' },
	{ "spolpars", required_argument, 0, 'd' },
	{ "uncertainties", no_argument, 0, 'u' },	
	{ "struebins", required_argument, 0, 't' },
	{ "useResoBiasCorrection", no_argument, 0, 'j' },
	{ "useErrorsInChi2", no_argument, 0, 'w' },
	{ "likelihoodFit", no_argument, 0, 'L' },	
	{ "useFitRange", no_argument, 0, 'm' },	
	{ "lgSfit-min", required_argument, 0, 'f' },
	{ "lgSfit-max", required_argument, 0, 'F' },
	{ "freq", required_argument, 0, 'n' },
	{ "sindex", required_argument, 0, 'N' },	
	{ "area", required_argument, 0, 'A' },	
	{ "interactive", no_argument, 0, 'I' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string histoName= "";
std::string outputFileName= "output.root";
int verbosity= 4;//INFO level
bool gRunInteractively= false;
double gSpectrumModelIndex1= 1;
double gSpectrumModelIndex2= 2.7;
double gSpectrumModelIndex3= 3.5;
double gSpectrumModelLgSBreak= -2.88;
double gSpectrumModelLgSBreak_min= -3.2;
double gSpectrumModelLgSBreak_max= -2.2;
double gSpectrumModelLgSCutoff= -0.67;
double gSpectrumModelLgSCutoff_min= -1;
double gSpectrumModelLgSCutoff_max= -0.3;
double gSpectrumModelCutoffSlope= 0.5;
double gSpectrumModelCutoffSlope_min= -10000;
double gSpectrumModelCutoffSlope_max= 10000;
std::vector<double> gSpectrumModelPolPars {224.546,188.565,53.5395,5.12536};
//std::vector<double> gEfficiencyPars {9.86885e-01,-2.85922e+00,4.14620e-01};
std::vector<double> gEfficiencyPars {9.86970e-01,-2.87361e+00,3.94529e-01};
std::vector<double> gEfficiencyParErrors {1.24084e-03,2.39855e-03,2.98135e-03};
std::vector<std::vector<double>> gEfficiencyParCovMatrix 
{
	{1.54e-06,1.091e-06,1.079e-06},
  {1.091e-06,5.753e-06,4.599e-06}, 
  {1.079e-06,4.599e-06,8.888e-06}
};
std::vector<double> gFluxBiasPars {1.43202e-06,3.64812e+00,-2.70949e-04,-3};
std::vector<double> gFluxBiasParErrors {3.12949e-07,6.99724e-02,1.07134e-04,0};
std::vector<std::vector<double>> gFluxBiasParCovMatrix 
{
	{9.794e-14,-2.177e-08,-3.524e-12,0}, 
  {-2.177e-08,0.004896,7.324e-07,0}, 
  {-3.524e-12,7.324e-07,1.148e-08,0},
	{0,0,0,0}
};

std::vector<double> gFluxResoPars {4.48925e-04,1.81314e+00,-3};
std::vector<double> gFluxResoParErrors {2.12782e-05,1.74581e-02,0};
std::vector<std::vector<double>> gFluxResoParCovMatrix 
{
	{4.528e-10,-3.631e-07,0},
  {-3.631e-07,0.0003048,0},
	{0,0,0}
};

int gBiasModel= eCONST_RESO;
int gResoModel= eCONST_BIAS;
int gEfficiencyModel= eCONST_EFF;
int gSpectrumModel= eFlat;
//int gSpectrumModel_fit= eBrokenPowerLaws;
int gSpectrumModel_fit= ePol3;
SpectrumPars* gSpectrumPars= 0;
double gConstEfficiency= 1;
double gConstBias= 0;
double gConstReso= 0;
bool gUseFitRange= false;
//double gLgSMin_fit= -3.;
double gLgSMin_fit= -4;
double gLgSMax_fit= 1.5;
bool gUseErrorsInChi2= false;
bool gComputeUncertainties= false;
bool gUseLikelihoodFit= false;
int gNRandSamples= 1000;
double gArea= 37.7;//in deg^2
double gNormSpectralIndex= 2.5;
double gDataFrequency= 1.4;//GHz
double gDataSpectralIndex= -0.9;
double gBmaj= 24;//arcsec
double gBmin= 20;//arcsec
double gNoiseRms= 300.e-3;//mJy
double gSThr= 1.5;//mJy
double gResolvedSourceThrP0= 1.08;
double gResolvedSourceThrP1= 2.03;
double gResolvedSourceThrSlope= 1;
double gPhiMedianSlope= 0.30;
double gPhiMedianSbreak= 1;
double gPhiMedianScaleFactor= 1;
bool gUseSourceCountsResoBiasCorr= true;
double gSourceCountsResoBiasCorrErr= 0.1;
double gExpSourceCountDataStatErr= 0.2;//taken from Windhorst+90 (section 3)

//Globar vars
TFile* inputFile= 0;
TFile* outputFile= 0;
TApplication* app= 0;
TH1D* gSourceCounts= 0;
std::vector<double> gRecBins;
std::vector<double> gTrueBins;
TF1* gTrueFitModelFcn= 0;
TH2D* gResponseMatrix= 0;
TF2* gResponseModelFcn= 0;
TH2D* gResponseMatrix_fluctuated= 0;
TF1* gResoModelFcn= 0;
TF1* gBiasModelFcn= 0;
ForwardFolder* gForwardFolder= 0;
TH1D* gUnfoldedSpectrum= 0;
TH1D* gUnfoldedSpectrum_totErrors= 0;
TH1D* gFFSpectrum= 0;
TH1D* gDiffSourceCounts= 0;
TH1D* gUnfoldedDiffSpectrum= 0;
TH1D* gUnfoldedDiffSpectrum_totErrors= 0;
TH1D* gUnfoldedDiffSpectrum_resoBiasCorrected= 0;
TH1D* gFFDiffSpectrum= 0;
TF1* gSourceCountsResoBiasCorrFcn= 0;
TH1D* gSourceCountsFitResidualHisto= 0;
TH2D* gRandSigmaResponseMatrix= 0;
TH2D* gRandMeanResponseMatrix= 0;
std::vector<TH2D*> gRandResponseMatrices;
TGraph* gExpDataFitAreaGraph= 0;
TGraph* gExpDataFitAreaGraph_min= 0;
TGraph* gExpDataFitAreaGraph_max= 0;
TGraph* gExpDataFitAreaGraph_mean= 0;
TF1* gExpDataFitFcn_Katgert= 0;
TF1* gExpDataFitFcn_Hopkins= 0;

std::vector<double> gSourceCountsExpDataFitPars_Hopkins {0.859,0.508,0.376,-0.049,-0.121,0.057,-0.008};
std::vector<double> gSourceCountsExpDataFitPars_Katgert {0.908,0.619,0.190,-0.075};

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadData(std::string filename);
int SetResponseMatrix();
int ComputeResponseMatrixErrors();
int ComputeSourceCountsModels();
int UnfoldData();
void Save();
int Init();
void PrintSourceCountTable();
void Draw();
void SetStyle();
void ClearData();
double Deg2ToSr(double deg2)
{
	double f= pow(TMath::Pi()/180.,2);
	double sr= deg2*f;
	return sr;
}

int main(int argc, char *argv[])
{
	//================================
	//== PARSE CMD LINE OPTIONS
	//================================
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Failed to parse command line options!");
		#endif
		return -1;
	}
	
	//=======================
	//== INITIALIZE
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing data...");
	#endif
	if(Init()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to initialize data!");
		#endif
		ClearData();
		return -1;
	}


	//=======================
	//== READ DATA
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading input data from file "<<fileName<<" ...");
	#endif
	if(ReadData(fileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of input data failed!");
		#endif
		ClearData();
		return -1;
	}

	//=======================
	//== SET RESPONSE MATRIX
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting response matrix ...");
	#endif
	if(SetResponseMatrix()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Setting of response matrix failed!");
		#endif
		ClearData();
		return -1;
	}
	
	//==========================
	//== UNFOLD
	//==========================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Unfolding source counts ...");
	#endif
	if(UnfoldData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Source counts unfolding failed!");
		#endif
		ClearData();
		return -1;
	}

	//==================================
	//== COMPUTE SOURCE COUNTS MODELS
	//==================================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing source counts models ...");
	#endif	
	if(ComputeSourceCountsModels()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Computing of source counts model failed!");
		#endif
		ClearData();
		return -1;
	}


	//=======================
	//== DRAW PLOTS
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Drawing data ...");
	#endif
	Draw();

	//=======================
	//== PRINT TABLE
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Print source counts table ...");
	#endif
	PrintSourceCountTable();

	//=======================
	//== SAVE
	//=======================
	// - Save data to file or run interatively
	if(app && gRunInteractively) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("End application run");	
		#endif
		app->Run();
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving data to file ...");
		#endif
		Save();
	}

	return 0;

}//close main



int ParseOptions(int argc, char *argv[])
{
	//## Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Parse options
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hi:H:o:v:Ig:G:l:c:C:e:b:r:E:B:R:S:s:d:t:uwLn:N:jmf:F:A:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
    	case 'i':	
			{
				fileName= std::string(optarg);	
				break;	
			}
			case 'H':	
			{
				histoName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'I':	
			{
				gRunInteractively= true;
				break;	
			}
			case 'g':	
			{
				gSpectrumModelIndex1= atof(optarg);
				break;	
			}
			case 'G':	
			{
				gSpectrumModelIndex2= atof(optarg);
				break;	
			}
			case 'l':	
			{
				gSpectrumModelIndex3= atof(optarg);
				break;	
			}
			case 'c':	
			{
				gSpectrumModelLgSBreak= atof(optarg);
				break;	
			}
			case 'C':	
			{
				gSpectrumModelLgSCutoff= atof(optarg);
				break;
			}
			case 'e':	
			{
				std::string effPars_str= std::string(optarg);
				std::vector<std::string> effPars_strs= CodeUtils::SplitStringOnPattern(effPars_str,',');
				gEfficiencyPars.clear();
				for(size_t i=0;i<effPars_strs.size();i++){					
					double effPar= atof(effPars_strs[i].c_str());
					cout<<"effPar_str="<<effPars_strs[i]<<", effPar="<<effPar<<endl;
					gEfficiencyPars.push_back(effPar);
				}
				break;	
			}
			case 'b':	
			{
				std::string biasPars_str= std::string(optarg);
				std::vector<std::string> biasPars_strs= CodeUtils::SplitStringOnPattern(biasPars_str,',');
				gFluxBiasPars.clear();
				for(size_t i=0;i<biasPars_strs.size();i++){
					double biasPar= atof(biasPars_strs[i].c_str());
					cout<<"biasPar_str="<<biasPars_strs[i]<<", biasPar="<<biasPar<<endl;
					gFluxBiasPars.push_back(biasPar);
				}
				break;	
			}
			case 'r':	
			{
				std::string resoPars_str= std::string(optarg);
				std::vector<std::string> resoPars_strs= CodeUtils::SplitStringOnPattern(resoPars_str,',');
				gFluxResoPars.clear();
				for(size_t i=0;i<resoPars_strs.size();i++){
					double resoPar= atof(resoPars_strs[i].c_str());
					cout<<"resoPar_str="<<resoPars_strs[i]<<", resoPar="<<resoPar<<endl;
					gFluxResoPars.push_back(resoPar);
				}
				break;	
			}
			case 'E':	
			{
				std::string effModel_str= std::string(optarg);
				if(effModel_str=="const") gEfficiencyModel= eCONST_EFF;
				else if(effModel_str=="sigmoid") gEfficiencyModel= eSIGMOID_EFF;
				else gEfficiencyModel= eCONST_EFF;
				break;	
			}
			case 'B':	
			{
				std::string biasModel_str= std::string(optarg);
				if(biasModel_str=="const") gBiasModel= eCONST_BIAS;
				else if(biasModel_str=="exp") gBiasModel= eEXP_BIAS;
				else if(biasModel_str=="expstep") gBiasModel= eEXP_STEP_BIAS;
				else gBiasModel= eCONST_EFF;
				break;	
			}
			case 'R':	
			{
				std::string resoModel_str= std::string(optarg);
				if(resoModel_str=="const") gResoModel= eCONST_RESO;
				else if(resoModel_str=="exp") gResoModel= eEXP_RESO;
				else if(resoModel_str=="expstep") gResoModel= eEXP_STEP_RESO;
				else gResoModel= eCONST_EFF;
				break;	
			}
			case 'S':	
			{
				std::string spectrumModel_str= std::string(optarg);
				if(spectrumModel_str=="flat") gSpectrumModel= eFlat;
				else if(spectrumModel_str=="powerlaw") gSpectrumModel= ePowerLaw;
				else if(spectrumModel_str=="brokenpowerlaw") gSpectrumModel= eBrokenPowerLaws;
				else if(spectrumModel_str=="2brokenpowerlaw") gSpectrumModel= eTwoBrokenPowerLaws;
				else if(spectrumModel_str=="smoothbrokenpowerlaw") gSpectrumModel= eSmoothBrokenPowerLaws;	
				else if(spectrumModel_str=="cutoffpowerlaw") gSpectrumModel= ePowerLawWithCutoff;
				else if(spectrumModel_str=="pol3") gSpectrumModel= ePol3;
				else if(spectrumModel_str=="pol6") gSpectrumModel= ePol6;
				else gSpectrumModel= eFlat;
				break;	
			}
			case 's':	
			{
				std::string spectrumModel_str= std::string(optarg);
				if(spectrumModel_str=="flat") gSpectrumModel_fit= eFlat;
				else if(spectrumModel_str=="powerlaw") gSpectrumModel_fit= ePowerLaw;
				else if(spectrumModel_str=="brokenpowerlaw") gSpectrumModel_fit= eBrokenPowerLaws;	
				else if(spectrumModel_str=="2brokenpowerlaw") gSpectrumModel_fit= eTwoBrokenPowerLaws;
				else if(spectrumModel_str=="smoothbrokenpowerlaw") gSpectrumModel_fit= eSmoothBrokenPowerLaws;
				else if(spectrumModel_str=="cutoffpowerlaw") gSpectrumModel_fit= ePowerLawWithCutoff;
				else if(spectrumModel_str=="pol3") gSpectrumModel_fit= ePol3;
				else if(spectrumModel_str=="pol6") gSpectrumModel_fit= ePol6;
				else gSpectrumModel_fit= ePol3;
				break;	
			}
			case 'd':	
			{
				std::string polPars_str= std::string(optarg);
				std::vector<std::string> polPars_strs= CodeUtils::SplitStringOnPattern(polPars_str,',');
				gSpectrumModelPolPars.clear();
				for(size_t i=0;i<polPars_strs.size();i++){
					double polPar= atof(polPars_strs[i].c_str());
					gSpectrumModelPolPars.push_back(polPar);
				}
				break;	
			}
			case 't':	
			{
				std::string STrueBins_str= std::string(optarg);
				std::vector<std::string> STrueBins_strs= CodeUtils::SplitStringOnPattern(STrueBins_str,',');
				gTrueBins.clear();
				for(size_t i=0;i<STrueBins_strs.size();i++){					
					double STrueBin= atof(STrueBins_strs[i].c_str());
					gTrueBins.push_back(STrueBin);
				}
				break;	
			}
			case 'u':	
			{
				gComputeUncertainties= true;
				break;	
			}
			case 'j':
			{
				gUseSourceCountsResoBiasCorr= true;
				break;
			}
			case 'w':	
			{
				gUseErrorsInChi2= true;
				break;	
			}
			case 'L':	
			{
				gUseLikelihoodFit= true;
				break;
			}

			case 'm':
			{
				gUseFitRange= true;
				break;
			}
			case 'f':
			{
				gLgSMin_fit= atof(optarg);
				break;
			}
			case 'F':
			{
				gLgSMax_fit= atof(optarg);
				break;
			}
			case 'n':	
			{	
				gDataFrequency= atof(optarg);
				break;
			}
			case 'N':	
			{	
				gDataSpectralIndex= atof(optarg);
				break;
			}
			case 'A':	
			{	
				gArea= atof(optarg);
				break;
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
 
	
	//=======================
	//== INIT LOGGER
	//=======================
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	//=======================
	//== CHECK ARGS 
	//=======================
	//Check input file name
	if(fileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty input file name given!");
		#endif
		return -1;
	}
	
	//Check histo name
	if(histoName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty histo name given!");
		#endif
		return -1;
	}

	//Check output file name
	if(outputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output file name given!");
		#endif
		return -1;
	}

	//Check bias pars
	if(gFluxBiasPars.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty bias pars given!");
		#endif
		return -1;
	}
	if(gBiasModel==eEXP_BIAS && gFluxBiasPars.size()!=3){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid bias pars given (3 expected)!");
		#endif
		return -1;
	}
	if(gBiasModel==eEXP_STEP_BIAS && gFluxBiasPars.size()!=4){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid bias pars given (4 expected)!");
		#endif
		return -1;
	}
	if(gBiasModel==eCONST_BIAS && gFluxBiasPars.size()!=1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid bias pars given (1 expected)!");
		#endif
		return -1;
	}

	//Check reso pars
	if(gResoModel==eEXP_RESO && gFluxResoPars.size()!=2){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid reso pars given (2 pars required)!");
		#endif
		return -1;
	}
	if(gResoModel==eEXP_STEP_RESO && gFluxResoPars.size()!=3){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid reso pars given (3 pars required)!");
		#endif
		return -1;
	}
	if(gResoModel==eCONST_RESO && gFluxResoPars.size()!=1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid reso pars given (1 par required)!");
		#endif
		return -1;
	}

	//Check efficiency pars
	if(gEfficiencyModel==eSIGMOID_EFF && gEfficiencyPars.size()!=3){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid efficiency pars given (3 pars required)!");
		#endif
		return -1;
	}
	if(gEfficiencyModel==eCONST_EFF && gEfficiencyPars.size()!=1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid efficiency pars given (1 pars required)!");
		#endif
		return -1;
	}

	cout<<"DEBUG: Bias pars {";
	for(size_t i=0;i<gFluxBiasPars.size()-1;i++){
		cout<<gFluxBiasPars[i]<<",";
	}
	cout<<gFluxBiasPars[gFluxBiasPars.size()-1]<<"}"<<endl;

	cout<<"DEBUG: Reso pars {";
	for(size_t i=0;i<gFluxResoPars.size()-1;i++){
		cout<<gFluxResoPars[i]<<",";
	}
	cout<<gFluxResoPars[gFluxResoPars.size()-1]<<"}"<<endl;

	cout<<"DEBUG: Efficiency pars {";
	for(size_t i=0;i<gEfficiencyPars.size()-1;i++){
		cout<<gEfficiencyPars[i]<<",";
	}
	cout<<gEfficiencyPars[gEfficiencyPars.size()-1]<<"}"<<endl;

	//=======================
	//== INIT ROOT APP
	//=======================
	if(!app && gRunInteractively){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Initializing ROOT app...");
		#endif
		app= new TApplication("App",&argc,argv);
	}
	
	return 0;

}//close ParseOptions()

int Init()
{
	//- Set draw & graphics style
	SetStyle();

	//- Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//- Initialize forward folder
	if(!gForwardFolder) gForwardFolder= new ForwardFolder;	

	return 0;

}//close Init()


void ClearData()
{
	//Delete data
	if(gSourceCounts){
		delete gSourceCounts;
		gSourceCounts= 0;
	}

	if(gForwardFolder){
		delete gForwardFolder;
		gForwardFolder= 0;
	}

	//Delete app (if created)
	if(app){
		delete app;
		app= 0;
	}

	//Close file
	if(outputFile && outputFile->IsOpen()){
		outputFile->Close();
	}

}//close ClearData()

int UnfoldData()
{
	//====================================================
	//==         SET SPECTRUM PARS
	//=====================================================
	// - Power law spectrum pars
	double gamma1= gSpectrumModelIndex1;
	double gamma2= gSpectrumModelIndex2;
	double gamma3= gSpectrumModelIndex3;
	double gamma_min= 0;
	double gamma_max= 10;
	double gamma1_min= -2;
	double gamma1_max= -1;
	double gamma2_min= 0;
	double gamma2_max= 1;
	double gamma3_min= 0.5;
	double gamma3_max= 1.5;
		
	double lgS_break= gSpectrumModelLgSBreak;
	double lgS_break_min= gSpectrumModelLgSBreak_min;
	double lgS_break_max= gSpectrumModelLgSBreak_max;
	double lgS_cutoff= gSpectrumModelLgSCutoff;
	double lgS_cutoff_min= gSpectrumModelLgSCutoff_min;
	double lgS_cutoff_max= gSpectrumModelLgSCutoff_max;
	double lgS_min= -5;
	double lgS_max= 10;

	double Wc= gSpectrumModelCutoffSlope;
	double Wc_min= gSpectrumModelCutoffSlope_min;
	double Wc_max= gSpectrumModelCutoffSlope_max;
	
	FitPar normPar= FitPar("Norm",1);
	
	FitPar gammaPar1= FitPar("Gamma1",gamma1,gamma1_min,gamma1_max);
	FitPar gammaPar2= FitPar("Gamma2",gamma2,gamma2_min,gamma2_max);
	FitPar gammaPar3= FitPar("Gamma3",gamma3,gamma3_min,gamma3_max);
	FitPar breakPar= FitPar("Break",lgS_break,lgS_break_min,lgS_break_max);
	FitPar cutoffPar= FitPar("Cutoff",lgS_cutoff,lgS_cutoff_min,lgS_cutoff_max);
	
	FitPar cutoffSlopePar= FitPar("Wc",Wc,Wc_min,Wc_max);
	
	
	/*	
	FitPar gammaPar1= FitPar("Gamma1",gamma1);
	FitPar gammaPar2= FitPar("Gamma2",gamma2);
	FitPar gammaPar3= FitPar("Gamma3",gamma3);
	FitPar breakPar= FitPar("Break",lgS_break);
	FitPar cutoffPar= FitPar("Cutoff",lgS_cutoff);
	*/


	// - Define spectrum pars
	SpectrumPars* spectrumPars= 0;
	if(gSpectrumModel_fit==ePowerLaw){
		FitPar gammaPar= FitPar("gammaPar",gamma2);
		spectrumPars= new PowerLawPars(normPar,gammaPar);
	}
	else if(gSpectrumModel_fit==eFlat){
		spectrumPars= new FlatSpectrumPars(normPar);
	}
	else if(gSpectrumModel_fit==eBrokenPowerLaws){
		spectrumPars= new BrokenPowerLawsPars(normPar,gammaPar1, gammaPar2, gammaPar3,breakPar,cutoffPar);
	}
	else if(gSpectrumModel_fit==eTwoBrokenPowerLaws){
		FitPar slope1Par= FitPar("gamma1",gamma1);
		FitPar slope2Par= FitPar("gamma2",gamma2);
		FitPar lgBreakPar= FitPar("lgS_break",lgS_break);
		//spectrumPars= new TwoBrokenPowerLawsPars(normPar,gammaPar2,gammaPar3,breakPar);
		spectrumPars= new TwoBrokenPowerLawsPars(normPar,slope1Par,slope2Par,lgBreakPar);
	}
	else if(gSpectrumModel_fit==ePowerLawWithCutoff){
		spectrumPars= new PowerLawWithCutoff(normPar,gammaPar2,cutoffPar,cutoffSlopePar);
	}
	else if(gSpectrumModel_fit==ePol3){	
		FitPar p0= FitPar("p0",gSpectrumModelPolPars[0]);
		FitPar p1= FitPar("p1",gSpectrumModelPolPars[1]);
		FitPar p2= FitPar("p2",gSpectrumModelPolPars[2]);
		FitPar p3= FitPar("p3",gSpectrumModelPolPars[3]);

		spectrumPars= new Pol3SpectrumPars(p0,p1,p2,p3);
	}
	else if(gSpectrumModel_fit==ePol6){	
		FitPar p0= FitPar("p0",gSpectrumModelPolPars[0]);
		FitPar p1= FitPar("p1",gSpectrumModelPolPars[1]);
		FitPar p2= FitPar("p2",gSpectrumModelPolPars[2]);
		FitPar p3= FitPar("p3",gSpectrumModelPolPars[3]);
		FitPar p4= FitPar("p4",gSpectrumModelPolPars[4]);
		FitPar p5= FitPar("p5",gSpectrumModelPolPars[5]);
		FitPar p6= FitPar("p6",gSpectrumModelPolPars[6]);

		spectrumPars= new Pol6SpectrumPars(p0,p1,p2,p3,p4,p5,p6);
	}

	//Fix spectrum pars?
	//spectrumPars->FixPar("Break",lgS_break);
	//spectrumPars->FixPar("Gamma1",gamma1);
	//spectrumPars->FixPar("Gamma2",gamma2);
	//spectrumPars->FixPar("Gamma3",gamma3);
	//spectrumPars->FixPar("Cutoff",lgS_cutoff);
	
	//====================================================
	//==         SET OPTIONS
	//=====================================================
	//- Use errors in chi2
	gForwardFolder->UseErrorsInChi2(gUseErrorsInChi2);	

	// - Use log likelihood 
	gForwardFolder->UseLikelihoodFit(gUseLikelihoodFit);

	//====================================================
	//==         RUN UNFOLDING
	//=====================================================
	int status= gForwardFolder->RunUnfold(
		gSourceCounts, gResponseMatrix, *spectrumPars,
		gComputeUncertainties,gNRandSamples,
		gUseFitRange,gLgSMin_fit,gLgSMax_fit
	);
	
	if(status<0) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Unfolding run failed!");
		#endif
		return -1;
	}

	//- Retrieve unfolding results
	gUnfoldedSpectrum= gForwardFolder->GetUnfoldedSpectrum();
	gUnfoldedSpectrum_totErrors= gForwardFolder->GetTotErrUnfoldedSpectrum();
	gFFSpectrum= gForwardFolder->GetForwardFoldedSpectrum();

	//Add source counts reso bias?
	if(gUseSourceCountsResoBiasCorr){
		for(int i=0;i<gUnfoldedSpectrum_totErrors->GetNbinsX();i++){
			for(int j=0;j<gUnfoldedSpectrum_totErrors->GetNbinsY();j++){
				double binContent= gUnfoldedSpectrum_totErrors->GetBinContent(i+1,j+1);
				double binError_stats= gUnfoldedSpectrum->GetBinError(i+1,j+1);
				double binError_unfolding= gUnfoldedSpectrum_totErrors->GetBinError(i+1,j+1);
				double binError_resoBias= gSourceCountsResoBiasCorrErr*binContent;
				//double binError_tot= sqrt( binError_unfolding*binError_unfolding + binError_resoBias*binError_resoBias);
				double binError_tot= sqrt( binError_stats*binError_stats + binError_unfolding*binError_unfolding + binError_resoBias*binError_resoBias);
				gUnfoldedSpectrum_totErrors->SetBinError(i+1,j+1,binError_tot);
			}
		}
	}//close if
	
	//Retrieve fit model fcn
	std::vector<double> fitPars= gForwardFolder->GetFitPars();
	
	if(gSpectrumModel_fit==eBrokenPowerLaws){
		gTrueFitModelFcn= new TF1("trueFitModelFcn",SpectrumUtils::BrokenPowerLawSpectrum,gLgSMin_fit,gLgSMax_fit,6);
		gTrueFitModelFcn->SetParameters(fitPars.data());
	}
	else if(gSpectrumModel_fit==ePowerLawWithCutoff){
		gTrueFitModelFcn= new TF1("trueFitModelFcn",SpectrumUtils::PowerLawWithCutoffSpectrum,gLgSMin_fit,gLgSMax_fit,4);
		gTrueFitModelFcn->SetParameters(fitPars.data());
	}
	else if(gSpectrumModel_fit==eTwoBrokenPowerLaws){
		gTrueFitModelFcn= new TF1("trueFitModelFcn",SpectrumUtils::TwoBrokenPowerLawSpectrum,gLgSMin_fit,gLgSMax_fit,4);
		gTrueFitModelFcn->SetParameters(fitPars.data());
	}
	else if(gSpectrumModel_fit==ePowerLaw){
		gTrueFitModelFcn= new TF1("trueFitModelFcn",SpectrumUtils::PowerLawSpectrum,gLgSMin_fit,gLgSMax_fit,2);
		gTrueFitModelFcn->SetParameters(fitPars.data());
	}
	else if(gSpectrumModel_fit==ePol3){
		gTrueFitModelFcn= new TF1("trueFitModelFcn",SpectrumUtils::Pol3Spectrum,gLgSMin_fit,gLgSMax_fit,4);
		gTrueFitModelFcn->SetParameters(fitPars.data());
	}
	else if(gSpectrumModel_fit==ePol6){
		gTrueFitModelFcn= new TF1("trueFitModelFcn",SpectrumUtils::Pol6Spectrum,gLgSMin_fit,gLgSMax_fit,7);
		gTrueFitModelFcn->SetParameters(fitPars.data());
	}

		
	

	//====================================================
	//==         COMPUTE DIFF COUNTS
	//=====================================================
	if(!gDiffSourceCounts) gDiffSourceCounts= (TH1D*)gSourceCounts->Clone("DiffSourceCounts");
	gDiffSourceCounts->Reset();

	if(!gUnfoldedDiffSpectrum) gUnfoldedDiffSpectrum= (TH1D*)gUnfoldedSpectrum->Clone("UnfoldedDiffSpectrum");
	gUnfoldedDiffSpectrum->Reset();

	if(!gUnfoldedDiffSpectrum_totErrors) gUnfoldedDiffSpectrum_totErrors= (TH1D*)gUnfoldedSpectrum->Clone("UnfoldedDiffSpectrum_totErrors");
	gUnfoldedDiffSpectrum_totErrors->Reset();

	if(!gFFDiffSpectrum) gFFDiffSpectrum= (TH1D*)gFFSpectrum->Clone("FFDiffSpectrum");
	gFFDiffSpectrum->Reset();

	if(!gUnfoldedDiffSpectrum_resoBiasCorrected) gUnfoldedDiffSpectrum_resoBiasCorrected= (TH1D*)gUnfoldedDiffSpectrum->Clone("UnfoldedDiffSpectrum_resoBiasCorrected");
	gUnfoldedDiffSpectrum_resoBiasCorrected->Reset();

	if(!gSourceCountsFitResidualHisto) gSourceCountsFitResidualHisto= (TH1D*)gSourceCounts->Clone("SourceCountsResiduals");
	gSourceCountsFitResidualHisto->Reset();

	double gamma= gNormSpectralIndex;
	double A_sr= Deg2ToSr(gArea);

	for(int i=0;i<gSourceCounts->GetNbinsX();i++)
	{
		double N= gSourceCounts->GetBinContent(i+1);
		double NErr= gSourceCounts->GetBinError(i+1); 
		double N_unfolded= gUnfoldedSpectrum->GetBinContent(i+1);
		double NErr_unfolded= gUnfoldedSpectrum->GetBinError(i+1);
		double NTotErr_unfolded= gUnfoldedSpectrum_totErrors->GetBinError(i+1);
		double N_ff= gFFSpectrum->GetBinContent(i+1);
		double NErr_ff= gFFSpectrum->GetBinError(i+1);
		
		if(N>0){
			double fitRelResidual= N_ff/N-1;
			gSourceCountsFitResidualHisto->SetBinContent(i+1,fitRelResidual);
		}

		//if(N<=0) continue;
		double lgS= gSourceCounts->GetBinCenter(i+1);		
		double dlgS= gSourceCounts->GetBinWidth(i+1);
		double lgS_min= gSourceCounts->GetBinLowEdge(i+1);
		double lgS_max= lgS_min + dlgS;
		double S_min= pow(10,lgS_min);
		double S_max= pow(10,lgS_max);
		double dS= S_max-S_min;
		double S= pow(10,lgS);		
		double NormFactor= pow(S,-gamma);
			
		//Set differential source counts
		double counts= N/(dS*A_sr);//in units: Jy^-1 sr^-1 
		double counts_norm= counts/NormFactor;//in units: Jy^-2.5 sr^-1
		double countsErr= NErr/(dS*A_sr);
		double countsErr_norm= countsErr/NormFactor;

		gDiffSourceCounts->SetBinContent(i+1,counts_norm);
		gDiffSourceCounts->SetBinError(i+1,countsErr_norm);
			
		//Set differential unfolded source counts
		double counts_unfolded= N_unfolded/(dS*A_sr);//in units: Jy^-1 sr^-1 
		double counts_unfolded_norm= counts_unfolded/NormFactor;//in units: Jy^-2.5 sr^-1
		double countsErr_unfolded= NErr_unfolded/(dS*A_sr);
		double countsErr_unfolded_norm= countsErr_unfolded/NormFactor;
		double countsTotErr_unfolded= NTotErr_unfolded/(dS*A_sr);
		double countsTotErr_unfolded_norm= countsTotErr_unfolded/NormFactor;

		gUnfoldedDiffSpectrum->SetBinContent(i+1,counts_unfolded_norm);
		gUnfoldedDiffSpectrum->SetBinError(i+1,countsErr_unfolded_norm);

		gUnfoldedDiffSpectrum_totErrors->SetBinContent(i+1,counts_unfolded_norm);
		gUnfoldedDiffSpectrum_totErrors->SetBinError(i+1,countsTotErr_unfolded_norm);

		//Set differential forward folded source counts
		double counts_ff= N_ff/(dS*A_sr);//in units: Jy^-1 sr^-1 
		double counts_ff_norm= counts_ff/NormFactor;//in units: Jy^-2.5 sr^-1
		double countsErr_ff= NErr_ff/(dS*A_sr);
		double countsErr_ff_norm= countsErr_ff/NormFactor;

		gFFDiffSpectrum->SetBinContent(i+1,counts_ff_norm);
		gFFDiffSpectrum->SetBinError(i+1,countsErr_ff_norm);
		

		cout<<"INFO: bin "<<i+1<<": N="<<N<<" +- "<<NErr<<", N(unfolded)="<<N_unfolded<<" +- "<<NErr_unfolded<<"(stat) +- "<<NTotErr_unfolded<<"(stat+syst), S(mJy)="<<S*1000<<", S(mJy)=["<<S_min*1000<<","<<S_max*1000<<"], dS(mJy)="<<dS*1000.<<", lgS=["<<lgS_min<<","<<lgS_max<<"], normFactor="<<NormFactor<<", counts="<<counts<<" +- "<<countsErr<<", counts(norm)="<<counts_norm<<" +- "<<countsErr_norm<<", counts(unfolded)="<<counts_unfolded<<" +- "<<countsErr_unfolded<<"(stat) +- "<<countsTotErr_unfolded<<"(stat+syst), counts(unfolded_norm)="<<counts_unfolded_norm<<" +- "<<countsErr_unfolded_norm<<"(stat) +- "<<countsTotErr_unfolded_norm<<" (stat+syst)"<<endl;
	
		//Compute reso bias corrected counts
		double resoBiasCorr= gSourceCountsResoBiasCorrFcn->Eval(lgS);
		double counts_unfolded_norm_corr= resoBiasCorr* counts_unfolded_norm;
		double countsErr_unfolded_norm_corr= resoBiasCorr* countsErr_unfolded_norm;
		gUnfoldedDiffSpectrum_resoBiasCorrected->SetBinContent(i+1,counts_unfolded_norm_corr);
		gUnfoldedDiffSpectrum_resoBiasCorrected->SetBinError(i+1,countsErr_unfolded_norm_corr);
		

	}//end loop bins


	

	return 0;

}//close UnfoldData()



int SetResponseMatrix()
{	
	//Set flux bins
	double lgSmin_true= gTrueBins[0];
	double lgSmax_true= gTrueBins[gTrueBins.size()-1];
	double lgSmin_rec= gRecBins[0];
	double lgSmax_rec= gRecBins[gRecBins.size()-1];

	//====================================================
	//==         SET RESPONSE MODEL PARS
	//=====================================================		
	//- Set power law spectrum pars
	double gamma1= gSpectrumModelIndex1;
	double gamma2= gSpectrumModelIndex2;
	double gamma3= gSpectrumModelIndex3;
	double gamma_min= 0;
	double gamma_max= 5;
	double lgS_break= gSpectrumModelLgSBreak;
	double lgS_cutoff= gSpectrumModelLgSCutoff;
	double lgS_min= -5;
	double lgS_max= 10;
	FitPar normPar= FitPar("Norm",1);
	FitPar gammaPar1= FitPar("Gamma1",gamma1,gamma_min,gamma_max);
	FitPar gammaPar2= FitPar("Gamma2",gamma2,gamma_min,gamma_max);
	FitPar gammaPar3= FitPar("Gamma3",gamma3,gamma_min,gamma_max);
	FitPar breakPar= FitPar("Break",lgS_break,lgS_min,lgS_max);
	FitPar cutoffPar= FitPar("Cutoff",lgS_cutoff,lgS_min,lgS_max);

	// - Set 3rd deg pol pars
	

	// - Set spectrum pars
	SpectrumPars* spectrumPars= 0;
	if(gSpectrumModel==ePowerLaw){
		gSpectrumPars= new PowerLawPars(normPar,gammaPar2);
	}
	else if(gSpectrumModel==eFlat){
		gSpectrumPars= new FlatSpectrumPars(normPar);
	}
	else if(gSpectrumModel==eBrokenPowerLaws){
		gSpectrumPars= new BrokenPowerLawsPars(normPar,gammaPar1, gammaPar2, gammaPar3,breakPar,cutoffPar);
	}
	else if(gSpectrumModel==eTwoBrokenPowerLaws){
		gSpectrumPars= new TwoBrokenPowerLawsPars(normPar,gammaPar2, gammaPar3,breakPar);
	}
	else if(gSpectrumModel==ePol3){
		FitPar p0= FitPar("p0",gSpectrumModelPolPars[0]);
		FitPar p1= FitPar("p1",gSpectrumModelPolPars[1]);
		FitPar p2= FitPar("p2",gSpectrumModelPolPars[2]);
		FitPar p3= FitPar("p3",gSpectrumModelPolPars[3]);

		gSpectrumPars= new Pol3SpectrumPars(p0,p1,p2,p3);
	}
	else if(gSpectrumModel==ePol6){	
		FitPar p0= FitPar("p0",gSpectrumModelPolPars[0]);
		FitPar p1= FitPar("p1",gSpectrumModelPolPars[1]);
		FitPar p2= FitPar("p2",gSpectrumModelPolPars[2]);
		FitPar p3= FitPar("p3",gSpectrumModelPolPars[3]);
		FitPar p4= FitPar("p4",gSpectrumModelPolPars[4]);
		FitPar p5= FitPar("p5",gSpectrumModelPolPars[5]);
		FitPar p6= FitPar("p6",gSpectrumModelPolPars[6]);

		gSpectrumPars= new Pol6SpectrumPars(p0,p1,p2,p3,p4,p5,p6);
	}
	
	// - Set bias pars
	BiasPars* biasPars= 0;
	//if(gBiasModel==eCONST_BIAS) biasPars= new ConstBiasPars(gConstBias);
	if(gBiasModel==eCONST_BIAS) biasPars= new ConstBiasPars(gFluxBiasPars[0]);
	else if(gBiasModel==eEXP_BIAS) biasPars= new ExpBiasPars(gFluxBiasPars[0],gFluxBiasPars[1],gFluxBiasPars[2]);
	else if(gBiasModel==eEXP_STEP_BIAS) biasPars= new ExpStepBiasPars(gFluxBiasPars[0],gFluxBiasPars[1],gFluxBiasPars[2],gFluxBiasPars[3]);

	gBiasModelFcn= SpectrumUtils::ComputeBiasModel(*biasPars,lgSmin_true,lgSmax_true,200);

	// - Set resolution pars
	ResoPars* resoPars= 0;
	//if(gResoModel==eCONST_RESO) resoPars= new ConstResoPars(gConstReso);
	if(gResoModel==eCONST_RESO) resoPars= new ConstResoPars(gFluxResoPars[0]);
	else if(gResoModel==eEXP_RESO) resoPars= new ExpResoPars(gFluxResoPars[0],gFluxResoPars[1]);
	else if(gResoModel==eEXP_STEP_RESO) resoPars= new ExpStepResoPars(gFluxResoPars[0],gFluxResoPars[1],gFluxResoPars[2]);

	gResoModelFcn= SpectrumUtils::ComputeResoModel(*resoPars,lgSmin_true,lgSmax_true,200);

	// - Set efficiency pars
	EfficiencyPars* effPars= 0;
	//if(gEfficiencyModel==eCONST_EFF) effPars= new ConstEfficiencyPars(gConstEfficiency);
	if(gEfficiencyModel==eCONST_EFF) effPars= new ConstEfficiencyPars(gEfficiencyPars[0]);
	else if(gEfficiencyModel==eSIGMOID_EFF) effPars= new SigmoidEfficiencyPars(gEfficiencyPars[0],gEfficiencyPars[1],gEfficiencyPars[2]);

	
	
	//====================================================
	//==         BUILD RESPONSE MATRIX
	//=====================================================
	

	/*
	gResponseModelFcn= SpectrumUtils::ComputeResponseModel(
		*gSpectrumPars,
		*biasPars,
		*resoPars,
		*effPars,
		lgSmin_true,lgSmax_true,
		lgSmin_rec,lgSmax_rec
	);
	*/
	
	//gResponseMatrix= SpectrumUtils::ComputeParametricResponse(
	gResponseMatrix= SpectrumUtils::BuildResponseMatrix(
		*gSpectrumPars,
		*biasPars,
		*resoPars,
		*effPars,
		gTrueBins, gRecBins
	);

	if(!gResponseMatrix){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the response matrix!");
		#endif
		return -1;
	}


	//=========================================================
	//==  COMPUTE SOURCE COUNT RESOLUTION BIAS CORR FACTOR
	//=========================================================	
	std::vector<double> scResoBiasPars
	{
		gBmaj,
		gBmin,
		gNoiseRms,
		gSThr,
		gResolvedSourceThrP0,
		gResolvedSourceThrP1,
		gResolvedSourceThrSlope,
		gDataFrequency,
		gDataSpectralIndex,
		gPhiMedianSlope,	
		gPhiMedianSbreak,
		gPhiMedianScaleFactor
	};

	gSourceCountsResoBiasCorrFcn= new TF1(
		"SourceCountsResoBiasCorrFcn",
		SpectrumUtils::SourceCountsResoBiasCorrFactor,
		lgSmin_rec,lgSmax_rec,
		12
	);
	gSourceCountsResoBiasCorrFcn->SetParameters(scResoBiasPars.data());
		

	if(gUseSourceCountsResoBiasCorr){
		for(int i=0;i<gResponseMatrix->GetNbinsX();i++){
			double lgS= gResponseMatrix->GetXaxis()->GetBinCenter(i+1);
			double resoBiasCorr= gSourceCountsResoBiasCorrFcn->Eval(lgS);
			
			for(int j=0;j<gResponseMatrix->GetNbinsY();j++){
				double w_old= gResponseMatrix->GetBinContent(i+1,j+1);
				double w_new= w_old/resoBiasCorr;
				gResponseMatrix->SetBinContent(i+1,j+1,w_new);
			}//end loop bins y
		}//end loop bins x
	}//close useSourceCountsResoBiasCorr

	//Compute response matrix errors
	if(ComputeResponseMatrixErrors()<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute response matrix uncertainties ...");
		#endif
	}

	return 0;

}//close SetResponseMatrix()


int ComputeResponseMatrixErrors()
{
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing errors on response matrix elements...");
	#endif

	//## Generate random bias model pars
	// - Flux bias model
	#ifdef LOGGING_ENABLED
		INFO_LOG("Generating random flux bias model pars ...");
	#endif
	std::vector<std::vector<double>> fluxBiasRandPars;
	if(StatsUtils::GenerateFitParsAroundCovMatrix(fluxBiasRandPars,gFluxBiasPars,gFluxBiasParCovMatrix,gNRandSamples)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to generate random bias model pars!");
		#endif
		return -1;
	}
	
	// - Flux reso model
	#ifdef LOGGING_ENABLED
		INFO_LOG("Generating random flux reso model pars ...");
	#endif
	std::vector<std::vector<double>> fluxResoRandPars;
	if(StatsUtils::GenerateFitParsAroundCovMatrix(fluxResoRandPars,gFluxResoPars,gFluxResoParCovMatrix,gNRandSamples)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to generate random reso model pars!");
		#endif
		return -1;
	}

	// - Efficiency model
	#ifdef LOGGING_ENABLED
		INFO_LOG("Generating random efficiency model pars ...");
	#endif
	std::vector<std::vector<double>> efficiencyRandPars;
	if(StatsUtils::GenerateFitParsAroundCovMatrix(efficiencyRandPars,gEfficiencyPars,gEfficiencyParCovMatrix,gNRandSamples)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to generate random efficiency model pars!");
		#endif
		return -1;
	}

	//# Compute response matrix uncertainties
	std::vector<TH2D*> gRandResponseMatrices;
	gRandMeanResponseMatrix= (TH2D*)gResponseMatrix->Clone("randMeanResponseMatrix");
	gRandMeanResponseMatrix->Reset();

	for(int k=0;k<gNRandSamples;k++)
	{
		//Generate random bias model
		BiasPars* biasPars_rand= 0;
		if(gBiasModel==eCONST_BIAS) biasPars_rand= new ConstBiasPars(fluxBiasRandPars[k][0]);
		else if(gBiasModel==eEXP_BIAS) biasPars_rand= new ExpBiasPars(fluxBiasRandPars[k][0],fluxBiasRandPars[k][1],fluxBiasRandPars[k][2]);
		else if(gBiasModel==eEXP_STEP_BIAS) biasPars_rand= new ExpStepBiasPars(fluxBiasRandPars[k][0],fluxBiasRandPars[k][1],fluxBiasRandPars[k][2],fluxBiasRandPars[k][3]);
		
		//Generate random resolution pars
		ResoPars* resoPars_rand= 0;
		if(gResoModel==eCONST_RESO) resoPars_rand= new ConstResoPars(fluxResoRandPars[k][0]);
		else if(gResoModel==eEXP_RESO) resoPars_rand= new ExpResoPars(fluxResoRandPars[k][0],fluxResoRandPars[k][1]);
		else if(gResoModel==eEXP_STEP_RESO) resoPars_rand= new ExpStepResoPars(fluxResoRandPars[k][0],fluxResoRandPars[k][1],fluxResoRandPars[k][2]);

		//Generate random efficiency pars
		EfficiencyPars* effPars_rand= 0;
		if(gEfficiencyModel==eCONST_EFF) effPars_rand= new ConstEfficiencyPars(efficiencyRandPars[k][0]);
		else if(gEfficiencyModel==eSIGMOID_EFF) effPars_rand= new SigmoidEfficiencyPars(efficiencyRandPars[k][0],efficiencyRandPars[k][1],efficiencyRandPars[k][2]);

		//Create response matrix
		TH2D* responseMatrix_rand= SpectrumUtils::BuildResponseMatrix(
			*gSpectrumPars,
			*biasPars_rand,
			*resoPars_rand,
			*effPars_rand,
			gTrueBins, gRecBins
		);
		
		//Add to matrix
		gRandResponseMatrices.push_back(responseMatrix_rand);
		gRandMeanResponseMatrix->Add(responseMatrix_rand);

		//Clear data
		if(biasPars_rand){
			delete biasPars_rand;
			biasPars_rand= 0;
		}
		if(resoPars_rand){
			delete resoPars_rand;
			resoPars_rand= 0;
		}
		if(effPars_rand){
			delete effPars_rand;
			effPars_rand= 0;
		}

	}//end loop rand samples

	//Set random response matrix in forward folder
	gForwardFolder->SetRandomResponseMatrices(gRandResponseMatrices);

	//Compute mean random response matrix
	gRandMeanResponseMatrix->Scale(1./gNRandSamples);

	//Compute sigma random response matrix
	gRandSigmaResponseMatrix= (TH2D*)gResponseMatrix->Clone("randSigmaResponseMatrix");
	gRandSigmaResponseMatrix->Reset();
	
	for(size_t k=0;k<gRandResponseMatrices.size();k++)
	{
		for(int i=0;i<gRandResponseMatrices[k]->GetNbinsX();i++){
			for(int j=0;j<gRandResponseMatrices[k]->GetNbinsY();j++){
				int binId= gRandSigmaResponseMatrix->GetBin(i+1,j+1);
				double w= gRandResponseMatrices[k]->GetBinContent(i+1,j+1);
				double w_mean= gRandSigmaResponseMatrix->GetBinContent(i+1,j+1);
				double s2= (w-w_mean)*(w-w_mean);
				gRandSigmaResponseMatrix->AddBinContent(binId,s2);
			}//end loop bins y
		}//end loop bins
	}//end loop rand samples

	gRandSigmaResponseMatrix->Scale(1./(gNRandSamples-1));

	for(int i=0;i<gRandSigmaResponseMatrix->GetNbinsX();i++){
		for(int j=0;j<gRandSigmaResponseMatrix->GetNbinsY();j++){
			double w= gRandSigmaResponseMatrix->GetBinContent(i+1,j+1);
			double wsqr= sqrt(w);
			gRandSigmaResponseMatrix->SetBinContent(i+1,j+1,wsqr);
		}
	}
	
	
	return 0;

}//close ComputeResponseMatrixErrors()


int ComputeSourceCountsModels()
{
	//Set reference frequency
	double refFreq= 1.4;//GHz
	double freqScaleFactor= pow(refFreq/gDataFrequency,gDataSpectralIndex);
	double lgFreqScaleFactor= log10(freqScaleFactor);

	double alphaEG_min= -0.9;
	double alphaEG_max= -0.7;
	double freqScaleFactor_min= pow(refFreq/gDataFrequency,alphaEG_min);
	double freqScaleFactor_max= pow(refFreq/gDataFrequency,alphaEG_max);
	double lgFreqScaleFactor_min= log10(freqScaleFactor_min);
	double lgFreqScaleFactor_max= log10(freqScaleFactor_max);

	double galSourceOverDensity= 0.2;
	double galSourceOverDensityFactor= 1. + galSourceOverDensity;

	//cout<<"refFreq(GHz)="<<refFreq<<", freq(GHz)="<<gDataFrequency<<", alpha="<<gDataSpectralIndex<<", scaleFactor="<<freqScaleFactor<<", lgFreqScaleFactor="<<lgFreqScaleFactor<<endl;

	//=======================================
	//===    KATGERT+90 SOURCE COUNTS
	//=======================================
	gExpDataFitFcn_Katgert= new TF1(
		"expDataFitFcn_Katgert",
		//"pow(10,([0] + [1]*(x+3) + [2]*pow(x+3,2) + [3]*pow(x+3,3)) )",
		//"pow(10,([0] + [1]*([4]+x+3) + [2]*pow([4]+x+3,2) + [3]*pow([4]+x+3,3)) )",
		"pow(10,([0] + [1]*([4]+x+3) + [2]*pow([4]+x+3,2) + [3]*pow([4]+x+3,3)) -2.5*[4] )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	//gExpDataFitFcn_Katgert->SetParameters(gSourceCountsExpDataFitPars_Katgert.data());
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Katgert.size();k++) gExpDataFitFcn_Katgert->SetParameter(k,gSourceCountsExpDataFitPars_Katgert[k]);
	gExpDataFitFcn_Katgert->SetParameter(gSourceCountsExpDataFitPars_Katgert.size(),lgFreqScaleFactor);
	gExpDataFitFcn_Katgert->SetLineColor(kBlack);
	gExpDataFitFcn_Katgert->SetLineStyle(kDashed);
	//gExpDataFitFcn_Katgert->Draw("l same");

	TF1* expDataFitFcn_fluxScaledMin_Katgert= new TF1(
		"expDataFitFcn_fluxScaledMin_Katgert",
		"pow(10,([0] + [1]*([4]+x+3) + [2]*pow([4]+x+3,2) + [3]*pow([4]+x+3,3)) -2.5*[4] )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Katgert.size();k++) expDataFitFcn_fluxScaledMin_Katgert->SetParameter(k,gSourceCountsExpDataFitPars_Katgert[k]);
	expDataFitFcn_fluxScaledMin_Katgert->SetParameter(gSourceCountsExpDataFitPars_Katgert.size(),lgFreqScaleFactor_min);

	TF1* expDataFitFcn_fluxScaledMax_Katgert= new TF1(
		"expDataFitFcn_fluxScaledMax_Katgert",
		"pow(10,([0] + [1]*([4]+x+3) + [2]*pow([4]+x+3,2) + [3]*pow([4]+x+3,3)) -2.5*[4] )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Katgert.size();k++) expDataFitFcn_fluxScaledMax_Katgert->SetParameter(k,gSourceCountsExpDataFitPars_Katgert[k]);
	expDataFitFcn_fluxScaledMax_Katgert->SetParameter(gSourceCountsExpDataFitPars_Katgert.size(),lgFreqScaleFactor_max);

	//Shift EG source counts model by +20% (to account galactic source overdensity)
	
	TF1* expDataFitFcn_fluxScaled_overdensity_Katgert= new TF1(
		"expDataFitFcn_fluxScaled_overdensity_Katgert",
		"pow(10,([0] + [1]*([4]+x+3) + [2]*pow([4]+x+3,2) + [3]*pow([4]+x+3,3)) -2.5*[4] + log10([5]) )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Katgert.size();k++) expDataFitFcn_fluxScaled_overdensity_Katgert->SetParameter(k,gSourceCountsExpDataFitPars_Katgert[k]);
	expDataFitFcn_fluxScaled_overdensity_Katgert->SetParameter(gSourceCountsExpDataFitPars_Katgert.size(),lgFreqScaleFactor);
	expDataFitFcn_fluxScaled_overdensity_Katgert->SetParameter(gSourceCountsExpDataFitPars_Katgert.size()+1,galSourceOverDensityFactor);
	expDataFitFcn_fluxScaled_overdensity_Katgert->SetLineColor(kRed);
	expDataFitFcn_fluxScaled_overdensity_Katgert->SetLineStyle(kDashed);
	//expDataFitFcn_fluxScaled_overdensity_Katgert->Draw("l same");
	

	//=======================================
	//===      HOPKINS+03 SOURCE COUNTS
	//=======================================
	gExpDataFitFcn_Hopkins= new TF1(
		"expDataFitFcn_Hopkins",
		//"pow(10,([0] + [1]*(x+3) + [2]*pow(x+3,2) + [3]*pow(x+3,3) + [4]*pow(x+3,4) + [5]*pow(x+3,5) + [6]*pow(x+3,6) ) )",
		"pow(10,([0] + [1]*([7]+x+3) + [2]*pow([7]+x+3,2) + [3]*pow([7]+x+3,3) + [4]*pow([7]+x+3,4) + [5]*pow([7]+x+3,5) + [6]*pow([7]+x+3,6) ) -2.5*[7] )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	//gExpDataFitFcn_Hopkins->SetParameters(gSourceCountsExpDataFitPars_Hopkins.data());
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Hopkins.size();k++) gExpDataFitFcn_Hopkins->SetParameter(k,gSourceCountsExpDataFitPars_Hopkins[k]);
	gExpDataFitFcn_Hopkins->SetParameter(gSourceCountsExpDataFitPars_Hopkins.size(),lgFreqScaleFactor);
	gExpDataFitFcn_Hopkins->SetLineColor(kBlack);
	gExpDataFitFcn_Hopkins->SetLineStyle(kDotted);
	//gExpDataFitFcn_Hopkins->Draw("l same");

	
	//Compute source count models shifted to input data frequency assuming a min/max spectral index for EG sources
	TF1* expDataFitFcn_fluxScaledMin_Hopkins= new TF1(
		"expDataFitFcn_fluxScaledMin_Hopkins",
		"pow(10,([0] + [1]*([7]+x+3) + [2]*pow([7]+x+3,2) + [3]*pow([7]+x+3,3) + [4]*pow([7]+x+3,4) + [5]*pow([7]+x+3,5) + [6]*pow([7]+x+3,6) ) -2.5*[7] )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Hopkins.size();k++) expDataFitFcn_fluxScaledMin_Hopkins->SetParameter(k,gSourceCountsExpDataFitPars_Hopkins[k]);
	expDataFitFcn_fluxScaledMin_Hopkins->SetParameter(gSourceCountsExpDataFitPars_Hopkins.size(),lgFreqScaleFactor_min);

	TF1* expDataFitFcn_fluxScaledMax_Hopkins= new TF1(
		"expDataFitFcn_fluxScaledMax_Hopkins",
		"pow(10,([0] + [1]*([7]+x+3) + [2]*pow([7]+x+3,2) + [3]*pow([7]+x+3,3) + [4]*pow([7]+x+3,4) + [5]*pow([7]+x+3,5) + [6]*pow([7]+x+3,6) ) -2.5*[7] )",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Hopkins.size();k++) expDataFitFcn_fluxScaledMax_Hopkins->SetParameter(k,gSourceCountsExpDataFitPars_Hopkins[k]);
	expDataFitFcn_fluxScaledMax_Hopkins->SetParameter(gSourceCountsExpDataFitPars_Hopkins.size(),lgFreqScaleFactor_max);

	TF1* expDataFitFcn_fluxScaled_overdensity_Hopkins= new TF1(
		"expDataFitFcn_fluxScaled_overdensity_Hopkins",
		"pow(10,([0] + [1]*([7]+x+3) + [2]*pow([7]+x+3,2) + [3]*pow([7]+x+3,3) + [4]*pow([7]+x+3,4) + [5]*pow([7]+x+3,5) + [6]*pow([7]+x+3,6) ) -2.5*[7] + log10([8]))",
		gTrueBins[0],gTrueBins[gTrueBins.size()-1]
	);
	for(size_t k=0;k<gSourceCountsExpDataFitPars_Hopkins.size();k++) expDataFitFcn_fluxScaled_overdensity_Hopkins->SetParameter(k,gSourceCountsExpDataFitPars_Hopkins[k]);
	expDataFitFcn_fluxScaled_overdensity_Hopkins->SetParameter(gSourceCountsExpDataFitPars_Hopkins.size(),lgFreqScaleFactor);
	expDataFitFcn_fluxScaled_overdensity_Hopkins->SetParameter(gSourceCountsExpDataFitPars_Hopkins.size()+1,galSourceOverDensityFactor);
	expDataFitFcn_fluxScaled_overdensity_Hopkins->SetLineColor(kRed);
	expDataFitFcn_fluxScaled_overdensity_Hopkins->SetLineStyle(kDotted);
	//expDataFitFcn_fluxScaled_overdensity_Hopkins->Draw("l same");

	//=======================================
	//===      SOURCE COUNTS MODEL
	//=======================================
	gExpDataFitAreaGraph_min= new TGraph;	
	gExpDataFitAreaGraph_max= new TGraph;
	gExpDataFitAreaGraph_mean= new TGraph;	
	double xmin = gDiffSourceCounts->GetXaxis()->GetXmin();
	double xmax = gDiffSourceCounts->GetXaxis()->GetXmax();
	int npx= expDataFitFcn_fluxScaledMin_Katgert->GetNpx();
	double dx = (xmax-xmin)/npx;
	double x= xmin + 0.5*dx; 
	int npoints= 0;
	std::vector<double> xlist;
	std::vector<double> yminlist;
	std::vector<double> ymaxlist;
	while (x <= xmax)
	{
		xlist.push_back(x);
		double y_Katgert_min= expDataFitFcn_fluxScaledMin_Katgert->Eval(x);
		double y_Katgert_max= expDataFitFcn_fluxScaledMax_Katgert->Eval(x);
		double y_Hopkins_min= expDataFitFcn_fluxScaledMin_Hopkins->Eval(x);
		double y_Hopkins_max= expDataFitFcn_fluxScaledMax_Hopkins->Eval(x);		
		double ymodel_min= std::min({y_Katgert_min,y_Katgert_max,y_Hopkins_min,y_Hopkins_max});  	
		double ymodel_max= std::max({y_Katgert_min,y_Katgert_max,y_Hopkins_min,y_Hopkins_max});  	
		double ymodel_mean= ymodel_min + 0.5*(ymodel_max-ymodel_min);
		
		//Enlarge the boundary by the stat/syst error in source count model (~20%)
		ymodel_min*= (1-gExpSourceCountDataStatErr);
		ymodel_max*= (1+gExpSourceCountDataStatErr);
		
		yminlist.push_back(ymodel_min);
		ymaxlist.push_back(ymodel_max);
		gExpDataFitAreaGraph_min->SetPoint(npoints,x,ymodel_min);
		gExpDataFitAreaGraph_max->SetPoint(npoints,x,ymodel_max);
		gExpDataFitAreaGraph_mean->SetPoint(npoints,x,ymodel_mean);
		npoints++;
		x+= dx;
		
	}//end loop points
		
	gExpDataFitAreaGraph= new TGraph(npoints*2);
	for (int i=0;i<npoints;i++) {
  	gExpDataFitAreaGraph->SetPoint(i,xlist[i],ymaxlist[i]);
    gExpDataFitAreaGraph->SetPoint(npoints+i,xlist[npoints-i-1],yminlist[npoints-i-1]);
  }
	

	/*
	// - Compute true differential source counts
	TH1D* trueDiffModelCounts= new TH1D("trueDiffModelCounts","",10,gLgSMin_fit,gLgSMax_fit);
	//TH1D* trueDiffModelCounts= (TH1D*)gSourceCounts->Clone("trueDiffModelCounts");
	double gamma= gNormSpectralIndex;
	double A_sr= Deg2ToSr(gArea);

	for(int i=0;i<trueDiffModelCounts->GetNbinsX();i++){
		double lgS= trueDiffModelCounts->GetBinCenter(i+1);
		double dlgS= trueDiffModelCounts->GetBinWidth(i+1);
		double lgS_min= trueDiffModelCounts->GetBinLowEdge(i+1);
		double lgS_max= lgS_min + dlgS;
		double S_min= pow(10,lgS_min);
		double S_max= pow(10,lgS_max);
		double dS= S_max-S_min;
		double S= pow(10,lgS);
		double NormFactor= pow(S,-gamma);

		double N= gTrueFitModelFcn->Integral(lgS_min,lgS_max);
		//double N= gTrueFitModelFcn->Eval(lgS);
		//trueDiffModelCounts->SetBinContent(i+1,N);
		//if(N>0) trueDiffModelCounts->SetBinError(i+1,sqrt(N));
		
		double counts= N/(dS*A_sr);//in units: Jy^-1 sr^-1
		//double counts= N/(A_sr);//in units: Jy^-1 sr^-1 
		//double counts= N/A_sr;
		//double counts= N/(dS*A_sr*fabs(dlgS));//in units: Jy^-1 sr^-1
		//double counts= trueDiffModelCounts->GetBinContent(i+1)/(dS*A_sr);//in units: Jy^-1 sr^-1 
		 	
		double counts_norm= counts/NormFactor;//in units: Jy^-2.5 sr^-1

		cout<<"lgS="<<lgS<<", dS="<<dS<<", Nmodel="<<N<<", counts="<<counts<<endl;
		trueDiffModelCounts->SetBinContent(i+1,counts_norm);
	}
	
	trueDiffModelCounts->SetLineColor(kMagenta);
	//trueDiffModelCounts->Draw("L hist same");
	*/

	//## Clear data
	if(expDataFitFcn_fluxScaledMin_Katgert) expDataFitFcn_fluxScaledMin_Katgert->Delete();
	if(expDataFitFcn_fluxScaledMax_Katgert) expDataFitFcn_fluxScaledMax_Katgert->Delete();
	if(expDataFitFcn_fluxScaled_overdensity_Katgert) expDataFitFcn_fluxScaled_overdensity_Katgert->Delete();
	if(expDataFitFcn_fluxScaledMin_Hopkins) expDataFitFcn_fluxScaledMin_Hopkins->Delete();
	if(expDataFitFcn_fluxScaledMax_Hopkins) expDataFitFcn_fluxScaledMax_Hopkins->Delete();
	if(expDataFitFcn_fluxScaled_overdensity_Hopkins) expDataFitFcn_fluxScaled_overdensity_Hopkins->Delete();

	return 0;

}//close ComputeSourceCountsModels()


void Draw()
{
	//=======================================
	//==    RESPONSE MATRIX
	//=======================================
	// - Draw response function 2D
	if(gResponseModelFcn){
		TCanvas* ResponseModelPlot= new TCanvas("ResponseModelPlot","ResponseModelPlot");
		ResponseModelPlot->cd();

		gResponseModelFcn->Draw("COLZ");
	}

	// - Draw response matrix
	TF1* diagFcn= new TF1("diagFcn","x",gTrueBins[0],gTrueBins[gTrueBins.size()-1]);
	diagFcn->SetLineColor(kWhite);

	if(gResponseMatrix){
		gStyle->SetPalette(kRainBow);

		TCanvas* ResponseMatrixPlot= new TCanvas("ResponseMatrixPlot","ResponseMatrixPlot");
		ResponseMatrixPlot->cd();

		gResponseMatrix->SetStats(0);
		gResponseMatrix->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
		gResponseMatrix->GetYaxis()->SetTitle("log_{10}(S_{meas}/Jy)");
		gResponseMatrix->GetXaxis()->CenterTitle(true);
		gResponseMatrix->GetYaxis()->CenterTitle(true);
		gResponseMatrix->GetZaxis()->SetRangeUser(-0.0001,1);
		gResponseMatrix->Draw("COLZ");

		//Draw mJy axis
		double ymin_h= gResponseMatrix->GetYaxis()->GetXmin();
		double ymax_h= gResponseMatrix->GetYaxis()->GetXmax();
		double xmin_h= gResponseMatrix->GetXaxis()->GetXmin();
		double xmax_h= gResponseMatrix->GetXaxis()->GetXmax();
		double xmin_shifted= xmin_h + 3;
		double xmax_shifted= xmax_h + 3;
		double ymin_shifted= ymin_h + 3;
		double ymax_shifted= ymax_h + 3;	

		TGaxis* ResponseMatrixPlotXAxis= new TGaxis(
			xmin_h,ymin_h,
			xmax_h,ymin_h,
			xmin_shifted,xmax_shifted,510,"+L"
		);
		ResponseMatrixPlotXAxis->SetLineColor(kBlack);
		ResponseMatrixPlotXAxis->SetTitleFont(52);
		ResponseMatrixPlotXAxis->SetTitleSize(0.06);	
		ResponseMatrixPlotXAxis->SetTitleOffset(0.85);
		ResponseMatrixPlotXAxis->SetLabelFont(42);
		ResponseMatrixPlotXAxis->SetLabelSize(0.04);
		ResponseMatrixPlotXAxis->CenterTitle(true);
		ResponseMatrixPlotXAxis->SetTitle("log_{10}(S_{gen}/mJy)");

		TGaxis* ResponseMatrixPlotYAxis= new TGaxis(
			xmin_h,ymin_h,
			xmin_h,ymax_h,
			ymin_shifted,ymax_shifted,510,""
		);
		ResponseMatrixPlotYAxis->SetLineColor(kBlack);
		ResponseMatrixPlotYAxis->SetTitleFont(52);
		ResponseMatrixPlotYAxis->SetTitleSize(0.06);	
		ResponseMatrixPlotYAxis->SetTitleOffset(0.87);
		ResponseMatrixPlotYAxis->SetLabelFont(42);
		ResponseMatrixPlotYAxis->SetLabelSize(0.04);
		ResponseMatrixPlotYAxis->CenterTitle(true);
		ResponseMatrixPlotYAxis->SetTitle("log_{10}(S_{meas}/mJy)");

		gResponseMatrix->GetXaxis()->SetTitle("");
  	gResponseMatrix->GetXaxis()->SetLabelOffset(999);
  	gResponseMatrix->GetXaxis()->SetLabelSize(0);
		ResponseMatrixPlotXAxis->Draw();

		gResponseMatrix->GetYaxis()->SetTitle("");
  	gResponseMatrix->GetYaxis()->SetLabelOffset(999);
  	gResponseMatrix->GetYaxis()->SetLabelSize(0);
		ResponseMatrixPlotYAxis->Draw();
	
		diagFcn->Draw("lsame");
	}

	if(gRandMeanResponseMatrix){
		gStyle->SetPaintTextFormat(".2f");

		TCanvas* RandMeanResponseMatrixPlot= new TCanvas("RandMeanResponseMatrixPlot","RandMeanResponseMatrixPlot");
		RandMeanResponseMatrixPlot->cd();

		gRandMeanResponseMatrix->SetStats(0);
		gRandMeanResponseMatrix->SetMarkerColor(kWhite);
		gRandMeanResponseMatrix->SetMarkerSize(0.9);
		gRandMeanResponseMatrix->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
		gRandMeanResponseMatrix->GetYaxis()->SetTitle("log_{10}(S_{meas}/Jy)");
		gRandMeanResponseMatrix->GetZaxis()->SetRangeUser(-0.0001,1);
		gRandMeanResponseMatrix->Draw("COLZ TEXT");
	
		diagFcn->Draw("lsame");
	}

	if(gRandSigmaResponseMatrix){
		gStyle->SetPaintTextFormat(".2f");

		TCanvas* RandSigmaResponseMatrixPlot= new TCanvas("RandSigmaResponseMatrixPlot","RandSigmaResponseMatrixPlot");
		RandSigmaResponseMatrixPlot->cd();

		gRandSigmaResponseMatrix->SetStats(0);
		gRandSigmaResponseMatrix->SetMarkerColor(kWhite);
		gRandSigmaResponseMatrix->SetMarkerSize(0.9);
		gRandSigmaResponseMatrix->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
		gRandSigmaResponseMatrix->GetYaxis()->SetTitle("log_{10}(S_{meas}/Jy)");
		gRandSigmaResponseMatrix->GetZaxis()->SetRangeUser(-0.0001,1);
		gRandSigmaResponseMatrix->Draw("COLZ TEXT");
	
		diagFcn->Draw("lsame");
	}

	// - Draw flux bias model
	if(gBiasModelFcn){
		TCanvas* BiasModelPlot= new TCanvas("BiasModelPlot","BiasModelPlot");
		BiasModelPlot->cd();

		gBiasModelFcn->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
		gBiasModelFcn->GetYaxis()->SetTitle("bias");
		gBiasModelFcn->Draw("l");
	}
	

	// - Draw flux resolution model
	if(gResoModelFcn){
		TCanvas* ResoModelPlot= new TCanvas("ResoModelPlot","ResoModelPlot");
		ResoModelPlot->cd();

		gResoModelFcn->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
		gResoModelFcn->GetYaxis()->SetTitle("reso");
		gResoModelFcn->Draw("l");
	}
	
	//=======================================
	//==    SOURCE COUNTS RESO BIAS
	//=======================================
	// - Draw resolution bias model
	if(gSourceCountsResoBiasCorrFcn){
		TCanvas* ResoBiasPlot= new TCanvas("ResoBiasPlot","ResoBiasPlot");
		ResoBiasPlot->cd();

		TH2D* ResoBiasPlotBkg= new TH2D("ResoBiasPlotBkg","",100,gSourceCounts->GetXaxis()->GetXmin(),gSourceCounts->GetXaxis()->GetXmax(),100,gSourceCountsResoBiasCorrFcn->GetMinimum(),gSourceCountsResoBiasCorrFcn->GetMaximum());
		ResoBiasPlotBkg->SetStats(0);
		ResoBiasPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
		ResoBiasPlotBkg->GetYaxis()->SetTitle("reso bias corr");
		ResoBiasPlotBkg->Draw();

		gSourceCountsResoBiasCorrFcn->SetNpx(200);
		gSourceCountsResoBiasCorrFcn->Draw("l same");
	}
	
	//=======================================
	//==    RAW SOURCE COUNTS
	//=======================================
	// - Draw raw source counts
	gStyle->SetOptLogy(true);
	gStyle->SetPadRightMargin(0.05);

	TCanvas* Plot= new TCanvas("Plot","Plot");
	Plot->cd();

	TH2D* PlotBkg= new TH2D("PlotBkg","",100,gSourceCounts->GetXaxis()->GetXmin(),gSourceCounts->GetXaxis()->GetXmax(),100,1.,1.e+4);
	PlotBkg->SetTitle(0);
	PlotBkg->SetStats(0);
	PlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	PlotBkg->GetXaxis()->SetTitleSize(0.06);
	PlotBkg->GetXaxis()->SetTitleOffset(0.8);
	PlotBkg->GetYaxis()->SetTitle("#sources");
	PlotBkg->GetYaxis()->SetTitleSize(0.06);
	PlotBkg->GetYaxis()->SetTitleOffset(1.3);
	PlotBkg->Draw();

	gSourceCounts->SetMarkerColor(kBlack);
	gSourceCounts->SetLineColor(kBlack);
	gSourceCounts->SetMarkerStyle(8);
	gSourceCounts->SetMarkerSize(1.3);
	gSourceCounts->Draw("EPZ same");
	
	gUnfoldedSpectrum->SetMarkerColor(kBlack);
	gUnfoldedSpectrum->SetLineColor(kBlack);
	gUnfoldedSpectrum->SetMarkerSize(1.3);
	gUnfoldedSpectrum->SetMarkerStyle(24);
	gUnfoldedSpectrum->Draw("EPZ same");

	gFFSpectrum->SetMarkerColor(kBlack);
	gFFSpectrum->SetLineColor(kBlack);
	gFFSpectrum->SetMarkerStyle(23);
	//gFFSpectrum->Draw("L hist same");
	gFFSpectrum->Draw("hist same");

	gTrueFitModelFcn->SetLineColor(kRed);
	gTrueFitModelFcn->Draw("l same");	


	TLegend* PlotLegend= new TLegend(0.6,0.6,0.8,0.8);	
	PlotLegend->SetTextSize(0.04);
	PlotLegend->SetTextFont(52);
	//PlotLegend->SetHeader("Source Counts");
	PlotLegend->AddEntry(gSourceCounts,"measured","PL");
	PlotLegend->AddEntry(gUnfoldedSpectrum,"unfolded","PL");
	PlotLegend->AddEntry(gFFSpectrum,"fit","L");
	PlotLegend->AddEntry(gTrueFitModelFcn,"model","L");
	PlotLegend->Draw("same");

	// - Draw fit relative residuals
	gStyle->SetOptLogy(false);

	TCanvas* FitResidualPlot= new TCanvas("FitResidualPlot","FitResidualPlot");
	FitResidualPlot->cd();

	gSourceCountsFitResidualHisto->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	gSourceCountsFitResidualHisto->GetYaxis()->SetTitle("N_{ff}/N-1");
	gSourceCountsFitResidualHisto->Draw("hist same");

	//=======================================
	//==    DIFFERENTIAL SOURCE COUNTS
	//=======================================
	//- Draw normalized differential counts
	gStyle->SetOptLogy(true);

	TCanvas* DiffCountsPlot= new TCanvas("DiffCountsPlot","DiffCountsPlot");
	DiffCountsPlot->cd();

	//double DiffCountsPlotBkg_xmin= gDiffSourceCounts->GetXaxis()->GetXmin();
	//double DiffCountsPlotBkg_xmax= gDiffSourceCounts->GetXaxis()->GetXmax();
	double DiffCountsPlotBkg_xmin= -3.45;
	double DiffCountsPlotBkg_xmax= 0.7;
	double DiffCountsPlotBkg_ymin= 1.e-3;
	double DiffCountsPlotBkg_ymax= 1.e+4;

	TH2D* DiffCountsPlotBkg= new TH2D("DiffCountsPlotBkg","",100,DiffCountsPlotBkg_xmin,DiffCountsPlotBkg_xmax,100,DiffCountsPlotBkg_ymin,DiffCountsPlotBkg_ymax);
	DiffCountsPlotBkg->SetTitle(0);
	DiffCountsPlotBkg->SetStats(0);
	DiffCountsPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	DiffCountsPlotBkg->GetXaxis()->SetTitleSize(0.06);
	DiffCountsPlotBkg->GetXaxis()->SetTitleOffset(0.8);
	DiffCountsPlotBkg->GetYaxis()->SetTitle("S^{2.5} x dN/dS (Jy^{1.5}sr^{-1})");
	DiffCountsPlotBkg->GetYaxis()->SetTitleSize(0.06);
	DiffCountsPlotBkg->GetYaxis()->SetTitleOffset(1.2);
	DiffCountsPlotBkg->GetXaxis()->CenterTitle(true);
	DiffCountsPlotBkg->GetYaxis()->CenterTitle(true);
	DiffCountsPlotBkg->Draw();

	gExpDataFitAreaGraph->SetFillStyle(1001);
  gExpDataFitAreaGraph->SetFillColor(kAzure-8);
  gExpDataFitAreaGraph->Draw("F");
		
	gDiffSourceCounts->SetMarkerColor(kBlack);
	gDiffSourceCounts->SetLineColor(kBlack);
	gDiffSourceCounts->SetMarkerStyle(8);
	gDiffSourceCounts->SetMarkerSize(1.3);
	gDiffSourceCounts->Draw("EPZ same");
	
	gUnfoldedDiffSpectrum->SetMarkerColor(kBlack);
	gUnfoldedDiffSpectrum->SetLineColor(kBlack);
	gUnfoldedDiffSpectrum->SetMarkerSize(1.1);
	gUnfoldedDiffSpectrum->SetMarkerStyle(24);
	//gUnfoldedDiffSpectrum->Draw("EPZ same");

	TGraphErrors* unfoldedDiffSpectrumGraph_totErrors= new TGraphErrors;
	TGraphErrors* unfoldedDiffSpectrumGraph= new TGraphErrors;
	int npts= 0;
	int npts2= 0;
	double minBinSize= 1.e+99;
	for(int i=0;i<gUnfoldedDiffSpectrum_totErrors->GetNbinsX();i++){
		double binWidth= gUnfoldedDiffSpectrum_totErrors->GetXaxis()->GetBinWidth(i+1);
		if(binWidth<minBinSize) minBinSize= binWidth;
	}
	//double offset= minBinSize/3.;
	double offset= 0;

	for(int i=0;i<gUnfoldedDiffSpectrum_totErrors->GetNbinsX();i++){
		double lgS= gUnfoldedDiffSpectrum_totErrors->GetXaxis()->GetBinCenter(i+1);
		double binContent= gUnfoldedDiffSpectrum_totErrors->GetBinContent(i+1);
		double binError= gUnfoldedDiffSpectrum_totErrors->GetBinError(i+1);
		if(binContent>0){
			unfoldedDiffSpectrumGraph_totErrors->SetPoint(npts,lgS+offset,binContent);
			unfoldedDiffSpectrumGraph_totErrors->SetPointError(npts,0,binError);
			npts++;
		}
		
		lgS= gUnfoldedDiffSpectrum->GetXaxis()->GetBinCenter(i+1);
		binContent= gUnfoldedDiffSpectrum->GetBinContent(i+1);
		binError= gUnfoldedDiffSpectrum->GetBinError(i+1);
		if(binContent>0){
			unfoldedDiffSpectrumGraph->SetPoint(npts2,lgS+offset,binContent);
			unfoldedDiffSpectrumGraph->SetPointError(npts2,0,binError);
			npts2++;
		}
	}

	unfoldedDiffSpectrumGraph->SetMarkerColor(kBlack);
	unfoldedDiffSpectrumGraph->SetLineColor(kBlack);
	unfoldedDiffSpectrumGraph->SetMarkerSize(1.1);
	unfoldedDiffSpectrumGraph->SetMarkerStyle(24);
	//unfoldedDiffSpectrumGraph->Draw("PZ");

	unfoldedDiffSpectrumGraph_totErrors->SetMarkerColor(kBlack);
	unfoldedDiffSpectrumGraph_totErrors->SetLineColor(kBlack);
	unfoldedDiffSpectrumGraph_totErrors->SetMarkerSize(1.1);
	unfoldedDiffSpectrumGraph_totErrors->SetMarkerStyle(24);
	//unfoldedDiffSpectrumGraph_totErrors->Draw("[]");
	unfoldedDiffSpectrumGraph_totErrors->Draw("PZ");

	gUnfoldedDiffSpectrum_resoBiasCorrected->SetMarkerColor(kBlue);
	gUnfoldedDiffSpectrum_resoBiasCorrected->SetLineColor(kBlue);
	gUnfoldedDiffSpectrum_resoBiasCorrected->SetMarkerSize(1.1);
	gUnfoldedDiffSpectrum_resoBiasCorrected->SetMarkerStyle(21);
	//gUnfoldedDiffSpectrum_resoBiasCorrected->Draw("EPZ same");
	
	gFFDiffSpectrum->SetMarkerColor(kBlack);
	gFFDiffSpectrum->SetLineColor(kBlack);
	gFFDiffSpectrum->SetMarkerStyle(23);
	//gFFDiffSpectrum->Draw("L hist same");

	//Draw mJy axis
	double ymin_h= DiffCountsPlotBkg->GetYaxis()->GetXmin();
	double ymax_h= DiffCountsPlotBkg->GetYaxis()->GetXmax();
	double xmin_h= DiffCountsPlotBkg->GetXaxis()->GetXmin();
	double xmax_h= DiffCountsPlotBkg->GetXaxis()->GetXmax();
	double xmin_shifted= xmin_h + 3;
	double xmax_shifted= xmax_h + 3;

	TGaxis* DiffCountsPlotBkgXAxis= new TGaxis(
		xmin_h,ymin_h,
		xmax_h,ymin_h,
		xmin_shifted,xmax_shifted,510,"+L"
	);
	DiffCountsPlotBkgXAxis->SetLineColor(kBlack);
	DiffCountsPlotBkgXAxis->SetTitleFont(52);
	DiffCountsPlotBkgXAxis->SetTitleSize(0.06);	
	DiffCountsPlotBkgXAxis->SetTitleOffset(0.87);
	DiffCountsPlotBkgXAxis->SetLabelFont(42);
	DiffCountsPlotBkgXAxis->SetLabelSize(0.04);
	DiffCountsPlotBkgXAxis->CenterTitle(true);
	DiffCountsPlotBkgXAxis->SetTitle("log_{10}(S/mJy)");

	DiffCountsPlotBkg->GetXaxis()->SetTitle("");
  DiffCountsPlotBkg->GetXaxis()->SetLabelOffset(999);
  DiffCountsPlotBkg->GetXaxis()->SetLabelSize(0);
	DiffCountsPlotBkgXAxis->Draw();


	TLegend* DiffCountsPlotLegend= new TLegend(0.6,0.6,0.8,0.8);	
	DiffCountsPlotLegend->SetTextSize(0.04);
	DiffCountsPlotLegend->SetTextFont(52);
	DiffCountsPlotLegend->AddEntry(gDiffSourceCounts,"measured","PL");
	DiffCountsPlotLegend->AddEntry(gUnfoldedDiffSpectrum,"unfolded","PL");
	DiffCountsPlotLegend->AddEntry(gExpDataFitAreaGraph,Form("1.4 GHz data shifted to %1.0f MHz",gDataFrequency*1000),"F");
	//DiffCountsPlotLegend->AddEntry(gFFDiffSpectrum,"fit","L");
	//DiffCountsPlotLegend->AddEntry(trueDiffModelCounts,"model","L");
	DiffCountsPlotLegend->Draw("same");


	// - Compute difference between observed counts and model
	gStyle->SetOptLogy(0);

	TCanvas* ObsModelCountDiffPlot= new TCanvas("ObsModelCountDiffPlot","ObsModelCountDiffPlot");
	ObsModelCountDiffPlot->cd();

	TH2D* ObsModelCountDiffPlotBkg= new TH2D("ObsModelCountDiffPlot","",100,gDiffSourceCounts->GetXaxis()->GetXmin(),gDiffSourceCounts->GetXaxis()->GetXmax(),100,-1,1);
	ObsModelCountDiffPlotBkg->SetTitle(0);
	ObsModelCountDiffPlotBkg->SetStats(0);
	ObsModelCountDiffPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	ObsModelCountDiffPlotBkg->GetXaxis()->SetTitleSize(0.06);
	ObsModelCountDiffPlotBkg->GetXaxis()->SetTitleOffset(0.8);
	ObsModelCountDiffPlotBkg->GetYaxis()->SetTitle("dN/dS_{obs}/dN/dS_{model}-1");
	ObsModelCountDiffPlotBkg->GetYaxis()->SetTitleSize(0.06);
	ObsModelCountDiffPlotBkg->GetYaxis()->SetTitleOffset(1.3);
	ObsModelCountDiffPlotBkg->Draw();

	TGraph* sourceCountsObsToModelDiffGraph_min= new TGraph;
	TGraph* sourceCountsObsToModelDiffGraph_mean= new TGraph;
	int np= 0;

	for(int i=0;i<gUnfoldedDiffSpectrum->GetNbinsX();i++){
		double binCenter= gUnfoldedDiffSpectrum->GetBinCenter(i+1);
		double counts_obs= gUnfoldedDiffSpectrum->GetBinContent(i+1);
		double counts_model_min= gExpDataFitAreaGraph_min->Eval(binCenter);
		double counts_model_mean= gExpDataFitAreaGraph_mean->Eval(binCenter);
				
		if(counts_obs<=0 || counts_model_min<=0 || counts_model_mean<=0) continue;
		
		double pull_mean= counts_obs/counts_model_mean-1;
		double pull_min= counts_obs/counts_model_min-1;
		sourceCountsObsToModelDiffGraph_mean->SetPoint(np,binCenter,pull_mean);
		sourceCountsObsToModelDiffGraph_min->SetPoint(np,binCenter,pull_min);
		np++;
	}//end loop bins
	
	sourceCountsObsToModelDiffGraph_min->SetMarkerColor(kBlack);
	sourceCountsObsToModelDiffGraph_min->SetLineColor(kBlack);
	sourceCountsObsToModelDiffGraph_min->SetMarkerStyle(8);
	sourceCountsObsToModelDiffGraph_min->SetMarkerSize(1.3);
	sourceCountsObsToModelDiffGraph_min->Draw("PLZ");

	sourceCountsObsToModelDiffGraph_mean->SetMarkerColor(kRed);
	sourceCountsObsToModelDiffGraph_mean->SetLineColor(kRed);
	sourceCountsObsToModelDiffGraph_mean->SetMarkerStyle(21);
	sourceCountsObsToModelDiffGraph_mean->SetMarkerSize(1.1);
	sourceCountsObsToModelDiffGraph_mean->Draw("PLZ");


}//close Draw()


void PrintSourceCountTable()
{
	std::string delimiter= "&";	
	
	cout<<"== SOURCE COUNTS TABLE =="<<endl;
	for(int i=0;i<gSourceCounts->GetNbinsX();i++){
		double lgS= gSourceCounts->GetBinCenter(i+1);
		double dlgS= gSourceCounts->GetBinWidth(i+1);	
		double lgSmin= gSourceCounts->GetBinLowEdge(i+1);	
		double lgSmax= lgSmin + dlgS;
		double S= pow(10,lgS)*1.e+3;//mJy
		double Smin= pow(10,lgSmin);
		double Smax= pow(10,lgSmax);
		double dS= (Smax-Smin)*1.e+3;//mJy

		double N= gSourceCounts->GetBinContent(i+1);
		double N_unfold= gUnfoldedSpectrum->GetBinContent(i+1);
		double counts= gDiffSourceCounts->GetBinContent(i+1);
		double counts_unfold= gUnfoldedDiffSpectrum->GetBinContent(i+1);
		double countsErr_unfold= gUnfoldedDiffSpectrum->GetBinError(i+1);
		double countsTotErr_unfold= gUnfoldedDiffSpectrum_totErrors->GetBinError(i+1);
		

		if(N<=0) continue;
		
		printf("%.2f %s %.2f %s %d %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f\n",dS,delimiter.c_str(),S,delimiter.c_str(),(int)(N),delimiter.c_str(),N_unfold,delimiter.c_str(),counts,delimiter.c_str(),counts_unfold,delimiter.c_str(),countsErr_unfold,delimiter.c_str(),countsTotErr_unfold);
	}
	cout<<"========================="<<endl;


}//close PrintSourceCountTable()

int ReadData(std::string filename)
{
	//Open files
	inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filename<<"!");
		#endif
		return -1;
	}

	//Get access to source count histo
	gSourceCounts= (TH1D*)inputFile->Get(histoName.c_str());
	if(!gSourceCounts){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get access to source counts histo with name "<<histoName<<" in file "<<filename<<"!");	
		#endif
		return -1;
	}

	//Set rec histogram binning
	double binLowEdge= 0;
	for(int i=0;i<gSourceCounts->GetNbinsX();i++){
		binLowEdge= gSourceCounts->GetBinLowEdge(i+1);
		gRecBins.push_back(binLowEdge);
	}
	gRecBins.push_back( binLowEdge + gSourceCounts->GetBinWidth( gSourceCounts->GetNbinsX() ) );

	//Set true histogram binning (it should be equal or larger than RecBins)
	//Assume here same rec binning for simplicity
	if(gTrueBins.empty()){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Set true bins to rec bins...");	
		#endif
		gTrueBins.assign(gRecBins.begin(),gRecBins.end());
	}

	std::stringstream ss;
	ss<<"RecBins(";
	for(size_t i=0;i<gRecBins.size()-1;i++){
		ss<<gRecBins[i]<<",";
	}
	ss<<gRecBins[gRecBins.size()-1]<<")";
	
	#ifdef LOGGING_ENABLED
		INFO_LOG(ss.str());	
	#endif

	ss.str("");
	ss<<"TrueBins(";
	for(size_t i=0;i<gTrueBins.size()-1;i++){
		ss<<gTrueBins[i]<<",";
	}
	ss<<gTrueBins[gTrueBins.size()-1]<<")";
	
	#ifdef LOGGING_ENABLED
		INFO_LOG(ss.str());	
	#endif

	return 0;

}//close ReadSourceData()


void Save()
{
	// - Save objects to file
	if(outputFile && outputFile->IsOpen())
	{
		outputFile->cd();		
		
		//Write source counts
		if(gSourceCounts) gSourceCounts->Write();		
		if(gDiffSourceCounts) gDiffSourceCounts->Write();
		if(gUnfoldedDiffSpectrum) gUnfoldedDiffSpectrum->Write();
		if(gUnfoldedDiffSpectrum_totErrors) gUnfoldedDiffSpectrum_totErrors->Write();
		if(gUnfoldedDiffSpectrum_resoBiasCorrected) gUnfoldedDiffSpectrum_resoBiasCorrected->Write();
		if(gFFDiffSpectrum) gFFDiffSpectrum->Write();

		//Write response matrix
		if(gResponseMatrix) gResponseMatrix->Write();
		if(gResponseModelFcn) {
			gResponseModelFcn->SetNameTitle("responseModelFcn","responseModelFcn");
			gResponseModelFcn->Write();
		}

		//Close file
		outputFile->Close();
	}

	
}//close Save()

void SetStyle()
{
	TStyle* myStyle= new TStyle("myStyle","myStyle");

	//## CANVAS & PAD
	myStyle->SetCanvasDefH(700); 
  myStyle->SetCanvasDefW(700); 
	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
	myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.16);
  myStyle->SetPadRightMargin(0.1);

	//## TITLE
	myStyle->SetOptTitle(0);
	myStyle->SetTitleX(0.1f);
	myStyle->SetTitleW(0.8f);                 
  myStyle->SetTitleXOffset(0.8);
  myStyle->SetTitleYOffset(1.1);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(0);//border size of Title PavelLabel
	myStyle->SetTitleSize(0.06,"X");
  myStyle->SetTitleSize(0.06,"Y");
	myStyle->SetTitleSize(0.06,"Z");
	
	//## STAT
	myStyle->SetOptStat("eMR");
	//myStyle->SetOptStat(1);
  myStyle->SetStatColor(0);
	myStyle->SetStatY(0.975);                
  myStyle->SetStatX(0.95);                
  myStyle->SetStatW(0.35);//0.2                
  myStyle->SetStatH(0.10);//0.15
  myStyle->SetStatBorderSize(1);

	myStyle->SetTitleFont(52,"X");
  myStyle->SetTitleFont(52,"Y");
  myStyle->SetTitleFont(52,"Z");
  myStyle->SetLabelFont(42,"X");
  myStyle->SetLabelFont(42,"Y");
  myStyle->SetLabelFont(42,"Z");   
	
	//## OTHER
  myStyle->SetOptFit(1);
	myStyle->SetOptLogx(0);
	myStyle->SetOptLogy(0);
  //myStyle->SetPalette(1,0);
  myStyle->SetMarkerStyle(8);
  myStyle->SetMarkerSize(0.6);
  myStyle->SetFuncWidth(1.); 
  myStyle->SetErrorX(0.);

	myStyle->SetNumberContours(999);
	//myStyle->SetPalette(kColorPrintableOnGrey);
	//myStyle->SetPalette(kTemperatureMap);
	//myStyle->SetPalette(kBlueGreenYellow);
	myStyle->SetPalette(kRainBow);
	
	gROOT->SetStyle("myStyle");
	gStyle= myStyle;
	myStyle->cd();

}//close SetStyle()


std::string GetStringLogLevel(int verbosity)
{
	std::string slevel= "";
	if(verbosity<=0) slevel= "FATAL";
	else if(verbosity==1) slevel= "FATAL";
	else if(verbosity==2) slevel= "ERROR";
	else if(verbosity==3) slevel= "WARN";
	else if(verbosity==4) slevel= "INFO";
	else if(verbosity>=5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

