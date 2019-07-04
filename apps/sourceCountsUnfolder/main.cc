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
	cout<<"-S, --smodel=[SPECTRUM_MODEL] \t Spectrum model used for response matrix creation (flat,powerlaw,brokenpowerlaw)"<<endl;
	cout<<"-s, --smodel-fit=[SPECTRUM_MODEL] \t Spectrum model used as fit model in unfolding (flat,powerlaw,brokenpowerlaw)"<<endl;
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
double gSpectrumModelLgSBreak_min= -3.0;
double gSpectrumModelLgSBreak_max= -2.6;
double gSpectrumModelLgSCutoff= -0.67;
double gSpectrumModelLgSCutoff_min= -0.8;
double gSpectrumModelLgSCutoff_max= -0.4;
std::vector<double> gEfficiencyPars {9.86885e-01,-2.85922e+00,4.14620e-01};
std::vector<double> gFluxBiasPars {1.97524e-06,3.82100e+00,-4.54680e-04};
std::vector<double> gFluxResoPars {6.60898e-04,2.04734e+00};
int gBiasModel= eCONST_RESO;
int gResoModel= eCONST_BIAS;
int gEfficiencyModel= eCONST_EFF;
int gSpectrumModel= eFlat;
int gSpectrumModel_fit= eBrokenPowerLaws;
double gConstEfficiency= 1;
double gConstBias= 0;
double gConstReso= 0.2;
bool gUseFitRange= true;
double gLgSMin_fit= -3.6;
double gLgSMax_fit= 1.5;
bool gComputeUncertainties= false;
int gNRandSamples= 100;
	

//Globar vars
TFile* inputFile= 0;
TFile* outputFile= 0;
TApplication* app= 0;
TH1D* gSourceCounts= 0;
std::vector<double> gRecBins;
std::vector<double> gTrueBins;
TH2D* gResponseMatrix= 0;
TF2* gResponseModelFcn= 0;
ForwardFolder* gForwardFolder= 0;
TH1D* gUnfoldedSpectrum= 0;
TH1D* gFFSpectrum= 0;


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadData(std::string filename);
int SetResponseMatrix();
int UnfoldData();
void Save();
int Init();
void Draw();
void SetStyle();
void ClearData();

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

	//=======================
	//== DRAW PLOTS
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Drawing data ...");
	#endif
	Draw();

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

	while((c = getopt_long(argc, argv, "hi:H:o:v:Ig:G:l:c:C:e:b:r:E:B:R:S:s:",options_tab, &option_index)) != -1) {
    
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
				else gResoModel= eCONST_EFF;
				break;	
			}
			case 'S':	
			{
				std::string spectrumModel_str= std::string(optarg);
				if(spectrumModel_str=="flat") gSpectrumModel= eFlat;
				else if(spectrumModel_str=="powerlaw") gSpectrumModel= ePowerLaw;
				else if(spectrumModel_str=="brokenpowerlaw") gSpectrumModel= eBrokenPowerLaws;
				else if(spectrumModel_str=="smoothbrokenpowerlaw") gSpectrumModel= eSmoothBrokenPowerLaws;
				else gSpectrumModel= ePowerLaw;
				break;	
			}
			case 's':	
			{
				std::string spectrumModel_str= std::string(optarg);
				if(spectrumModel_str=="flat") gSpectrumModel_fit= eFlat;
				else if(spectrumModel_str=="powerlaw") gSpectrumModel_fit= ePowerLaw;
				else if(spectrumModel_str=="brokenpowerlaw") gSpectrumModel_fit= eBrokenPowerLaws;
				else if(spectrumModel_str=="smoothbrokenpowerlaw") gSpectrumModel_fit= eSmoothBrokenPowerLaws;
				else gSpectrumModel_fit= ePowerLaw;
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

	//Check reso pars
	if(gFluxResoPars.empty() || gFluxResoPars.size()!=2){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid reso pars given (2 pars required)!");
		#endif
		return -1;
	}

	//Check efficiency pars
	if(gEfficiencyPars.empty() || gEfficiencyPars.size()!=3){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid efficiency pars given (3 pars required)!");
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
	double gamma1= gSpectrumModelIndex1;
	double gamma2= gSpectrumModelIndex2;
	double gamma3= gSpectrumModelIndex3;
	double gamma_min= 0;
	double gamma_max= 10;
	double lgS_break= gSpectrumModelLgSBreak;
	double lgS_break_min= gSpectrumModelLgSBreak_min;
	double lgS_break_max= gSpectrumModelLgSBreak_max;
	double lgS_cutoff= gSpectrumModelLgSCutoff;
	double lgS_cutoff_min= gSpectrumModelLgSCutoff_min;
	double lgS_cutoff_max= gSpectrumModelLgSCutoff_max;
	double lgS_min= -5;
	double lgS_max= 10;
	FitPar normPar= FitPar("Norm",1);
	FitPar gammaPar1= FitPar("Gamma1",gamma1,gamma_min,gamma_max);
	FitPar gammaPar2= FitPar("Gamma2",gamma2,gamma_min,gamma_max);
	FitPar gammaPar3= FitPar("Gamma3",gamma3,gamma_min,gamma_max);
	FitPar breakPar= FitPar("Break",lgS_break,lgS_break_min,lgS_break_max);
	FitPar cutoffPar= FitPar("Cutoff",lgS_cutoff,lgS_cutoff_min,lgS_cutoff_max);


	SpectrumPars* spectrumPars= 0;
	if(gSpectrumModel_fit==ePowerLaw){
		spectrumPars= new PowerLawPars(normPar,gammaPar2);
	}
	else if(gSpectrumModel_fit==eFlat){
		spectrumPars= new FlatSpectrumPars(normPar);
	}
	else if(gSpectrumModel_fit==eBrokenPowerLaws){
		spectrumPars= new BrokenPowerLawsPars(normPar,gammaPar1, gammaPar2, gammaPar3,breakPar,cutoffPar);
	}

	//Fix spectrum pars?
	//spectrumPars.FixPar("Gamma1",gamma1);
	//spectrumPars.FixPar("Break",lgS_break);
	

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
	gFFSpectrum= gForwardFolder->GetForwardFoldedSpectrum();
		
	
	return 0;

}//close SelectSources()



int SetResponseMatrix()
{
	//====================================================
	//==         SET RESPONSE MODEL PARS
	//=====================================================		
	//- Set spectrum pars
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


	SpectrumPars* spectrumPars= 0;
	if(gSpectrumModel==ePowerLaw){
		spectrumPars= new PowerLawPars(normPar,gammaPar2);
	}
	else if(gSpectrumModel==eFlat){
		spectrumPars= new FlatSpectrumPars(normPar);
	}
	else if(gSpectrumModel==eBrokenPowerLaws){
		spectrumPars= new BrokenPowerLawsPars(normPar,gammaPar1, gammaPar2, gammaPar3,breakPar,cutoffPar);
	}
	
	
	// - Set bias pars
	BiasPars* biasPars= 0;
	if(gBiasModel==eCONST_BIAS) biasPars= new ConstBiasPars(gConstBias);
	else if(gBiasModel==eEXP_BIAS) biasPars= new ExpBiasPars(gFluxBiasPars[0],gFluxBiasPars[1],gFluxBiasPars[2]);
	else if(gBiasModel==eEXP_STEP_BIAS) biasPars= new ExpStepBiasPars(gFluxBiasPars[0],gFluxBiasPars[1],gFluxBiasPars[2],gFluxBiasPars[3]);

	
	// - Set resolution pars
	ResoPars* resoPars= 0;
	if(gResoModel==eCONST_RESO) resoPars= new ConstResoPars(gConstReso);
	else if(gResoModel==eEXP_RESO) resoPars= new ExpResoPars(gFluxResoPars[0],gFluxResoPars[1]);

	
	// - Set efficiency pars
	EfficiencyPars* effPars= 0;
	if(gEfficiencyModel==eCONST_EFF) effPars= new ConstEfficiencyPars(gConstEfficiency);
	else if(gEfficiencyModel==eSIGMOID_EFF) effPars= new SigmoidEfficiencyPars(gEfficiencyPars[0],gEfficiencyPars[1],gEfficiencyPars[2]);

	
	
	//====================================================
	//==         BUILD RESPONSE MATRIX
	//=====================================================
	double lgSmin_true= gTrueBins[0];
	double lgSmax_true= gTrueBins[gTrueBins.size()-1];
	double lgSmin_rec= gRecBins[0];
	double lgSmax_rec= gRecBins[gRecBins.size()-1];

	gResponseModelFcn= SpectrumUtils::ComputeResponseModel(
		*spectrumPars,
		*biasPars,
		*resoPars,
		*effPars,
		lgSmin_true,lgSmax_true,
		lgSmin_rec,lgSmax_rec
	);

	
	gResponseMatrix= SpectrumUtils::ComputeParametricResponse(
		*spectrumPars,
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
	

	return 0;

}//close SetResponseMatrix()


void Draw()
{
	// - Draw response function 2D
	TCanvas* ResponseModelPlot= new TCanvas("ResponseModelPlot","ResponseModelPlot");
	ResponseModelPlot->cd();

	if(gResponseModelFcn){
		gResponseModelFcn->Draw("COLZ");
	}

	// - Draw response matrix
	TCanvas* ResponseMatrixPlot= new TCanvas("ResponseMatrixPlot","ResponseMatrixPlot");
	ResponseMatrixPlot->cd();

	if(gResponseMatrix){
		gResponseMatrix->SetStats(0);
		gResponseMatrix->Draw("COLZ");
	}

	TF1* diagFcn= new TF1("diagFcn","x",gTrueBins[0],gTrueBins[gTrueBins.size()-1]);
	diagFcn->Draw("lsame");

	
	//--> Results Plot
	TCanvas* Plot= new TCanvas("Plot","Plot");
	Plot->cd();

	gSourceCounts->SetTitle(0);
	gSourceCounts->SetStats(0);
	gSourceCounts->GetXaxis()->SetTitle("lg(S/Jy)");
	gSourceCounts->GetXaxis()->SetTitleSize(0.06);
	gSourceCounts->GetXaxis()->SetTitleOffset(0.8);
	gSourceCounts->GetYaxis()->SetTitle("nentries");
	gSourceCounts->GetYaxis()->SetTitleSize(0.06);
	gSourceCounts->GetYaxis()->SetTitleOffset(1.3);
	gSourceCounts->SetMarkerColor(kBlack);
	gSourceCounts->SetLineColor(kBlack);
	gSourceCounts->SetMarkerStyle(8);
	gSourceCounts->Draw("ep");
	
	gUnfoldedSpectrum->SetMarkerColor(kRed);
	gUnfoldedSpectrum->SetLineColor(kRed);
	gUnfoldedSpectrum->SetMarkerStyle(21);
	gUnfoldedSpectrum->Draw("ep same");
	
	gFFSpectrum->SetMarkerColor(kGreen+1);
	gFFSpectrum->SetLineColor(kGreen+1);
	gFFSpectrum->SetMarkerStyle(23);
	gFFSpectrum->Draw("ep same");

	TLegend* PlotLegend= new TLegend(0.6,0.6,0.8,0.8);	
	PlotLegend->SetTextSize(0.045);
	PlotLegend->SetTextFont(52);
	PlotLegend->AddEntry(gSourceCounts,"rec","PL");
	PlotLegend->AddEntry(gUnfoldedSpectrum,"unfolded","PL");
	PlotLegend->AddEntry(gFFSpectrum,"fitted","PL");
	PlotLegend->Draw("same");

}//close Draw()

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
	gTrueBins.assign(gRecBins.begin(),gRecBins.end());
	

	std::stringstream ss;
	ss<<"RecBins(";
	for(size_t i=0;i<gRecBins.size()-1;i++){
		ss<<gRecBins[i]<<",";
	}
	ss<<gRecBins[gRecBins.size()-1]<<")";
	
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
	myStyle->SetPalette(55);

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

