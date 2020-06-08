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
/**
* @file SpectrumUtils.cc
* @class SpectrumUtils
* @brief Spectrum utilities
* 
* @author S. Riggi
* @date 03/07/2019
*/


#include <SpectrumUtils.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>


#include <Math/GSLRndmEngines.h>
#include <TDecompChol.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>

using namespace std;
using namespace ROOT::Math;
using namespace Caesar;

//==========================================
//      STATIC VARS
//==========================================
int ResoPars::nPars;
int BiasPars::nPars;
int EfficiencyPars::nPars;
int SpectrumPars::nPars;
int SourceCountsResoBiasPars::nPars;


//==========================================
//      SPECTRUM UTILS
//==========================================
SpectrumUtils::SpectrumUtils()
{
	
}//close constructor


SpectrumUtils::~SpectrumUtils()
{
	
}//close destructor


TH1D* SpectrumUtils::GetFoldedSpectrum(TH1D* trueSpectrum,TH2D* responseMatrix)
{
	//## Forward-folding of the input histogram with the response matrix
	//Check input args
	if(!trueSpectrum || !responseMatrix){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input spectrum and/or response matrix!");
		#endif	
		return nullptr;
	}

	int nTrueBins= responseMatrix->GetNbinsX();
	int nRecBins= responseMatrix->GetNbinsY();

	if(trueSpectrum->GetNbinsX() != nTrueBins){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Number of bins of passed histogram differs from that of the response matrix!");
		#endif	
		return nullptr;
	}

	//Check consistency between the true input spectrum and response matrix (xaxis)
	if(!CheckSpectrumBinnings(trueSpectrum->GetXaxis(),responseMatrix->GetXaxis())){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Bin mismatch between response matrix and true spectrum to be folded!");
		#endif
		return nullptr;
	}

	//Create the folded spectrum
	bool isVariableBin= responseMatrix->GetYaxis()->IsVariableBinSize();
	TH1D* foldSpectrum= 0;
	if(isVariableBin){
		const double* ybins= responseMatrix->GetYaxis()->GetXbins()->GetArray();
		foldSpectrum= new TH1D("foldSpectrum","foldSpectrum",nRecBins,ybins);
	}
	else{
		foldSpectrum= new TH1D("foldSpectrum","foldSpectrum",nRecBins,responseMatrix->GetYaxis()->GetXmin(),responseMatrix->GetYaxis()->GetXmax());
	}
	foldSpectrum->Reset();

	//Fill the folded spectrum
	for(int i=0;i<nRecBins;i++)
	{
		double nFolded= 0.;		

		for(int j=0;j<nTrueBins;j++){
			double binCenterX= responseMatrix->GetXaxis()->GetBinCenter(j+1);
			double w= responseMatrix->GetBinContent(j+1,i+1);

			//Find corresponding bin in true spectrum 
			int binId= trueSpectrum->FindBin(binCenterX);
			double nTrue= trueSpectrum->GetBinContent(binId);
			
			nFolded+= w*nTrue;
		}//end loop rec bins

		foldSpectrum->SetBinContent(i+1,nFolded);
		foldSpectrum->SetBinError(i+1,0.);

	}//end loop true bins

	return foldSpectrum;

}//close GetFoldedSpectrum()


TH1D* SpectrumUtils::GetFluctuatedSpectrum(TH1D* inputSpectrum)
{
	//## Takes an histogram as input and returns an histogram with entries
	//## fluctuated according to a multinomial distribution 
	if(!inputSpectrum){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null pointer to input histogram!");
		#endif
		return nullptr;
	}
	
	if(inputSpectrum->GetEntries()<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Input histogram is empty...returning same histo!");
		#endif
		return ((TH1D*)inputSpectrum->Clone("fluctuatedSpectrum"));
	}

	//### INIT RANDOM GENERATOR
	//## Generate a suitable seed
	unsigned long int seed1= time(0);
	unsigned long int seed2= clock();
	timespec timeStruct;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timeStruct);
	
	unsigned long int seed3= timeStruct.tv_nsec;
	
	ROOT::Math::GSLRandomEngine randomEngine;
	randomEngine.Initialize(); 
	randomEngine.SetSeed(seed1+seed2+seed3);

	//## Compute vector with probabilities
	//## e.g. given the input histogram, regarded as the expected counts, p= n_i/n_tot 
	std::vector<double> multinomialProbability;
	multinomialProbability.clear();
	multinomialProbability.assign(inputSpectrum->GetNbinsX(),0.);//init

	std::vector<double> poissonNEvents;
	poissonNEvents.clear();
	poissonNEvents.assign(inputSpectrum->GetNbinsX(),0.);//init
	double ntot= inputSpectrum->Integral();
	for(int i=0;i<inputSpectrum->GetNbinsX();i++){
		double n= inputSpectrum->GetBinContent(i+1);
		multinomialProbability[i]= n/ntot;

		poissonNEvents[i]= randomEngine.Poisson(n);
	}

	std::vector<unsigned int> fluctuatedEntries = randomEngine.Multinomial(inputSpectrum->GetEntries(), multinomialProbability);
	
	TH1D* fluctuatedSpectrum= (TH1D*)inputSpectrum->Clone("fluctuatedSpectrum");
	fluctuatedSpectrum->Reset();
	
	for(int i=0;i<fluctuatedSpectrum->GetNbinsX();i++){
		fluctuatedSpectrum->SetBinContent(i+1,fluctuatedEntries[i]);
		fluctuatedSpectrum->SetBinError(i+1,sqrt(fluctuatedEntries[i]));
	}
	
	return fluctuatedSpectrum;
	
}//close GetFluctuatedSpectrum()


int SpectrumUtils::CopySpectrumContent(TH1D* h1,TH1D* h2)
{
	//Check inputs
	if(!h1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr histogram given!");
		#endif
	  return -1;
	}

	h2->Reset();
	for(int i=0;i<h1->GetNbinsX();i++){
		h2->SetBinContent(i+1,h1->GetBinContent(i+1));
		h2->SetBinError(i+1,h1->GetBinError(i+1));
	}
	
	return 0;

}//close CopyHistoContent()


bool SpectrumUtils::CheckSpectrumBinnings(const TAxis* a1, const TAxis* a2)
{
	//Axis 1 shall be contained inside the Axis 2 AND binning shall match
	int Nbins_1= a1->GetNbins();
	int Nbins_2= a2->GetNbins();
	float xmin_1= a1->GetXmin();
	float xmax_1= a1->GetXmax();
	float xmin_2= a2->GetXmin();
	float xmax_2= a2->GetXmax();

	if(xmin_1<xmin_2 || xmax_1>xmax_2){
		#ifdef LOGGING_ENABLED
			WARN_LOG("First histo axis (["<<xmin_1<<","<<xmax_1<<"]) is not contained inside the second histo (["<<xmin_2<<","<<xmax_2<<"])!");
		#endif
		return false;
	}

	//Check binnings
	bool isFailed= false;
	for(int i=0;i<Nbins_1;i++){
		double binCenter_1= a1->GetBinCenter(i+1);
		double binLowEdge_1= a1->GetBinLowEdge(i+1);
		double binUpEdge_1= binLowEdge_1 + a1->GetBinWidth(i+1);
		
		//Find bin id in axis 2
		int binId_2= a2->FindBin(binCenter_1);
		double binCenter_2= a2->GetBinCenter(binId_2);
		double binLowEdge_2= a2->GetBinLowEdge(binId_2);
		double binUpEdge_2= binLowEdge_2+a2->GetBinWidth(binId_2);
		
		if( !TMath::AreEqualRel(binLowEdge_1,binLowEdge_2,1E-10) || !TMath::AreEqualRel(binUpEdge_1,binUpEdge_2,1E-10)){
			isFailed= true;
			break;
		}
	}//end loop internal histo bins

	if(isFailed){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Bin mismatch between the two histograms!");
		#endif
		return false;
	}

	return true;

}//close CheckSpectrumBinnings()


int SpectrumUtils::GetModelSpectrum(TH1D& spectrum,TF1* spectrumModel,bool integrateBins)
{
	//Check inputs
	if(!spectrumModel) return -1;

	//Fill histo from function
	spectrum.Reset();
	if(integrateBins){
		for(int i=0;i<spectrum.GetNbinsX();i++){
			double xmin= spectrum.GetBinLowEdge(i+1);
			double xmax= xmin+spectrum.GetBinWidth(i+1);
			double w= spectrumModel->Integral(xmin,xmax)/(xmax-xmin);
			if(!TMath::IsNaN(w) && fabs(w)!=TMath::Infinity() ){
				spectrum.SetBinContent(i+1,w);
				spectrum.SetBinError(i+1,0);
			}
		}
	}//close if	
	else{
		for(int i=0;i<spectrum.GetNbinsX();i++){
			double x= spectrum.GetBinCenter(i+1);
			double w= spectrumModel->Eval(x);
			if(!TMath::IsNaN(w) && fabs(w)!=TMath::Infinity() ){
				spectrum.SetBinContent(i+1,w);
				spectrum.SetBinError(i+1,0);
			}
		}
	}//close else

	return 0;

}//close GetModelSpectrum()


TF1* SpectrumUtils::ComputeSpectrumModel(SpectrumPars& pars,double xmin,double xmax,int npts)
{
	//Create model function
	TF1* SpectrumModel= 0;
	int nSpectrumPars= pars.GetNPars();
	int spectrumModel= pars.GetModel();
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("nSpectrumPars="<<nSpectrumPars<<", spectrumModel="<<spectrumModel);
	#endif

	if(spectrumModel==ePowerLaw){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::PowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==eFlat){
		SpectrumModel= new TF1("SpectrumModel","[0]",xmin,xmax);
	}
	else if(spectrumModel==eBrokenPowerLaws){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::BrokenPowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==eTwoBrokenPowerLaws){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::TwoBrokenPowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::SmoothCutoffPowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==ePowerLawWithCutoff){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::PowerLawWithCutoffSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==ePol3){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::Pol3Spectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==ePol6){
		SpectrumModel= new TF1("SpectrumModel",SpectrumUtils::Pol6Spectrum,xmin,xmax,nSpectrumPars);
	}	
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid model selected!");
		#endif
		return nullptr;
	}

	//Set model parameters
	SpectrumModel->SetNpx(npts);
	for(int i=0;i<nSpectrumPars;i++){
		double parValue= pars.GetPar(i)->GetValue();
		SpectrumModel->SetParameter(i,parValue);
	}	
	
	//Compute integral and normalize to 1
	double integral= pars.GetIntegral(xmin,xmax);
	if(spectrumModel==eFlat){
		integral= (xmax-xmin);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		integral= SpectrumModel->Integral(xmin,xmax);
	}
	else if(spectrumModel==ePowerLawWithCutoff){
		integral= SpectrumModel->Integral(xmin,xmax);
	}
	else if(spectrumModel==ePol3){
		integral= SpectrumModel->Integral(xmin,xmax);
	}	
	else if(spectrumModel==ePol6){
		integral= SpectrumModel->Integral(xmin,xmax);
	}

	//Check if integral is 0 or nan
	if(!std::isfinite(integral) || std::isnan(integral)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Spectrum model integral is inf or nan!");
		#endif
		return nullptr;
	}

	SpectrumModel->SetParameter(0,1./integral);

	return SpectrumModel;

}//close ComputeSpectrumModel()


TF1* SpectrumUtils::ComputeEfficiencyModel(EfficiencyPars& pars,double xmin,double xmax,int npts)
{
	//Create model function
	TF1* efficiencyModel= 0;
	int nEffPars= pars.GetNPars();
	int effModel= pars.GetModel();
	std::vector<double> effParList= pars.GetPars();
	
	if(effModel==eCONST_EFF){
		efficiencyModel= new TF1("EfficiencyModel","[0]",xmin,xmax);
	}
	else if(effModel==eSIGMOID_EFF){
		efficiencyModel= new TF1("EfficiencyModel",SpectrumUtils::EfficiencyModel,xmin,xmax,nEffPars);
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid model selected!");
		#endif
		return nullptr;
	}

	//Set model parameters
	efficiencyModel->SetNpx(npts);
	efficiencyModel->SetParameters(effParList.data());	
	
	return efficiencyModel;

}//close ComputeEfficiencyModel()


TF1* SpectrumUtils::ComputeBiasModel(BiasPars& biasPars,double xmin,double xmax,int npts)
{
	//Create model function
	TF1* biasModelFcn= 0;
	int nBiasPars= biasPars.GetNPars();
	int biasModel= biasPars.GetModel();
	std::vector<double> biasParList= biasPars.GetPars();

	//Compute bias
	if(biasModel==eCONST_BIAS){
		biasModelFcn= new TF1("BiasModel","[0]",xmin,xmax);
	}	
	else if(biasModel==eEXP_BIAS){
		biasModelFcn= new TF1("BiasModel",SpectrumUtils::ExpBiasModel,xmin,xmax,nBiasPars);
	}	
	else if(biasModel==eEXP_STEP_BIAS){
		biasModelFcn= new TF1("BiasModel",SpectrumUtils::StepExpBiasModel,xmin,xmax,nBiasPars);
	}

	//Set model parameters
	biasModelFcn->SetNpx(npts);
	biasModelFcn->SetParameters(biasParList.data());	

	return biasModelFcn;

}//close ComputeBiasModel()


TF1* SpectrumUtils::ComputeResoModel(ResoPars& resoPars,double xmin,double xmax,int npts)
{
	//Create model function
	TF1* resoModelFcn= 0;
	int nResoPars= resoPars.GetNPars();
	int resoModel= resoPars.GetModel();
	std::vector<double> resoParList= resoPars.GetPars();

	//Compute bias
	if(resoModel==eCONST_RESO){
		resoModelFcn= new TF1("ResoModel","[0]",xmin,xmax);
	}	
	else if(resoModel==eEXP_RESO){
		resoModelFcn= new TF1("ResoModel",SpectrumUtils::ExpResolutionModel,xmin,xmax,nResoPars);
	}	
	else if(resoModel==eEXP_STEP_RESO){
		resoModelFcn= new TF1("ResoModel",SpectrumUtils::ExpStepResolutionModel,xmin,xmax,nResoPars);
	}

	//Set model parameters
	resoModelFcn->SetNpx(npts);
	resoModelFcn->SetParameters(resoParList.data());	

	return resoModelFcn;

}//close ComputeResoModel()


TF2* SpectrumUtils::ComputeResponseModel(SpectrumPars& spectrumPars,BiasPars& biasPars,ResoPars& resoPars,EfficiencyPars& effPars,double xmin,double xmax,double ymin,double ymax,int npts)
{
	//Get spectrum pars
	int nSpectrumPars= spectrumPars.GetNPars();
	int spectrumModel= spectrumPars.GetModel();
	FitPars spectrumParList= spectrumPars.GetPars();

	//Get bias pars
	int nBiasPars= biasPars.GetNPars();
	int biasModel= biasPars.GetModel();
	std::vector<double> biasParList= biasPars.GetPars();

	//Get reso pars
	int nResoPars= resoPars.GetNPars();
	int resoModel= resoPars.GetModel();
	std::vector<double> resoParList= resoPars.GetPars();

	//Get trigger pars
	int nEffPars= effPars.GetNPars();
	int effModel= effPars.GetModel();
	std::vector<double> effParList= effPars.GetPars();

	int nTotPars= (nBiasPars+1) + (nResoPars+1) + (nEffPars+1) + (nSpectrumPars+1);

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("nBiasPars="<<nBiasPars<<", nResoPars="<<nResoPars<<", nEffPars="<<nEffPars<<", nSpectrumPars="<<nSpectrumPars<<", nTotPars="<<nTotPars<<", spectrumModel="<<spectrumModel<<", xmin/xmax="<<xmin<<"/"<<xmax<<", ymin/ymax="<<ymin<<"/"<<ymax);
	#endif

	//Init response fcn
	TF2* ResponseModelFcn= new TF2("ResponseModel",SpectrumUtils::ResponseModel,xmin,xmax,ymin,ymax,nTotPars);
	ResponseModelFcn->SetNpx(npts);
	ResponseModelFcn->SetNpy(npts);

	//Set spectrum pars
	#ifdef LOGGING_ENABLED
		INFO_LOG("Set spectrum pars (spectrumModel="<<spectrumModel<<") ...");
	#endif

	TF1* spectrumModelFcn= ComputeSpectrumModel(spectrumPars,xmin,xmax);
	if(!spectrumModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute spectrum model (hint: check integral calculation), returning nullptr!");
		#endif
		return nullptr;
	}

	double integral= spectrumPars.GetIntegral(xmin,xmax);
	if(spectrumModel==eFlat){
		integral= (xmax-xmin);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		integral= spectrumModelFcn->Integral(xmin,xmax);
	}	
	else if(spectrumModel==ePol3){
		integral= spectrumModelFcn->Integral(xmin,xmax);
	}	
	else if(spectrumModel==ePol6){
		integral= spectrumModelFcn->Integral(xmin,xmax);
	}	

	if(!std::isfinite(integral) || std::isnan(integral)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Spectrum integral is inf or nan, returning nullptr!");
		#endif
		return nullptr;
	}


	int par_counter= 0;
	ResponseModelFcn->SetParameter(par_counter,spectrumModel);
	par_counter++;
	for(int i=0;i<spectrumParList.GetNPars();i++){
		double parValue= spectrumParList.GetPar(i)->GetValue();
		if(i==0) parValue= 1./integral;
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Spectrum par no. "<<i<<"="<<parValue);
		#endif
		ResponseModelFcn->SetParameter(par_counter,parValue);
		par_counter++;
	}

	//Set bias pars
	#ifdef LOGGING_ENABLED
		INFO_LOG("Set bias pars (biasModel="<<biasModel<<") ...");
	#endif
	ResponseModelFcn->SetParameter(par_counter,biasModel);
	par_counter++;
	for(size_t i=0;i<biasParList.size();i++){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Bias par no. "<<i<<"="<<biasParList[i]);
		#endif
		ResponseModelFcn->SetParameter(par_counter,biasParList[i]);
		par_counter++;
	}

	//Set reso pars
	#ifdef LOGGING_ENABLED
		INFO_LOG("Set reso pars (resoModel="<<resoModel<<") ...");
	#endif
	ResponseModelFcn->SetParameter(par_counter,resoModel);
	par_counter++;
	for(size_t i=0;i<resoParList.size();i++){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reso par no. "<<i<<"="<<resoParList[i]);
		#endif
		ResponseModelFcn->SetParameter(par_counter,resoParList[i]);
		par_counter++;
	}

	
	//Set efficiency pars
	#ifdef LOGGING_ENABLED
		INFO_LOG("Set efficiency pars (effModel="<<effModel<<") ...");
	#endif
	ResponseModelFcn->SetParameter(par_counter,effModel);
	par_counter++;
	for(size_t i=0;i<effParList.size();i++){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Efficiency par no. "<<i<<"="<<effParList[i]);
		#endif
		ResponseModelFcn->SetParameter(par_counter,effParList[i]);
		par_counter++;
	}


	//Check given pars
	for(int i=0;i<ResponseModelFcn->GetNpar();i++){
		double parValue= ResponseModelFcn->GetParameter(i);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Response par["<<i<<"]="<<parValue);
		#endif
	}

	if(spectrumModelFcn) spectrumModelFcn->Delete();

	return ResponseModelFcn;

}//close ComputeResponseModel()


TH2D* SpectrumUtils::ComputeParametricResponse(SpectrumPars& spectrumPars,BiasPars& biasPars,ResoPars& resoPars,EfficiencyPars& effPars,std::vector<double>& TrueBins, std::vector<double>& RecBins,int npts)
{
	//## Check bins
	int NTrueBins= (int)TrueBins.size()-1;
	int NRecBins= (int)RecBins.size()-1;
	if(NTrueBins<=0 || NRecBins<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of bins specified!");
		#endif
		return 0;
	}
	double BinEdge_Rec[NRecBins+1];
	double BinEdge_True[NTrueBins+1];

	for(int s=0;s<NRecBins;s++) {
		BinEdge_Rec[s]= RecBins[s];
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("BinEdge_Rec["<<s<<"]="<<BinEdge_Rec[s]);
		#endif
	}
	BinEdge_Rec[NRecBins]= RecBins[NRecBins];		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("BinEdge_Rec["<<NRecBins<<"]="<<BinEdge_Rec[NRecBins]);
	#endif
	
	for(int s=0;s<NTrueBins;s++) {	
		BinEdge_True[s]= TrueBins[s];
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("BinEdge_True["<<s<<"]="<<BinEdge_True[s]);
		#endif
	}
	BinEdge_True[NTrueBins]= TrueBins[NTrueBins];		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("BinEdge_True["<<NTrueBins<<"]="<<BinEdge_True[NTrueBins]);
	#endif
	
	double LgSmin_true= BinEdge_True[0];
	double LgSmax_true= BinEdge_True[NTrueBins];
	double LgSmin_rec= BinEdge_Rec[0];
	double LgSmax_rec= BinEdge_Rec[NRecBins];


	//## Init spectrum model fcn
	TF1* SpectrumModelFcn= ComputeSpectrumModel(spectrumPars,LgSmin_true,LgSmax_true);
	if(!SpectrumModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the spectrum model!");
		#endif
		return 0;
	}

	//## Init response model
	TF2* ResponseModelFcn= ComputeResponseModel(spectrumPars,biasPars,resoPars,effPars,LgSmin_true,LgSmax_true,LgSmin_rec,LgSmax_rec,npts);
	if(!ResponseModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the response model!");
		#endif
		if(SpectrumModelFcn) SpectrumModelFcn->Delete();
		return 0;
	}

	//## Init response matrix
	TH2D* ResponseMat= new TH2D("","",NTrueBins,BinEdge_True,NRecBins,BinEdge_Rec);
	ResponseMat->Sumw2();


	//## Fill matrix
	for(int i=0;i<ResponseMat->GetNbinsX();i++){
		double lgS_true= ResponseMat->GetXaxis()->GetBinCenter(i+1);
		double binWidth_true= ResponseMat->GetXaxis()->GetBinWidth(i+1);	
		double lgSMin_true= ResponseMat->GetXaxis()->GetBinLowEdge(i+1);
		double lgSMax_true= lgSMin_true + binWidth_true;
		double ProbNorm= SpectrumModelFcn->Integral(lgSMin_true,lgSMax_true);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("binWidth_true="<<binWidth_true<<", lgSMin_true="<<lgSMin_true<<", lgSMax_true="<<lgSMax_true<<", ProbNorm="<<ProbNorm);
		#endif
		
		for(int j=0;j<ResponseMat->GetNbinsY();j++){
			double lgS_rec= ResponseMat->GetYaxis()->GetBinCenter(j+1);
			double binWidth_rec= ResponseMat->GetYaxis()->GetBinWidth(j+1);	
			double lgSMin_rec= ResponseMat->GetYaxis()->GetBinLowEdge(j+1);
			double lgSMax_rec= lgSMin_rec + binWidth_rec;

			#ifdef LOGGING_ENABLED
				DEBUG_LOG("binWidth_rec="<<binWidth_rec<<", lgSMin_rec="<<lgSMin_rec<<", lgSMax_rec="<<lgSMax_rec);
			#endif

			double Rji_noNorm= ResponseModelFcn->Integral(lgSMin_true,lgSMax_true,lgSMin_rec,lgSMax_rec);
			//double Rji_noNorm= ResponseModelFcn->Eval(lgS_true,lgS_rec);
			double Rji= Rji_noNorm/ProbNorm;

			#ifdef LOGGING_ENABLED
				INFO_LOG("lgS_true="<<lgS_true<<" ["<<lgSMin_true<<","<<lgSMax_true<<"], lgS_rec="<<lgS_rec<<" ["<<lgSMin_rec<<","<<lgSMax_rec<<"], ProbNorm="<<ProbNorm<<", Rji_noNorm="<<Rji_noNorm<<" Rji="<<Rji);
			#endif
			
			ResponseMat->SetBinContent(i+1,j+1,Rji);
			ResponseMat->SetBinError(i+1,j+1,0.);
			
		}//end loop rec bins
	}//end loop true bins

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting spectrum model...");
	#endif
	if(SpectrumModelFcn) SpectrumModelFcn->Delete();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting response model...");
	#endif
	if(ResponseModelFcn) ResponseModelFcn->Delete();	
	
	return ResponseMat;

}//close ComputeParametricResponse()





TH2D* SpectrumUtils::BuildResponseMatrix(SpectrumPars& spectrumPars,BiasPars& biasPars,ResoPars& resoPars,EfficiencyPars& effPars,std::vector<double>& TrueBins, std::vector<double>& RecBins,int npts)
{
	//## Check bins
	int NTrueBins= (int)TrueBins.size()-1;
	int NRecBins= (int)RecBins.size()-1;
	if(NTrueBins<=0 || NRecBins<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of bins specified!");
		#endif
		return 0;
	}
	double BinEdge_Rec[NRecBins+1];
	double BinEdge_True[NTrueBins+1];

	for(int s=0;s<NRecBins;s++) {
		BinEdge_Rec[s]= RecBins[s];
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("BinEdge_Rec["<<s<<"]="<<BinEdge_Rec[s]);
		#endif
	}
	BinEdge_Rec[NRecBins]= RecBins[NRecBins];		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("BinEdge_Rec["<<NRecBins<<"]="<<BinEdge_Rec[NRecBins]);
	#endif
	
	for(int s=0;s<NTrueBins;s++) {	
		BinEdge_True[s]= TrueBins[s];
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("BinEdge_True["<<s<<"]="<<BinEdge_True[s]);
		#endif
	}
	BinEdge_True[NTrueBins]= TrueBins[NTrueBins];		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("BinEdge_True["<<NTrueBins<<"]="<<BinEdge_True[NTrueBins]);
	#endif
	
	double LgSmin_true= BinEdge_True[0];
	double LgSmax_true= BinEdge_True[NTrueBins];
	double LgSmin_rec= BinEdge_Rec[0];
	double LgSmax_rec= BinEdge_Rec[NRecBins];


	//## Init spectrum model fcn
	TF1* SpectrumModelFcn= ComputeSpectrumModel(spectrumPars,LgSmin_true,LgSmax_true,npts);
	if(!SpectrumModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the spectrum model!");
		#endif
		return 0;
	}

	//## Init efficiency model fcn
	TF1* EfficiencyModelFcn= ComputeEfficiencyModel(effPars,LgSmin_true,LgSmax_true,npts);
	if(!EfficiencyModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the efficiency model!");
		#endif
		return 0;
	}

	//## Init bias model fcn
	TF1* BiasModelFcn= ComputeBiasModel(biasPars,LgSmin_true,LgSmax_true,npts);
	if(!BiasModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the bias model!");
		#endif
		return 0;
	}

	//## Init reso model fcn
	TF1* ResoModelFcn= ComputeResoModel(resoPars,LgSmin_true,LgSmax_true,npts);
	if(!ResoModelFcn){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute the reso model!");
		#endif
		return 0;
	}

	//## Init response model 1D fcn
	TF1* ResponseModelFcn= new TF1("ResponseModel",SpectrumUtils::ResponseModel1D,LgSmin_true,LgSmax_true,4);
	ResponseModelFcn->SetNpx(npts);
	

	//## Init response matrix
	TH2D* ResponseMat= new TH2D("","",NTrueBins,BinEdge_True,NRecBins,BinEdge_Rec);
	ResponseMat->Sumw2();

	
	//## Fill matrix
	for(int i=0;i<ResponseMat->GetNbinsX();i++){
		double lgS_true= ResponseMat->GetXaxis()->GetBinCenter(i+1);
		double binWidth_true= ResponseMat->GetXaxis()->GetBinWidth(i+1);	
		double lgSMin_true= ResponseMat->GetXaxis()->GetBinLowEdge(i+1);
		double lgSMax_true= lgSMin_true + binWidth_true;
		double N= SpectrumModelFcn->Integral(lgSMin_true,lgSMax_true);

		double eff= EfficiencyModelFcn->Eval(lgS_true);
		double bias= BiasModelFcn->Eval(lgS_true);
		double reso= ResoModelFcn->Eval(lgS_true);
		//ResponseModelFcn->SetParameters(N,bias,reso);
		ResponseModelFcn->SetParameters(1,lgS_true,bias,reso);

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("binWidth_true="<<binWidth_true<<", lgSMin_true="<<lgSMin_true<<", lgSMax_true="<<lgSMax_true<<", bias="<<bias<<", reso="<<reso<<", eff="<<eff);
		#endif
		
		for(int j=0;j<ResponseMat->GetNbinsY();j++){
			double lgS_rec= ResponseMat->GetYaxis()->GetBinCenter(j+1);
			double binWidth_rec= ResponseMat->GetYaxis()->GetBinWidth(j+1);	
			double lgSMin_rec= ResponseMat->GetYaxis()->GetBinLowEdge(j+1);
			double lgSMax_rec= lgSMin_rec + binWidth_rec;

			#ifdef LOGGING_ENABLED
				DEBUG_LOG("binWidth_rec="<<binWidth_rec<<", lgSMin_rec="<<lgSMin_rec<<", lgSMax_rec="<<lgSMax_rec);
			#endif
			
			//double response= ResponseModelFcn->Eval(lgS_rec);
			double response= ResponseModelFcn->Integral(lgSMin_rec,lgSMax_rec);
			double Rji= eff*response;
			
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("lgS_true="<<lgS_true<<" ["<<lgSMin_true<<","<<lgSMax_true<<"], lgS_rec="<<lgS_rec<<" ["<<lgSMin_rec<<","<<lgSMax_rec<<"], N="<<N<<", eff="<<eff<<", bias="<<bias<<", reso="<<reso<<", Rji="<<Rji);
			#endif
			
			ResponseMat->SetBinContent(i+1,j+1,Rji);
			ResponseMat->SetBinError(i+1,j+1,0.);
			
		}//end loop rec bins
	}//end loop true bins
	

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting spectrum model...");
	#endif
	if(SpectrumModelFcn) SpectrumModelFcn->Delete();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting efficiency model...");
	#endif
	if(EfficiencyModelFcn) EfficiencyModelFcn->Delete();	

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting bias model...");
	#endif
	if(BiasModelFcn) BiasModelFcn->Delete();	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting reso model...");
	#endif
	if(ResoModelFcn) ResoModelFcn->Delete();	

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting response model...");
	#endif
	if(ResponseModelFcn) ResponseModelFcn->Delete();	
	

	return ResponseMat;

}//close BuildResponseMatrix()

double SpectrumUtils::PowerLawSpectrum(double* x, double* par)
{
	double lgS= x[0];
	double S= pow(10,lgS);
	double norm= par[0];
	double gamma= par[1]; 
	double spectrum= norm*log(10)*pow(S,-gamma+1);
	
	return spectrum;

}//close PowerLawSpectrum()

double SpectrumUtils::PowerLawWithCutoffSpectrum(double* x, double* par)
{
	double lgS= x[0];	
	double S= pow(10,lgS);
	double norm= par[0];
	double gamma= par[1];
	double lgS0= par[2];
	double Wc= par[3];

	double arg= (lgS-lgS0)/Wc;
	double cutoff= 2./(1.+TMath::Exp(arg));

	double spectrum_before= norm*log(10)*pow(S,-gamma+1);
	double spectrum_after= norm*log(10)*pow(S,-gamma+1)*cutoff;
	double spectrum= 0.;	

	if(lgS<lgS0){
		spectrum= spectrum_before;
	}
	else{
		spectrum= spectrum_after;
	}
	
	double fval= spectrum;

	return fval;

}//close PowerLawWithCutoffSpectrum()


double SpectrumUtils::BrokenPowerLawSpectrum(double* x, double* par)
{
	double lgS= x[0];	
	
	double norm= par[0];
	double gamma1= par[1];
	double gamma2= par[2];
	double gamma3= par[3];
	double LgSb1= par[4];
	double LgSb2= par[5];
	
	double norm1= 1;
	double norm2= norm1* pow(pow(10,LgSb1),gamma2-gamma1);
	double norm3= norm2* pow(pow(10,LgSb2),gamma3-gamma2);
	
	double spectrum1Par[2]= {norm1,gamma1};
	double spectrum2Par[2]= {norm2,gamma2};
	double spectrum3Par[2]= {norm3,gamma3};
		
	double spectrum1= PowerLawSpectrum(&lgS,spectrum1Par);
	double spectrum2= PowerLawSpectrum(&lgS,spectrum2Par);
	double spectrum3= PowerLawSpectrum(&lgS,spectrum3Par);
	
	double spectrum= 0.;
	if(lgS<LgSb1){
		spectrum= spectrum1;
	}
	if(lgS>=LgSb1 && lgS<LgSb2){
		spectrum= spectrum2;
	}
	if(lgS>=LgSb2){
		spectrum= spectrum3;
	}

	return norm*spectrum;

}//close BrokenPowerLawSpectrum()


double SpectrumUtils::TwoBrokenPowerLawSpectrum(double* x, double* par)
{
	double lgS= x[0];	
	
	double norm= par[0];
	double gamma1= par[1];
	double gamma2= par[2];
	double lgS0= par[3];
	
	double norm1= 1;
	double norm2= norm1* pow(pow(10,lgS0),gamma2-gamma1);

	double spectrum1Par[2]= {norm1,gamma1};
	double spectrum2Par[2]= {norm2,gamma2};
		
	double spectrum1= PowerLawSpectrum(&lgS,spectrum1Par);
	double spectrum2= PowerLawSpectrum(&lgS,spectrum2Par);
	
	double spectrum= 0.;
	if(lgS<lgS0){
		spectrum= spectrum1;
	}
	else{
		spectrum= spectrum2;
	}
	
	return norm*spectrum;

}//close TwoBrokenPowerLawSpectrum()

double SpectrumUtils::SmoothCutoffPowerLawSpectrum(double* x, double* par)
{
	double lgS= x[0];
	double S= pow(10,lgS);

	double norm= pow(10,par[0]);
	double gamma1= par[1];
	double gamma2= par[2];
	double LgSb1= par[3];
	double LgSb2= par[4];
	double Wc= par[5];
	
	double cutoffArg= (lgS-LgSb2)/Wc;
	double Cutoff= 1./(1.+TMath::Exp(cutoffArg));

	double Sb1= pow(10,LgSb1);
	//double Sb2= pow(10,LgSb2);
	
	double NormBeforeBreak= norm* pow(Sb1,gamma1-gamma2)*Cutoff;
	double SpectrumAfterBreak= norm* log(10)*pow(S,-gamma2+1)*Cutoff;
	double SpectrumBeforeBreak= NormBeforeBreak* log(10)*pow(S,-gamma1+1);
	
	double Spectrum= 0.;	

	if(lgS<LgSb1){
		Spectrum= SpectrumBeforeBreak;
	}
	else{
		Spectrum= SpectrumAfterBreak;
	}
	
	double fval= Spectrum;

	return fval;

}//close SmoothCutoffPowerLawSpectrum()


double SpectrumUtils::Pol3Spectrum(double* x, double* par)
{
	double lgS= x[0];
	double p0= par[0];
	double p1= par[1];
	double p2= par[2];
	double p3= par[3];
	double fval= p0 + p1*lgS + p2*lgS*lgS + p3*lgS*lgS*lgS;
	//double counts= p0 + p1*(lgS+3) + p2*pow(lgS+3,2) + p3*pow(lgS+3,3);
	//double fval= pow(10,counts);

	return fval;

}//close Pol3Spectrum()


double SpectrumUtils::Pol6Spectrum(double* x, double* par)
{
	double lgS= x[0];
	double p0= par[0];
	double p1= par[1];
	double p2= par[2];
	double p3= par[3];
	double p4= par[4];
	double p5= par[5];
	double p6= par[6];
	double fval= p0 + p1*lgS + p2*pow(lgS,2) + p3*pow(lgS,3) + p4*pow(lgS,4) + p5*pow(lgS,5) + p6*pow(lgS,6);
	//double counts= p0 + p1*(lgS+3) + p2*pow(lgS+3,2) + p3*pow(lgS+3,3) + p4*pow(lgS+3,4) + p5*pow(lgS+3,5) + p6*pow(lgS+3,6);
	//double fval= pow(10,counts);

	return fval;

}//close Pol6Spectrum()

double SpectrumUtils::GetPowerLawIntegral(double gamma,double lgSMin, double lgSMax)
{
	double SpectrumInt= 0;
	if(gamma==1){
		SpectrumInt= (lgSMax-lgSMin);
	}
	else{
		double SMin= pow(10,lgSMin);
		double SMax= pow(10,lgSMax);
		SpectrumInt= (pow(SMax,-gamma+1)-pow(SMin,-gamma+1))/(-gamma+1);
	}

	return SpectrumInt;

}//close GetPowerLawIntegral()


double SpectrumUtils::GetBrokenPowerLawIntegral(double Gamma1,double Gamma2,double Gamma3,double Break,double Cutoff,double lgSMin, double lgSMax) 
{
	double norm1= 1;
	double norm2= norm1* pow(pow(10,Break),Gamma2-Gamma1);
	double norm3= norm2* pow(pow(10,Cutoff),Gamma3-Gamma2);
	double SpectrumInt1= norm1*SpectrumUtils::GetPowerLawIntegral(Gamma1,lgSMin,Break);
	double SpectrumInt2= norm2*SpectrumUtils::GetPowerLawIntegral(Gamma2,Break,Cutoff);
	double SpectrumInt3= norm3*SpectrumUtils::GetPowerLawIntegral(Gamma3,Cutoff,lgSMax);
	double SpectrumInt= SpectrumInt1+SpectrumInt2+SpectrumInt3;

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("gamma1="<<Gamma1<<", gamma2="<<Gamma2<<", gamma3="<<Gamma3<<", lgS_break="<<Break<<", lgS_cutoff="<<Cutoff<<", min/max="<<lgSMin<<"/"<<lgSMax<<", norm1="<<norm1<<", norm2="<<norm2<<", SpectrumInt1="<<SpectrumInt1<<", SpectrumInt2="<<SpectrumInt2<<", SpectrumInt3="<<SpectrumInt3<<", SpectrumInt="<<SpectrumInt);
	#endif

	return SpectrumInt;

}//close GetBrokenPowerLawIntegral()

double SpectrumUtils::GetTwoBrokenPowerLawIntegral(double Gamma1,double Gamma2,double Break,double lgSMin, double lgSMax) 
{
	double norm1= 1;
	double norm2= norm1* pow(pow(10,Break),Gamma2-Gamma1);
	double SpectrumInt1= norm1*SpectrumUtils::GetPowerLawIntegral(Gamma1,lgSMin,Break);
	double SpectrumInt2= norm2*SpectrumUtils::GetPowerLawIntegral(Gamma2,Break,lgSMax);
	double SpectrumInt= SpectrumInt1 + SpectrumInt2;

	return SpectrumInt;

}//close GetTwoBrokenPowerLawIntegral()



double SpectrumUtils::ExpResolutionModel(double* x, double* par)
{
	double lgS= x[0];
	double norm= par[0];
	double beta= par[1];
	double fval= norm*exp(-beta*lgS);
	return fval;

}//close ExpResolutionModel()

double SpectrumUtils::ExpStepResolutionModel(double* x, double* par)
{
	double lgS= x[0];
	double norm= par[0];
	double beta= par[1];
	double lgS_step= par[2];
	double fval= norm*exp(-beta*lgS);
	double fval_step= norm*exp(-beta*lgS_step);
	if(lgS<=lgS_step) fval= fval_step; 

	return fval;

}//close ExpStepResolutionModel()

double SpectrumUtils::ExpBiasModel(double* x, double* par)
{
	double lgS= x[0];
	double norm= par[0];
	double beta= par[1];
	double offset= par[2];
	double fval= offset + norm*exp(-beta*lgS);
	return fval;

}//close ExpBiasModel()

double SpectrumUtils::StepExpBiasModel(double* x, double* par)
{
	double lgS= x[0];
	double norm= par[0];
	double beta= par[1];
	double offset= par[2];
	double lgS_step= par[3];
	double fval= offset + norm*exp(-beta*lgS);
	double fval_step= offset + norm*exp(-beta*lgS_step);
	if(lgS<=lgS_step) fval= fval_step; 

	return fval;

}//close StepExpBiasModel()

double SpectrumUtils::EfficiencyModel(double* x, double* par)
{
	double lgS= x[0];
	double norm= par[0];
	double lgS0= par[1];
	double k= par[2];
	double fval= norm*0.5*(1. + TMath::Erf((lgS-lgS0)/k));
	return fval;

}//close EfficiencyModel()


double SpectrumUtils::ResponseModel1D(double* x, double* par)
{
	double lgS_rec= x[0];
	double norm= par[0];
	double lgS_true= par[1];
	double bias= par[2];
	double reso= par[3];
	
	//double arg = (lgS + bias)/reso;	
	double arg= (lgS_rec - (lgS_true + bias) )/reso;
	double gaussNorm= 1./(reso*sqrt(2.*TMath::Pi()));
	double fval = norm*gaussNorm*TMath::Exp(-0.5*arg*arg);
	return fval;

}//close ResponseModel1D()

double SpectrumUtils::ResponseModel(double* x, double* par)
{
	double lgS_true= x[0];
	double lgS_rec= x[1];
	double S_true= pow(10,lgS_true);
	double S_rec= pow(10,lgS_rec);
	int par_counter= 0;

	//## Set spectrum pars
	int spectrumModel= par[par_counter++];
	int nSpectrumPars= 0;
	
	if(spectrumModel==ePowerLaw){
		nSpectrumPars= PowerLawPars::GetParNumber();
	}
	else if(spectrumModel==eBrokenPowerLaws){
		nSpectrumPars= BrokenPowerLawsPars::GetParNumber();
	}
	else if(spectrumModel==eTwoBrokenPowerLaws){
		nSpectrumPars= TwoBrokenPowerLawsPars::GetParNumber();
	}
	else if(spectrumModel==eFlat){
		nSpectrumPars= FlatSpectrumPars::GetParNumber();
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		nSpectrumPars= SmoothCutoffPowerLaws::GetParNumber();
	}
	else if(spectrumModel==ePowerLawWithCutoff){
		nSpectrumPars= PowerLawWithCutoff::GetParNumber();
	}
	else if(spectrumModel==ePol3){
		nSpectrumPars= Pol3SpectrumPars::GetParNumber();
	}
	else if(spectrumModel==ePol6){
		nSpectrumPars= Pol6SpectrumPars::GetParNumber();
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unknown spectrum model given ("<<spectrumModel<<"), returning zero!");
		#endif
		return 0;
	}

	double spectrumPars[nSpectrumPars];
	for(int i=0;i<nSpectrumPars;i++){
		spectrumPars[i]= par[par_counter];
		par_counter++;
	}

	double spectrum= 0;
	
	if(spectrumModel==ePowerLaw){
		spectrum= PowerLawSpectrum(&lgS_true,spectrumPars);
	}
	else if(spectrumModel==eBrokenPowerLaws){
		spectrum= BrokenPowerLawSpectrum(&lgS_true,spectrumPars);
	}		
	else if(spectrumModel==eTwoBrokenPowerLaws){
		spectrum= TwoBrokenPowerLawSpectrum(&lgS_true,spectrumPars);
	}	
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		spectrum= SmoothCutoffPowerLawSpectrum(&lgS_true,spectrumPars);
	}	
	else if(spectrumModel==ePowerLawWithCutoff){
		spectrum= PowerLawWithCutoffSpectrum(&lgS_true,spectrumPars);
	}	
	else if(spectrumModel==eFlat){
		spectrum= 1;
	}
	else if(spectrumModel==ePol3){
		spectrum= Pol3Spectrum(&lgS_true,spectrumPars);
	}
	else if(spectrumModel==ePol6){
		spectrum= Pol6Spectrum(&lgS_true,spectrumPars);
	}

	//## Set bias pars
	int biasModel= par[par_counter++];
	int nBiasPars= 0;
	
	if(biasModel==eCONST_BIAS){
		nBiasPars= ConstBiasPars::GetParNumber();
	}	
	else if(biasModel==eEXP_BIAS){
		nBiasPars= ExpBiasPars::GetParNumber();
	}
	else if(biasModel==eEXP_STEP_BIAS){
		nBiasPars= ExpStepBiasPars::GetParNumber();
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unknown bias model given ("<<biasModel<<"), returning zero!");
		#endif
		return 0;
	}

	double biasPars[nBiasPars];
	for(int i=0;i<nBiasPars;i++){
		biasPars[i]= par[par_counter];	
		par_counter++;
	}

	double bias= 0;
	if(biasModel==eCONST_BIAS){
		bias= biasPars[0];
	}	
	else if(biasModel==eEXP_BIAS){
		bias= ExpBiasModel(&lgS_true,biasPars);
	}	
	else if(biasModel==eEXP_STEP_BIAS){
		bias= StepExpBiasModel(&lgS_true,biasPars);
	}

	//## Set reso pars
	int resoModel= par[par_counter++];
	int nResoPars= 0;
	if(resoModel==eCONST_RESO){
		nResoPars= ConstResoPars::GetParNumber();
	}	
	else if(resoModel==eEXP_RESO){
		nResoPars= ExpResoPars::GetParNumber();
	}
	else if(resoModel==eEXP_STEP_RESO){
		nResoPars= ExpStepResoPars::GetParNumber();
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unknown reso model given ("<<resoModel<<"), returning zero!");
		#endif
		return 0;
	}

	double resoPars[nResoPars];
	for(int i=0;i<nResoPars;i++){
		resoPars[i]= par[par_counter];
		//#ifdef LOGGING_ENABLED
		//	DEBUG_LOG("resoPars["<<i<<"]="<<resoPars[i]);
		//#endif
		par_counter++;
	}

	double reso= 0;
	if(resoModel==eCONST_RESO){
		reso= resoPars[0];
	}	
	else if(resoModel==eEXP_RESO){
		reso= ExpResolutionModel(&lgS_true,resoPars);
	}
	else if(resoModel==eEXP_STEP_RESO){
		reso= ExpStepResolutionModel(&lgS_true,resoPars);
	}

	//## Set efficiency pars
	int effModel= par[par_counter++];
	int nEffPars= 0;
	if(effModel==eCONST_EFF){
		nEffPars= ConstEfficiencyPars::GetParNumber();
	}	
	else if(effModel==eSIGMOID_EFF){
		nEffPars= SigmoidEfficiencyPars::GetParNumber();
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unknown efficiency model given ("<<effModel<<"), returning zero!");
		#endif
		return 0;
	}

	double effPars[nEffPars];
	for(int i=0;i<nEffPars;i++){
		effPars[i]= par[par_counter];
		par_counter++;
	}

	double efficiency= 1;
	
	if(effModel==eCONST_EFF){
		efficiency= effPars[0];
	}	
	else if(effModel==eSIGMOID_EFF){
		efficiency= EfficiencyModel(&lgS_rec,effPars);
	}
	
	//## Set response function
  double arg = (lgS_rec - (lgS_true + bias) )/reso;
	double gaussNorm= 1./(reso*sqrt(2.*TMath::Pi()));
	double gaussResponse= gaussNorm*TMath::Exp(-0.5*arg*arg);
  double response = spectrum * gaussResponse * efficiency;
	 
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("lgS_true="<<lgS_true<<", lgS_rec="<<lgS_rec<<", reso="<<reso<<", bias="<<bias<<", spectrum="<<spectrum<<", eff="<<efficiency<<", arg="<<arg<<", gaussNorm="<<gaussNorm<<", gaussResponse="<<gaussResponse<<", response="<<response);
	#endif
	

	return response;

}//close ResponseModel()


//==========================================
//      SOURCE COUNTS RESO BIAS MODEL
//==========================================
double SpectrumUtils::PhiMaxModel(double* x,double* pars)
{
	double S_Jy= x[0];
	double S= S_Jy*1.e+3;//mJy
	double bmaj= pars[0];
	double bmin= pars[1];
	double Sthr= pars[2];//in mJy
	double sourceAngSize= 0;
	if(S/Sthr>1) sourceAngSize= sqrt(bmaj*bmin)*sqrt(S/Sthr-1);
	
	return sourceAngSize;

}//close PhiMaxPDF()


double SpectrumUtils::PhiMinModel(double* x,double* pars)
{
	double Sp_Jy= x[0];
	double Sp= Sp_Jy*1.e+3;//mJy
	double bmaj= pars[0];
	double bmin= pars[1];
	double rms= pars[2];//mJy
	double p0= pars[3];
	double p1= pars[4];
	double beta= pars[5];
	
	double arg= p0+p1/pow(Sp/rms,beta);
	//double arg= p0+p1/pow(Sp/rms,beta)-1;
	double sourceAngSize= 0;
	if(arg>=0) sourceAngSize= sqrt(bmaj*bmin)*sqrt(arg);
	
	return sourceAngSize;

}//close PhiMinPDF()

double SpectrumUtils::PhiMedianModel(double* x,double* pars)
{
	double S_Jy= x[0];
	double S= S_Jy*1.e+3;//mJy
	double nu= pars[0];
	double alpha= pars[1];
	double beta= pars[2];
	double S_break= pars[3];//in mJy
	double scaleFactor= pars[4]; 

	double nu_ref= 1.4;//GHz
	double S_ref= S*pow(nu_ref/nu,alpha);
	double phi0= 2.0;//in arcsec
	double phi= 0;//in arcsec
	
	//if(S<S_break) phi= phi0;
	if(S_ref<S_break) phi= phi0;
	else phi= phi0*pow(S_ref,beta);
	

	phi*= scaleFactor;
	
	//cout<<"S="<<S<<", S_ref="<<S_ref<<", alpha="<<alpha<<", nu="<<nu<<", beta="<<beta<<", scaleFactor="<<scaleFactor<<endl;

	return phi;

}//close PhiMedianPDF()


double SpectrumUtils::SourceCountsResoBiasModel(double* x,double* pars)
{
	double S_Jy= x[0];
	double S= S_Jy*1.e+3;//mJy
	double bmaj= pars[0];
	double bmin= pars[1];
	double rms= pars[2];
	double Sthr= pars[3];
	double p0= pars[4];
	double p1= pars[5];
	double beta= pars[6];
	double nu= pars[7];
	double alpha= pars[8];
	double beta_phimedian= pars[9];
	double S_break= pars[10];
	double scaleFactor= pars[11];
	
	//Compute phi min/max
	double pars_phiMin[]= {bmaj,bmin,rms,p0,p1,beta};
	double phiMin= SpectrumUtils::PhiMinModel(x,pars_phiMin);

	double pars_phiMax[]= {bmaj,bmin,Sthr};
	double phiMax= SpectrumUtils::PhiMaxModel(x,pars_phiMax);
	
	//Compute phi limit
	double phiLimit= std::max(phiMin,phiMax);
	
	//Compute phi median
	double pars_phiMedian[]= {nu,alpha,beta_phimedian,S_break,scaleFactor};
	double phiMedian= SpectrumUtils::PhiMedianModel(x,pars_phiMedian);

	//Compute h fcn
	double h= exp(-log(2.)*pow(phiLimit/phiMedian,0.62));

	//cout<<"S="<<S<<", phi min/max/prime="<<phiMin<<"/"<<phiMax<<"/"<<phiMaxPrime<<", phiLimit="<<phiLimit<<", phiMedian="<<phiMedian<<", phiLimit/phiMedian="<<phiLimit/phiMedian<<", h="<<h<<endl; 
	
	return h;

}//close SourceCountsResoBiasModel()

double SpectrumUtils::SourceCountsResoBiasCorrFactor(double* x,double* pars)
{
	double lgS= x[0];
	double S= pow(10,lgS);
	double h= SpectrumUtils::SourceCountsResoBiasModel(&S,pars);
	double corrFactor= 1;
	if(h!=1) corrFactor= 1./(1.-h);

	return corrFactor;

}//close SourceCountsResoBiasCorrFactor()

//==========================================
//      SPECTRUM PARS
//==========================================

double FlatSpectrumPars::GetIntegral(double xmin,double xmax) 
{
	//double Gamma= (pars.GetPar(1))->GetValue();
	double Gamma= 1;
	return SpectrumUtils::GetPowerLawIntegral(Gamma,xmin,xmax);
}

double PowerLawPars::GetIntegral(double xmin,double xmax) 
{
	double Gamma= (pars.GetPar(1))->GetValue();
	return SpectrumUtils::GetPowerLawIntegral(Gamma,xmin,xmax);
}

double BrokenPowerLawsPars::GetIntegral(double xmin,double xmax) 
{
	double Gamma1= pars.GetPar(1)->GetValue();
	double Gamma2= pars.GetPar(2)->GetValue();
	double Gamma3= pars.GetPar(3)->GetValue();
	double Break= pars.GetPar(4)->GetValue();
	double Cutoff= pars.GetPar(5)->GetValue();
	return SpectrumUtils::GetBrokenPowerLawIntegral(Gamma1,Gamma2,Gamma3,Break,Cutoff,xmin,xmax);
}

double TwoBrokenPowerLawsPars::GetIntegral(double xmin,double xmax) 
{
	double Gamma1= pars.GetPar(1)->GetValue();
	double Gamma2= pars.GetPar(2)->GetValue();
	double Break= pars.GetPar(3)->GetValue();
	return SpectrumUtils::GetTwoBrokenPowerLawIntegral(Gamma1,Gamma2,Break,xmin,xmax);
}




