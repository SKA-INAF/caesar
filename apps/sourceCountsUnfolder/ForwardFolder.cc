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
* @file ForwardFolder.cc
* @class ForwardFolder
* @brief Forward folding of a count spectrum
* 
* @author S. Riggi
* @date 03/07/2019
*/

#include <ForwardFolder.h>
#include <SpectrumUtils.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TDecompChol.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>

using namespace std;
using namespace Caesar;


//==========================================
//      STATIC VARS
//==========================================
bool ForwardFolder::fUseFitRange;
double ForwardFolder::fLgEMin_fit;
double ForwardFolder::fLgEMax_fit;
bool ForwardFolder::fUseErrorsInChi2= false;
bool ForwardFolder::fIsLikelihoodFit= false;
SpectrumPars* ForwardFolder::fInitFitPars;

TF1* ForwardFolder::fTrueSpectrumModelFcn;
TH1D* ForwardFolder::fCurrentRecSpectrum;		
TH1D* ForwardFolder::fTrueSpectrum;
TH1D* ForwardFolder::fForwardFoldedSpectrum;
TH1D* ForwardFolder::fCurrentForwardFoldedSpectrum;
TH2D* ForwardFolder::fResponseMatrix;



ForwardFolder::ForwardFolder()
{
	//Initialize class data
	fTrueSpectrumModelFcn= 0;
	fUnfoldedSpectrum= 0;
	fTrueSpectrum= 0;
	fForwardFoldedSpectrum= 0;
	fCurrentForwardFoldedSpectrum= 0;
	fCurrentRecSpectrum= 0;
	fCovarianceMatrix= 0;

}//close costructor


ForwardFolder::~ForwardFolder()
{
	//Clear class data
	Clear();

}//close destructor

int ForwardFolder::Clear()
{
	//Clear data	
	if(fTrueSpectrumModelFcn) fTrueSpectrumModelFcn->Delete();
	if(fUnfoldedSpectrum) fUnfoldedSpectrum->Delete();
	if(fTrueSpectrum) fTrueSpectrum->Delete();
	if(fForwardFoldedSpectrum) fForwardFoldedSpectrum->Delete();
	if(fCurrentForwardFoldedSpectrum) fCurrentForwardFoldedSpectrum->Delete();
	if(fCurrentRecSpectrum) fCurrentRecSpectrum->Delete();
	if(fCovarianceMatrix) fCovarianceMatrix->Delete();

	return 0;

}//close Clear()

int ForwardFolder::Init()
{
	// - Clear existing data
	if(Clear()<0) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to clear data!");
		#endif
		return -1;
	}

	// - Set spectrum & matrix bins
	if(!fRecSpectrum || !fResponseMatrix){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Rec spectrum or response matrix are not available (this should not occur)!");
		#endif
		return -1;
	}

	fNTrueBins= fResponseMatrix->GetNbinsX();
	fNRecBins= fResponseMatrix->GetNbinsY();
	
	fLgEBins_true.assign(fNTrueBins+1,0);
	fLgEBins_rec.assign(fNRecBins+1,0);
	
	double BinEdge_RecSpectrum[fNRecBins+1];
	double BinEdge_TrueSpectrum[fNTrueBins+1];

	cout<<"== TRUE BINS =="<<endl;
	cout<<"BinEdge_TrueSpectrum={";
	for(int s=0;s<fNTrueBins;s++) {
		fLgEBins_true[s]= fResponseMatrix->GetXaxis()->GetBinLowEdge(s+1);
		BinEdge_TrueSpectrum[s]= fLgEBins_true[s];
		cout<<BinEdge_TrueSpectrum[s]<<",";
	}
	fLgEBins_true[fNTrueBins]= fLgEBins_true[fNTrueBins-1]+fResponseMatrix->GetXaxis()->GetBinWidth(fNTrueBins);
	BinEdge_TrueSpectrum[fNTrueBins]= fLgEBins_true[fNTrueBins];
	cout<<BinEdge_TrueSpectrum[fNTrueBins]<<"}"<<endl;
	cout<<endl;

	cout<<"== REC BINS =="<<endl;
	cout<<"BinEdge_RecSpectrum={";
	for(int s=0;s<fNRecBins;s++) {
		fLgEBins_rec[s]= fResponseMatrix->GetYaxis()->GetBinLowEdge(s+1);
		BinEdge_RecSpectrum[s]= fLgEBins_rec[s];
		cout<<BinEdge_RecSpectrum[s]<<",";
	}
	fLgEBins_rec[fNRecBins]= fLgEBins_rec[fNRecBins-1]+fResponseMatrix->GetXaxis()->GetBinWidth(fNRecBins);
	BinEdge_RecSpectrum[fNRecBins]= fLgEBins_rec[fNRecBins];
	cout<<BinEdge_RecSpectrum[fNRecBins]<<"}"<<endl;
	cout<<endl;

	// - Ensure rec spectrum and response matrix have no bin offsets
	bool hasSameBinLimits= SpectrumUtils::CheckSpectrumBinnings(fRecSpectrum->GetXaxis(),fResponseMatrix->GetXaxis());
	if(!hasSameBinLimits){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Response matrix and rec spectrum have different bin limits!");
		#endif
		return -1;
	}
	

	// - Initialize spectrum model
	int npts= 1000;
	if(fTrueSpectrumModelFcn) fTrueSpectrumModelFcn->Delete();
	int spectrumModel= fInitFitPars->GetModel();
	int nSpectrumPars= fInitFitPars->GetNPars();
	if(spectrumModel==ePowerLaw){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::PowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eBrokenPowerLaws){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::BrokenPowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eTwoBrokenPowerLaws){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::TwoBrokenPowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::SmoothCutoffPowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==ePowerLawWithCutoff){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::PowerLawWithCutoffSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eFlat){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn","[0]",fLgEMin_true,fLgEMax_true);
	}
	else if(spectrumModel==ePol3){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::Pol3Spectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==ePol6){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",SpectrumUtils::Pol6Spectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid spectrum model selected!");
		#endif
		return -1;
	}

	for(int i=0;i<nSpectrumPars;i++){
		double parValue= fInitFitPars->GetPar(i)->GetValue();
		fTrueSpectrumModelFcn->SetParameter(i,parValue);	
	}

	fTrueSpectrumModelFcn->SetNpx(npts);
	
	// - Initialize spectrum histos
	if(fUnfoldedSpectrum) fUnfoldedSpectrum->Delete();
	fUnfoldedSpectrum= new TH1D("UnfoldedSpectrum","UnfoldedSpectrum",fNTrueBins,BinEdge_TrueSpectrum);
	fUnfoldedSpectrum->Sumw2();
	fUnfoldedSpectrum->SetMarkerStyle(24);
	fUnfoldedSpectrum->SetMarkerSize(1.1);
	fUnfoldedSpectrum->SetMarkerColor(kBlack);
	fUnfoldedSpectrum->SetLineColor(kBlack);
	
	if(fTrueSpectrum) fTrueSpectrum->Delete();
	fTrueSpectrum= new TH1D("TrueSpectrum","TrueSpectrum",fNTrueBins,BinEdge_TrueSpectrum);
	fTrueSpectrum->Sumw2();
	fTrueSpectrum->SetMarkerStyle(8);
	fTrueSpectrum->SetMarkerSize(1.1);
	fTrueSpectrum->SetMarkerColor(kBlack);
	fTrueSpectrum->SetLineColor(kBlack);

	if(fForwardFoldedSpectrum) fForwardFoldedSpectrum->Delete();
	fForwardFoldedSpectrum= new TH1D("ForwardFoldedSpectrum","ForwardFoldedSpectrum",fNRecBins,BinEdge_RecSpectrum);
	fForwardFoldedSpectrum->Sumw2();
	fForwardFoldedSpectrum->SetMarkerStyle(8);
	fForwardFoldedSpectrum->SetMarkerSize(1.1);
	fForwardFoldedSpectrum->SetMarkerColor(kBlack);
	fForwardFoldedSpectrum->SetLineColor(kBlack);
	
	if(fCurrentForwardFoldedSpectrum) fCurrentForwardFoldedSpectrum->Delete();
	fCurrentForwardFoldedSpectrum= new TH1D("CurrentForwardFoldedSpectrum","CurrentForwardFoldedSpectrum",fNRecBins,BinEdge_RecSpectrum);
	fCurrentForwardFoldedSpectrum->Sumw2();
	fCurrentForwardFoldedSpectrum->SetMarkerStyle(8);
	fCurrentForwardFoldedSpectrum->SetMarkerSize(1.1);
	fCurrentForwardFoldedSpectrum->SetMarkerColor(kBlack);
	fCurrentForwardFoldedSpectrum->SetLineColor(kBlack);
	
	if(fCurrentRecSpectrum) fCurrentRecSpectrum->Delete();
	fCurrentRecSpectrum= new TH1D("CurrentRecSpectrum","CurrentRecSpectrum",fNRecBins,BinEdge_RecSpectrum);
	fCurrentRecSpectrum->Sumw2();
	fCurrentRecSpectrum->SetMarkerStyle(8);
	fCurrentRecSpectrum->SetMarkerSize(1.1);
	fCurrentRecSpectrum->SetMarkerColor(kBlack);
	fCurrentRecSpectrum->SetLineColor(kBlack);

	cout<<"== REC SPECTRUM =="<<endl;
	fCurrentRecSpectrum->Reset();
	for(int i=0;i<fRecSpectrum->GetNbinsX();i++) {
		cout<<"Bin "<<i+1<<": x="<<fRecSpectrum->GetBinCenter(i+1)<<", y="<<fRecSpectrum->GetBinContent(i+1)<<endl;
		fCurrentRecSpectrum->SetBinContent(i+1,fRecSpectrum->GetBinContent(i+1));
		fCurrentRecSpectrum->SetBinError(i+1,fRecSpectrum->GetBinError(i+1));
	}

	cout<<"== REC SPECTRUM CURRENT =="<<endl;
	for(int i=0;i<fCurrentRecSpectrum->GetNbinsX();i++) {
		cout<<"Bin "<<i+1<<": x="<<fCurrentRecSpectrum->GetBinCenter(i+1)<<", y="<<fCurrentRecSpectrum->GetBinContent(i+1)<<endl;
	}

	// - Initialize fit covariance matrix
	if(fCovarianceMatrix) fCovarianceMatrix->Delete();
	fCovarianceMatrix= new TMatrixD(nSpectrumPars,nSpectrumPars);

	// - Init random number generator
	delete gRandom;
	gRandom= new TRandom3;

	return 0;

}//close Init()

int ForwardFolder::RunUnfold(TH1D* recSpectrum,TH2D* responseMatrix,SpectrumPars& initFitPars,bool computeUncertainties,int nRandomSamples,bool useFitRange,double fitMin,double fitMax)
{
	// - Check input spectrum/matrix
	if(!recSpectrum || !responseMatrix){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input spectrum given!");
		#endif
		return -1;
	}
	fRecSpectrum= recSpectrum;
	fResponseMatrix= responseMatrix;
	fInitFitPars= initFitPars.Clone();

	fUseFitRange= useFitRange;
	fLgEMin_fit= fRecSpectrum->GetXaxis()->GetXmin();
	fLgEMax_fit= fRecSpectrum->GetXaxis()->GetXmax();
	if(useFitRange){
		fLgEMin_fit= fitMin;
		fLgEMax_fit= fitMax;
	}

	// - Initialize data
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initialize data ...");
	#endif
	if(Init()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to initialize class data!");
		#endif
		return -1;
	}

	// - Unfold spectrum
	#ifdef LOGGING_ENABLED
		INFO_LOG("Start spectrum unfolding ...");
	#endif
	fUnfoldedSpectrum= UnfoldSpectrum(initFitPars,"FIT");

	// - Copyng final forward folded spectrum 
	SpectrumUtils::CopySpectrumContent(fCurrentForwardFoldedSpectrum,fForwardFoldedSpectrum);

	//- Calculate syst uncertainties
	if(computeUncertainties) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("Calculating systematics on the unfolded spectrum ...");
		#endif
		ComputeUncertainties(initFitPars,nRandomSamples);
	}
	return 0;

}//close Unfold()


TH1D* ForwardFolder::UnfoldSpectrum(SpectrumPars& initFitPars,std::string runMode)
{
	//#######################
	//##    INIT FIT
	//#######################	
	double arglist[10];
  int ierflag = 0;
	double amin,edm,errdef;
  int nvpar,nparx,icstat;
	int nPar= initFitPars.GetNPars();	
	
	TMinuit* gMinuit = new TMinuit(nPar);
  gMinuit->SetPrintLevel(0);
	gMinuit->SetMaxIterations(100000);
	gMinuit->SetFCN(ForwardFolder::MinimizedFcn);
	gMinuit->mnexcm("SET NOW",arglist,0,ierflag);

	arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflag);

	//## Set starting values and step sizes for parameters
	for(int i=0;i<nPar;i++){
		FitPar* thisFitPar= initFitPars.GetPar(i);
		std::string parName= thisFitPar->GetParName();
		double parValue= thisFitPar->GetValue();
		double stepValueRel= thisFitPar->GetStepSize();
		double stepValue= fabs(stepValueRel*parValue);
		bool isLimited= thisFitPar->IsLimited();
		double parMinValue= thisFitPar->GetMinValue();
		double parMaxValue= thisFitPar->GetMaxValue();
		double isFixed= thisFitPar->IsFixed();
		if(isLimited) gMinuit->mnparm(i,parName.c_str(), parValue, stepValue, parMinValue, parMaxValue, ierflag);
		else gMinuit->mnparm(i,parName.c_str(), parValue, stepValue, 0, 0, ierflag);

		if(isFixed) gMinuit->FixParameter(i);
	}

	
	//## SET MINUIT STRATEGY
  // 0 ==> low level minimization but small number of FCN calls
  // 1 ==> intermediate
  // 2 ==> max level but many FCN calls
  arglist[0]= 1;
  gMinuit->mnexcm("SET STR",arglist,1,ierflag);

	arglist[0]= 1;//0.5 likelihood, 1 ChiSquare
	if(fIsLikelihoodFit) arglist[0]= 0.5;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflag);

	if(runMode=="CALL"){ //## Call Chi2 with initial values
		gMinuit->mnexcm("CALL FCN", arglist ,0,ierflag);
	}
	else{ //## Fit
  	arglist[0] = 100000;
  	arglist[1] = 0.1;
  	gMinuit->mnexcm("MINIMIZE", arglist ,2,ierflag);
	}

	//## Print results
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  
	fFitStatus= gMinuit->GetStatus();

	//## Get fitted parameters
	double p,ep;
	double FinalPars[nPar];
	double FinalParsErr[nPar];

	cout<<"*** FIT RESULTS ***"<<endl;
	for(int i=0;i<nPar;i++){
		gMinuit->GetParameter(i,p,ep);	
		FinalPars[i]= p;
		FinalParsErr[i]= ep;
		fFitPar[i]= p;
		fFitParErr[i]= ep;
		cout<<"FinalPar"<<i+1<<"="<<FinalPars[i]<<"  +- "<<FinalParsErr[i]<<endl;
	}
	

	//## Get covariance matrix
	double covMatrix[nPar][nPar];		
  gMinuit->mnemat(&covMatrix[0][0],nPar);
  	
  cout<<"*** FIT COV MATRIX ***"<<endl;  	
	TMatrixD CovarianceMatrix(nPar,nPar);
		
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){
		 	CovarianceMatrix(i,j)= covMatrix[i][j];		
			(*fCovarianceMatrix)(i,j)= covMatrix[i][j];	
		  if(j==nPar-1) cout<<covMatrix[i][j]<<endl;
		  else cout<<covMatrix[i][j]<<"  ";
		}//close for j
	}//close for i
	cout<<endl;

	//## Correction factor for rec flux
	//## Create unfolded histo and fill with results
	TH1D* aUnfoldedSpectrumHisto= (TH1D*)fUnfoldedSpectrum->Clone(); 
	aUnfoldedSpectrumHisto->Reset();
	fTrueSpectrumModelFcn->SetParameters(FinalPars);

	double nTot= fCurrentRecSpectrum->Integral();
	double nTrueTot= fTrueSpectrum->Integral();

	//for(int j=0;j<fNRecBins;j++){
	for(int j=0;j<fCurrentRecSpectrum->GetNbinsX();j++){
		double lgERec= fCurrentRecSpectrum->GetBinCenter(j+1);
		double nRec= fCurrentRecSpectrum->GetBinContent(j+1);
		double nRecErr= fCurrentRecSpectrum->GetBinError(j+1);
		int recBinId= fTrueSpectrum->FindBin(lgERec);
				
		double nTrue= fTrueSpectrum->GetBinContent(recBinId);
		double nTrueErr= fTrueSpectrum->GetBinError(recBinId);
		fTrueSpectrum->SetBinContent(recBinId,nTrue);
		fTrueSpectrum->SetBinError(recBinId,nTrueErr);

		//double nFF= fForwardFoldedSpectrum->GetBinContent(j+1);
		double nFF= fCurrentForwardFoldedSpectrum->GetBinContent(j+1);
		
		double correctionFactor= 0;
		if(nFF!=0) correctionFactor= (double)(nTrue)/(double)(nFF);
		double nUnfold= nRec*correctionFactor;
		double nUnfoldErr= nRecErr*correctionFactor;
		aUnfoldedSpectrumHisto->SetBinContent(recBinId,nUnfold);
		aUnfoldedSpectrumHisto->SetBinError(recBinId,nUnfoldErr);

		#ifdef LOGGING_ENABLED
			INFO_LOG("lgERec="<<lgERec<<", nRec="<<nRec<<", nTrue="<<nTrue<<", nFF="<<nFF<<", corr="<<correctionFactor<<", nUnfold="<<nUnfold<<", Unfold/nRec-1="<<nUnfold/nRec-1);
		#endif

	}//end loop rec bins

	return aUnfoldedSpectrumHisto;

}//close UnfoldSpectrum()


void ForwardFolder::MinimizedFcn(int& nPar, double* gin, double &f, double* par, int iflag)
{
	//## The parameters are the correction factors to be applied to the measured spectrum 
	//## Set correction in current true spectrum
	f = 0.;
	double LogLikelihood= 0.;
	double Chi2= 0.;
	double NDF= 0.;
	int nTotPar= fInitFitPars->GetNPars();//total number of parameters (nPar is the free number of pars)

	//## Update "true" spectrum given current fit parameters
	fTrueSpectrumModelFcn->SetParameters(par);
	//fTrueSpectrumModelFcn->SetParameter(0,pow(10,par[0]));
	
	SpectrumUtils::GetModelSpectrum(*fTrueSpectrum,fTrueSpectrumModelFcn,false);

	//## Forward fold the "true" spectrum with resolution model
	TH1D* FoldedSpectrum= SpectrumUtils::GetFoldedSpectrum(fTrueSpectrum,fResponseMatrix);


	//fForwardFoldedSpectrum->Reset();
	fCurrentForwardFoldedSpectrum->Reset();

	for(int i=0;i<fCurrentRecSpectrum->GetNbinsX();i++){
		double lgERec= fCurrentRecSpectrum->GetBinCenter(i+1);
		double lgERecMin= fCurrentRecSpectrum->GetBinLowEdge(i+1);
		double lgERecMax= lgERecMin + fCurrentRecSpectrum->GetBinWidth(i+1);
		
		double nData= fCurrentRecSpectrum->GetBinContent(i+1);
		double nDataErr= fCurrentRecSpectrum->GetBinError(i+1);
		if(nDataErr==0) nDataErr= 1;
		
		//Find corresponding bin in folded spectrum
		int binId= FoldedSpectrum->FindBin(lgERec);
		if(binId==-1){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot find bin (lgERec="<<lgERec<<") in folded spectrum, skip bin!");
			#endif
			continue;
		}
		double nExp= FoldedSpectrum->GetBinContent(binId);
		cout<<"INFO: lgERec="<<lgERec<<", nData="<<nData<<" +- "<<nDataErr<<", nExp="<<nExp<<endl;
	
		//fForwardFoldedSpectrum->SetBinContent(i+1,nExp);
		//fForwardFoldedSpectrum->SetBinError(i+1,sqrt(nExp));
		fCurrentForwardFoldedSpectrum->SetBinContent(i+1,nExp);
		fCurrentForwardFoldedSpectrum->SetBinError(i+1,sqrt(nExp));

		
		//Skip bins if requested
		if(nData==0) continue;
		if(fUseFitRange && (lgERecMin<fLgEMin_fit || lgERecMax>fLgEMax_fit) ) continue;

		
		if(fIsLikelihoodFit){
			double thisChi2= 0;
			if(nExp>0){
  			LogLikelihood+= -(nData*log(nExp)-nExp);
				if(nData>0){
					Chi2+= nExp-nData + nData*log(nData/nExp);
					thisChi2= nExp-nData + nData*log(nData/nExp);
				}
				else{
					Chi2+= nExp;
					thisChi2= nExp;
				}
				cout<<"INFO: lgERec="<<lgERec<<", nData="<<nData<<" +- "<<nDataErr<<", nExp="<<nExp<<", LL="<<LogLikelihood<<", Chi2="<<Chi2<<", DeltaChi2="<<thisChi2<<endl;
			}
			else{
				LogLikelihood+= 1.e+99;
			}
		}
		else{
			double diff= (nExp-nData);
			if(fUseErrorsInChi2 && nDataErr>0) diff/= nDataErr;
			double thisChi2= diff*diff;
			Chi2+= thisChi2;
			cout<<"INFO: lgERec="<<lgERec<<", nData="<<nData<<" +- "<<nDataErr<<", nExp="<<nExp<<", thisChi2="<<thisChi2<<", Chi2="<<Chi2<<endl;
		}

		NDF++;

	}//end loop rec bins

	
	if(fIsLikelihoodFit) Chi2*= 2;
	NDF-= nPar-2;

	if(fIsLikelihoodFit){ 
		f = LogLikelihood;
	}
	else{
		f= Chi2;
	}

	cout<<"*** CURRENT FIT ***"<<endl;
	cout<<"nPar="<<nPar<<" nTotPar="<<nTotPar<<endl;
	cout<<"LL="<<LogLikelihood<<endl;
	cout<<"Chi2="<<Chi2<<"  Ndf="<<NDF<<endl;
	cout<<"Chi2/Ndf="<<Chi2/NDF<<endl;
	cout<<"== CURRENT PARS =="<<endl;
	for(int i=0;i<nTotPar;i++){
		cout<<"Par["<<i<<"]="<<par[i]<<endl;
	}
	cout<<"*******************"<<endl;
	cout<<endl;
	
	FoldedSpectrum->Delete();

}//close MinimizedFcn()


int ForwardFolder::ComputeUncertainties(SpectrumPars& initFitPars,int nRandomSamples)
{
	//## Check covariance matrix
	if(!fCovarianceMatrix){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null pointer to covariance matrix, cannot compute uncertainties!");
		#endif
		return -1;
	}
	if(nRandomSamples<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of random sampled given!");
		#endif
		return -1;
	}

	cout<<"--> Fit Covariance Matrix"<<endl;
	fCovarianceMatrix->Print();

	int nPar= fCovarianceMatrix->GetNrows();
	int nFreePar= 0;
	std::vector<bool> isFreeFlags;
	std::vector<int> freeParIndex;
	for(int i=0;i<initFitPars.GetNPars();i++){
		FitPar* thisFitPar= initFitPars.GetPar(i);
		bool isFixed= thisFitPar->IsFixed();
		if(!isFixed) {
			nFreePar++;
			isFreeFlags.push_back(true);
		}	
		else{
			isFreeFlags.push_back(false);
			freeParIndex.push_back(i);
		}
	}

	TMatrixD C(nFreePar,nFreePar);
	double FittedPar[nPar];
	int ix= 0;
	int iy= 0;

	for(int i=0;i<nPar;i++){
		FittedPar[i]= fFitPar[i];

		for(int j=0;j<nPar;j++){
			if(isFreeFlags[i] && isFreeFlags[j]) {
				C(ix,iy)= (*fCovarianceMatrix)(i,j);
				if(j==nPar-1) iy=0;
				else iy++;
			}
		}//end loop y
		if(isFreeFlags[i]) ix++;
	}//end loop x

	cout<<"--> Covariance matrix for free pars (N="<<nFreePar<<")"<<endl;
	C.Print();


	/*
	TMatrixD C(nPar,nPar);
	double FittedPar[nPar];
	for(int i=0;i<nPar;i++){
		FittedPar[i]= fFitPar[i];
		for(int j=0;j<nPar;j++){
			C(i,j)= (*fCovarianceMatrix)(i,j);
		}
	}
	*/

	//## Find a Cholesky decomposition of the covariance matrix for error propagation
	cout<<"|C|="<<C.Determinant()<<endl;
	
	TDecompChol CCholDecomp(C);
  CCholDecomp.Decompose();
	cout<<"--> Choleski Decomposition..."<<endl;	
	CCholDecomp.Print();
	
	//TMatrixD* CCholDecompTriang= new TMatrixD(nPar,nPar);
	TMatrixD* CCholDecompTriang= new TMatrixD(nFreePar,nFreePar);
	CCholDecompTriang= (TMatrixD*)(&CCholDecomp.GetU());
	cout<<"--> Choleski Decomposition Triangular..."<<endl;	
	CCholDecompTriang->Print();

	//TMatrixD* CCholDecompTriangTransp= new TMatrixD(nPar,nPar);
	TMatrixD* CCholDecompTriangTransp= new TMatrixD(nFreePar,nFreePar);
	CCholDecompTriangTransp->Transpose(*CCholDecompTriang);


	//################################
	//##  STATISTICAL UNCERTAINTIES 
	//################################
	double StatSigma[fNTrueBins];
	std::vector< std::vector<double> > nUnfold_Stat;
	//double nUnfold_Stat[nRandomSamples][fNTrueBins];
	double nUnfoldMean_Stat[fNTrueBins];
	int nToys_Stat= 0;
	for(int i=0;i<fNTrueBins;i++) {
		StatSigma[i]= 0.;
		nUnfoldMean_Stat[i]= 0.;
		nUnfold_Stat.push_back( std::vector<double>() );
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("Propagate rec spectrum statistical uncertainties into the unfolded flux ...");
	#endif
	
	for(int k=0;k<nRandomSamples;k++){
		cout<<"== RAND SPECTRUM NO. "<<k+1<<" =="<<endl;	

		//Get a fluctuated version of the rec spectrum
		TH1D* thisFluctuatedHisto= SpectrumUtils::GetFluctuatedSpectrum(fRecSpectrum);
		if(!thisFluctuatedHisto) continue;

		//Copy content to tmp rec spectrum (used in the unfolding fit)
		fCurrentRecSpectrum->Reset();
		if(SpectrumUtils::CopySpectrumContent(thisFluctuatedHisto,fCurrentRecSpectrum)<0) continue;
			
		//Run the unfolding over the fluctuated spectrum
		TH1D* FluctUnfoldedSpectrum_Stat= UnfoldSpectrum(initFitPars,"FIT");
		if(fFitStatus!=0) {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Unfolding run not converged for rndom sample no. "<<k+1<<", skip it!");
			#endif
			continue;
		}

		//Update unfolded stats counts fr uncertainty calculation
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Updating mean counts for forward-folded random spectra for random sample no. "<<k+1<<" ...");
		#endif
		nToys_Stat++;

		for(int i=0;i<fNTrueBins;i++) {
			double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
			double nUnfoldToyMC_Stat= FluctUnfoldedSpectrum_Stat->GetBinContent(i+1);
			if(TMath::IsNaN(nUnfoldToyMC_Stat)){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Unfolded bin "<<i+1<<" for random spectrum sample no. "<<k<<" is invalid (NAN)!");
				#endif
			}
			nUnfoldMean_Stat[i]+= nUnfoldToyMC_Stat;
			//nUnfold_Stat[k][i]= nUnfoldToyMC_Stat;
			nUnfold_Stat[i].push_back(nUnfoldToyMC_Stat);
		}//end loop true bins
				
		if(FluctUnfoldedSpectrum_Stat) FluctUnfoldedSpectrum_Stat->Delete();
	}//end loop toys

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<nToys_Stat<<" unfolding over the rec spectrum performed ...");
	#endif
	
	if(nToys_Stat<=1){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No random unfolding succeeded, cannot compute uncertainty!");
		#endif
		if(CCholDecompTriang) CCholDecompTriang->Delete();
		if(CCholDecompTriangTransp) CCholDecompTriangTransp->Delete();
		return -1;
	}
	

	//## Calculate sigma syst true model
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing statistical uncertainties ...");
	#endif

	for(int i=0;i<fNTrueBins;i++) {
		double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
		nUnfoldMean_Stat[i]/= (double)(nToys_Stat);
			
		//for(int k=0;k<nRandomSamples;k++){
		for(size_t k=0;k<nUnfold_Stat[i].size();k++){
			//StatSigma[i]+= pow(nUnfold_Stat[k][i]-nUnfoldMean_Stat[i],2);	
			StatSigma[i]+= pow(nUnfold_Stat[i][k]-nUnfoldMean_Stat[i],2);
		}
		StatSigma[i]/= (double)(nToys_Stat-1);
		StatSigma[i]= sqrt(StatSigma[i]);
		
		double relBias_Stat= 0.;
		double relSyst_Stat= 0.;
		if(nUnfold>0) {
			relBias_Stat= nUnfoldMean_Stat[i]/nUnfold-1;	
		}
		if(nUnfoldMean_Stat[i]>0){
			relSyst_Stat= StatSigma[i]/nUnfoldMean_Stat[i];
		}
		#ifdef LOGGING_ENABLED
			INFO_LOG("Bin "<<i<<": StatErr="<<StatSigma[i]<<", Bias="<<relBias_Stat<<", RelStatErr="<<relSyst_Stat);
		#endif
		
	}//end loop true bins
		



	//##############################
	//##   FIT PARS SYSTEMATICS 
	//##############################
	#ifdef LOGGING_ENABLED
		INFO_LOG("Propagate true model systematics into the unfolded flux ...");
	#endif

	double SystSigmaTrueModel[fNTrueBins];
	double nUnfold_SystTrueModel[nRandomSamples][fNTrueBins];
	double nUnfoldMean_SystTrueModel[fNTrueBins];
	for(int i=0;i<fNTrueBins;i++) {
		SystSigmaTrueModel[i]= 0.;
		nUnfoldMean_SystTrueModel[i]= 0.;
	}

	for(int k=0;k<nRandomSamples;k++){
		cout<<"== RANDOM SPECTRUM NO. "<<k+1<<" =="<<endl;		

		//## Generate new true model around uncertainties		
		TMatrixD RandModelPar(nFreePar,1);
		for(int i=0;i<nFreePar;i++){
			RandModelPar(i,0)= gRandom->Gaus(0,1);
		}
		TMatrixD FluctModelPar(nFreePar,1);
		FluctModelPar= (*CCholDecompTriangTransp)*RandModelPar;
		

		double newModelPar[nPar];
		for(int i=0;i<nPar;i++) newModelPar[i]= FittedPar[i];
		for(unsigned int i=0;i<freeParIndex.size();i++){
			int parIndex= freeParIndex[i];	
			newModelPar[parIndex]+= FluctModelPar(i,0);
		}
		

		//Update pars
		SpectrumPars* modelFitPars= initFitPars.Clone();
		for(int i=0;i<nPar;i++){
			modelFitPars->SetParValue(i,newModelPar[i]);
		}

		TH1D* FluctUnfoldedSpectrum= UnfoldSpectrum(*modelFitPars,"CALL");

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Calculate mean counts for forward-folded toy spectra for random sample no. "<<k+1<<" ...");
		#endif

		for(int i=0;i<fNTrueBins;i++) {
			double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
			double nUnfoldToyMC= FluctUnfoldedSpectrum->GetBinContent(i+1);
			//cout<<"bin "<<i<<" nUnfold="<<nUnfold<<"  nUnfoldToyMC="<<nUnfoldToyMC<<endl;
			nUnfoldMean_SystTrueModel[i]+= nUnfoldToyMC;
			nUnfold_SystTrueModel[k][i]= nUnfoldToyMC;
			
		}//end loop true bins

		if(modelFitPars){
			delete modelFitPars;
			modelFitPars= 0;
		}

	}//end loop toys


	//## Calculate sigma syst true model
	#ifdef LOGGING_ENABLED
		INFO_LOG("Calculate sigma syst for true model...");
	#endif

	for(int i=0;i<fNTrueBins;i++) {
		double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
		nUnfoldMean_SystTrueModel[i]/= (double)(nRandomSamples);
		for(int k=0;k<nRandomSamples;k++){
			SystSigmaTrueModel[i]+= pow(nUnfold_SystTrueModel[k][i]-nUnfoldMean_SystTrueModel[i],2); 	
		}
		SystSigmaTrueModel[i]/= (double)(nRandomSamples-1);
		SystSigmaTrueModel[i]= sqrt(SystSigmaTrueModel[i]);
		double relBias= 0.;
		double relSyst= 0.;
		if(nUnfold>0) {
			relBias= nUnfoldMean_SystTrueModel[i]/nUnfold-1;
		}
		if(nUnfoldMean_SystTrueModel[i]>0){
			relSyst= SystSigmaTrueModel[i]/nUnfoldMean_SystTrueModel[i];
		}
		#ifdef LOGGING_ENABLED
			INFO_LOG("Bin "<<i<<": Syst(TrueModel)="<<SystSigmaTrueModel[i]<<", Bias(TrueModel)="<<relBias<<", RelSyst(TrueModel)="<<relSyst);
		#endif
	}
	
	
	//## Forward-folding results
	#ifdef LOGGING_ENABLED
		INFO_LOG("Final unfolding results...");
	#endif
	
	for(int i=0;i<fUnfoldedSpectrum->GetNbinsX();i++){
		double binCenter= fUnfoldedSpectrum->GetBinCenter(i+1);
		double binContent= fUnfoldedSpectrum->GetBinContent(i+1);
		int recBinId= fRecSpectrum->FindBin(binCenter);
		double nRec= fRecSpectrum->GetBinContent(recBinId);			
		double nRecErr= fRecSpectrum->GetBinError(recBinId);			
		
		double totErrorSqr= 0;

		//Add stat uncertainties
		double statError_rec= fRecSpectrum->GetBinError(recBinId);
		double statError= StatSigma[i];
		totErrorSqr+= statError*statError;
	
		//Add syst uncertainty due to unfolding fit
		double systError_UnfoldingAlgo= SystSigmaTrueModel[i];
		totErrorSqr+= systError_UnfoldingAlgo*systError_UnfoldingAlgo;

	
		double totError= sqrt(totErrorSqr);
		double binError= statError;	
			
		fUnfoldedSpectrum->SetBinContent(i+1,binContent);
		fUnfoldedSpectrum->SetBinError(i+1,binError);	
		cout<<"bin "<<i+1<<"  nRec="<<nRec<<"  +- "<<nRecErr<<"("<<nRecErr/nRec<<")  nUnfold="<<binContent<<"  +- "<<binError<<"("<<binError/binContent<<")  (recStatErr="<<statError_rec<<"), syst(unfoldAlgo)="<<systError_UnfoldingAlgo<<", totError="<<totError<<endl;
	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Bin "<<i+1<<": nRec="<<nRec<<"  +- "<<nRecErr<<"("<<nRecErr/nRec<<"). nUnfold="<<binContent<<"  +- "<<binError<<"("<<binError/binContent<<")  (recStatErr="<<statError_rec<<"), syst(unfoldAlgo)="<<systError_UnfoldingAlgo<<", totError="<<totError);
		#endif
	}//end loop true bins
	
	if(CCholDecompTriang) CCholDecompTriang->Delete();
	if(CCholDecompTriangTransp) CCholDecompTriangTransp->Delete();
	
	return 0;

}//close ComputeUncertainties()


