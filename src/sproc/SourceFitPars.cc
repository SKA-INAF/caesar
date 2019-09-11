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
* @file SourceFitPars.cc
* @class SourceFitPars
* @brief SourceFitPars
*
* Source fit parameters class
* @author S. Riggi
* @date 01/09/2017
*/

#include <SourceFitPars.h>
#include <SourceComponentPars.h>

#include <Image.h>
#include <Source.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <SysUtils.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <Consts.h>
#include <WCSUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TFitResult.h>
#include <TBackCompFitter.h>

//C++ headers
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
#include <queue>
#include <chrono>

using namespace std;
using namespace std::placeholders;

ClassImp(Caesar::SourceFitPars)

namespace Caesar {

//Constructor
SourceFitPars::SourceFitPars()
	: TObject()
{
	Init();

}//close costructor

//Parametric constructor
SourceFitPars::SourceFitPars(int N)
	: TObject()
{
	Init();
	SetNComponents(N);
}	

//Destructor
SourceFitPars::~SourceFitPars()
{
	
}//close destructor

//Copy constructor
SourceFitPars::SourceFitPars(const SourceFitPars& sourceFitPars)
{
	((SourceFitPars&)sourceFitPars).Copy(*this);
}

//Assignment operator
SourceFitPars& SourceFitPars::operator=(const SourceFitPars& sourceFitPars)
{
	// Operator =
 	if (this != &sourceFitPars) ((SourceFitPars&)sourceFitPars).Copy(*this);
  return *this;
}

//Copy method
void SourceFitPars::Copy(TObject& obj) const 
{
	// Copy this source to source obj	
  ((SourceFitPars&)obj).nComponents = nComponents;
	((SourceFitPars&)obj).chi2 = chi2;	
	((SourceFitPars&)obj).ndof = ndof;
	((SourceFitPars&)obj).npars = npars;	
	((SourceFitPars&)obj).npars_free = npars_free;	
	((SourceFitPars&)obj).npars_component= npars_component;
	((SourceFitPars&)obj).nfit_points = nfit_points;	
	((SourceFitPars&)obj).status = status;	
	((SourceFitPars&)obj).minimizer_status = minimizer_status;	
	((SourceFitPars&)obj).offset = offset;	
	((SourceFitPars&)obj).offset_err = offset_err;	
	((SourceFitPars&)obj).residualMean = residualMean;	
	((SourceFitPars&)obj).residualRMS = residualRMS;	
	((SourceFitPars&)obj).residualMedian = residualMedian;	
	((SourceFitPars&)obj).residualMAD = residualMAD;	
	((SourceFitPars&)obj).residualMin = residualMin;	
	((SourceFitPars&)obj).residualMax = residualMax;	

	((SourceFitPars&)obj).pars = pars;	
	((SourceFitPars&)obj).thetaFixed = thetaFixed;	
	((SourceFitPars&)obj).offsetFixed = offsetFixed;	
	((SourceFitPars&)obj).sigmaFixed = sigmaFixed;	
	((SourceFitPars&)obj).fluxDensity = fluxDensity;	
	((SourceFitPars&)obj).fluxDensityErr = fluxDensityErr;	

	((SourceFitPars&)obj).fitQuality = fitQuality;
		
	//Copy matrix
	int nRows= fitCovarianceMatrix.GetNrows();
	int nCols= fitCovarianceMatrix.GetNcols();
	((SourceFitPars&)obj).fitCovarianceMatrix.ResizeTo(nRows,nCols);
	for(int i=0;i<nRows;i++){
		for(int j=0;j<nCols;j++){
			double w= fitCovarianceMatrix(i,j);
			(((SourceFitPars&)obj).fitCovarianceMatrix)(i,j)= w;
		}
	}
		
	nRows= fluxDensityDerivMatrix.GetNrows();
	nCols= fluxDensityDerivMatrix.GetNcols();
	((SourceFitPars&)obj).fluxDensityDerivMatrix.ResizeTo(nRows,nCols);
	for(int i=0;i<nRows;i++){
		for(int j=0;j<nCols;j++){
			double w= fluxDensityDerivMatrix(i,j);
			(((SourceFitPars&)obj).fluxDensityDerivMatrix)(i,j)= w;
		}
	}
		
}//close Copy()


void SourceFitPars::Init()
{
	nComponents= 0;
	chi2= 0;
	ndof= 0;
	npars= 0;
	npars_free= 0;	
	npars_component= 6;
	nfit_points= 0;
	status= -1;
	minimizer_status= -1;
	offset= 0;
	offset_err= 0;
	residualMean= 0;
	residualRMS= 0;
	residualMedian= 0;
	residualMAD= 0;
	residualMin= 0;
	residualMax= 0;
	pars.clear();
	thetaFixed= false;
	offsetFixed= false;
	sigmaFixed= false;	
	fluxDensity= 0;
	fluxDensityErr= 0;
	fitQuality= eBadFit;
		
}//close Init()

int SourceFitPars::RemoveComponents(std::vector<int> componentIds)
{
	//Check components given
	int nComponentsToRemove= static_cast<int>(componentIds.size());
	if(nComponentsToRemove<=0 || nComponentsToRemove>nComponents){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No components given to remove or too many (exceeding number of current components)!");
		#endif
		return -1;
	}
	for(size_t i=0;i<componentIds.size();i++){
		int componentId= componentIds[i];
		if(componentId<0 || componentId>=nComponents) {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Component "<<componentId<<" does not exist!");
			#endif
			return -1;
		}
	}//end loop

	//Delete fit components from vector
	CodeUtils::DeleteItems(pars,componentIds);

	//Recompute nComponents
	nComponents= static_cast<int>(pars.size());

	return 0;

}//close RemoveComponents()

double SourceFitPars::GetComponentFluxDensity(int componentId)
{
	if(componentId<0 || componentId>=nComponents) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Component "<<componentId<<" does not exist, returning zero flux!");
		#endif
		return 0;
	}

	return pars[componentId].GetFluxDensity();

}//close GetComponentFluxDensity()

double SourceFitPars::GetComponentFluxDensityErr(int componentId)
{
	if(componentId<0 || componentId>=nComponents){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Component "<<componentId<<" does not exist, returning zero error flux!");
		#endif
		return 0;
	}
	//Get component flux density deriv matrix
	TMatrixD D;
	if(GetComponentFluxDerivMatrix(D,componentId)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute component flux derivative matrix, returning zero error flux!");
		#endif
		return 0;
	}

	TMatrixD D_t= TMatrixD(TMatrixD::kTransposed,D);
	TMatrixD VarMatrix= D*fitCovarianceMatrix*D_t;
	double Var= VarMatrix(0,0);
	if(Var<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Flux density variance for component "<<componentId<<" is negative (this should not occur, check for bugs or numerical roundoff errors!)");
		#endif
		return 0;
	}
	double Err= sqrt(Var);

	//Convert to Jy
	Err/= 1.e+3;

	return Err;			

}//close GetComponentFluxDensityErr()


}//close namespace
