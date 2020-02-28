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
* @file SourceFitPars.h
* @class SourceFitPars
* @brief SourceFitPars
*
* Source fit parameters class
* @author S. Riggi
* @date 01/09/2017
*/

#ifndef _SOURCE_FIT_PARS_h
#define _SOURCE_FIT_PARS_h 1

#include <SourceComponentPars.h>

#include <SysUtils.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <AstroUtils.h>
#include <WCSUtils.h>
#include <Consts.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TFitResultPtr.h>
#include <Fit/FitResult.h>
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TMinuitMinimizer.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>


namespace Caesar {


//========================================
//==         SOURCE FIT PARS
//========================================
class SourceFitPars : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFitPars();
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFitPars(int N);

		/**
		* \brief Copy constructor
		*/
		SourceFitPars(const SourceFitPars& sourceFitPars);

		/**
		* \brief Class destructor: free allocated memory
		*/		
		virtual ~SourceFitPars();
		
		
		/**
		* \brief Assignment Operator
		*/
		SourceFitPars& operator=(const SourceFitPars& sourceFitPars);

		/**
		* \brief Copy method
		*/
		void Copy(TObject& obj) const;

	public:
		/**
		* \brief Set offset par value 
		*/
		void SetOffsetPar(double value){offset=value;}
		/**
		* \brief Get offset par value 
		*/
		double GetOffsetPar(){return offset;}

		/**
		* \brief Set offset par error
		*/
		void SetOffsetParErr(double value){offset_err=value;}
		/**
		* \brief Get offset par error value 
		*/
		double GetOffsetParErr(){return offset_err;}

		/**
		* \brief Get fitted ellipses
		*/
		std::vector<TEllipse*> GetFittedEllipses(bool useFWHM=true){
			std::vector<TEllipse*> ellipses;
			for(size_t i=0;i<pars.size();i++){
				TEllipse* ellipse= pars[i].GetFitEllipse(useFWHM);
				ellipses.push_back(ellipse);
			}
			return ellipses;
		}

		/**
		* \brief Get pars
		*/
		std::vector<SourceComponentPars>const& GetPars() const {return pars;}

		/**
		* \brief Set pars for component i-th
		*/
		int SetComponentPars(int componentId,SourceComponentPars& componentPars){
			if(componentId<0 || componentId>=(signed)pars.size()) return -1;
			pars[componentId]= componentPars;
			return 0;
		}

		/**
		* \brief Set par value & error
		*/
		int SetParValueAndError(int componentId,std::string parName,double parValue,double parError){
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Invalid component id ("<<componentId<<" given, cannot find par to be set!");
				#endif
				return -1;
			}
			return pars[componentId].SetParValueAndError(parName,parValue,parError);
		}

		/** 
		\brief Get par value
 		*/
		double GetParValue(int componentId,std::string parName){
			if(componentId<0 || componentId>=nComponents) return -999;
			return pars[componentId].GetParValue(parName);
		}

		/** 
		\brief Get par error
 		*/
		double GetParError(int componentId,std::string parName){
			if(componentId<0 || componentId>=nComponents) return -999;
			return pars[componentId].GetParError(parName);
		}

		/**
		* \brief Initialize component pars
		*/
		void SetNComponents(int N){
			pars.clear();
			for(int i=0;i<N;i++) pars.push_back(SourceComponentPars());
			nComponents= N;
		}
	
		/**
		* \brief Get number of components
		*/
		int GetNComponents(){return nComponents;}

		/**
		* \brief Get number of selected components
		*/
		int GetNSelComponents()
		{
			int nSelComponents= 0;
			for(size_t i=0;i<pars.size();i++){
				bool isSelected= pars[i].IsSelected();
				if(isSelected) nSelComponents++;
			}
			return nSelComponents;
		}
	
		/**
		* \brief Check if component is selected
		*/
		bool IsSelectedComponent(int componentId){
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning false!");
				#endif
				return false;
			}
			return pars[componentId].IsSelected();
		}

		/**
		* \brief Set if component is selected
		*/
		int SetSelectedComponent(int componentId,bool selected)
		{
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning false!");
				#endif
				return -1;
			}
			pars[componentId].SetSelected(selected);
			return 0;
		}

		/**
		* \brief Remove fit components
		*/
		int RemoveComponents(std::vector<int> componentIds);

		/**
		* \brief Get component position
		*/
		int GetComponentPosition(double& xpos,double& ypos,int componentId)
		{
			xpos= 0;
			ypos= 0;
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning zero pos!");
				#endif
				return -1;
			}
			pars[componentId].GetPosition(xpos,ypos);
			return 0;
		}

		/**
		* \brief Get component peak flux
		*/
		double GetComponentPeakFlux(int componentId)
		{
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning zero peak flux!");
				#endif
				return 0;
			}
			return pars[componentId].GetPeakFlux();
		}

		/**
		* \brief Get component flux density
		*/
		double GetComponentFluxDensity(int componentId);

		/**
		* \brief Get component flux density error
		*/
		double GetComponentFluxDensityErr(int componentId);
		/**
		* \brief Get component fit beam ellipse pars
		*/
		int GetComponentBeamEllipsePars(int componentId,double& bmaj,double& bmin,double& pa)
		{
			//Init values
			bmaj= 0;
			bmin= 0;
			pa= 0;
			if(componentId<0 || componentId>=nComponents){	
				#ifdef LOGGING_ENABLED	
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetBeamEllipsePars(bmaj,bmin,pa);
		}

		/**
		* \brief Get component beam ellipse eccentricity
		*/
		double GetComponentBeamEllipseEccentricity(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning E=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetBeamEllipseEccentricity();
		}

		/**
		* \brief Get component beam ellipse area
		*/
		double GetComponentBeamEllipseArea(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning A=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetBeamEllipseArea();
		}

		/**
		* \brief Has component beam pars stored?
		*/
		bool HasComponentBeamEllipsePars(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){	
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning false!");
				#endif
				return 0;
			}
			return pars[componentId].HasBeamEllipsePars();		
		}

		/**
		* \brief Get component fit ellipse pars
		*/
		int GetComponentFitEllipsePars(int componentId,double& x0,double& y0,double& bmaj,double& bmin,double& pa)
		{
			//Init values
			x0= 0;
			y0= 0;
			bmaj= 0;
			bmin= 0;
			pa= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetEllipsePars(x0,y0,bmaj,bmin,pa);
		}

		/**
		* \brief Get component fit ellipse eccentricity
		*/
		double GetComponentFitEllipseEccentricity(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning E=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetEllipseEccentricity();
		}

		/**
		* \brief Get component fit ellipse area
		*/
		double GetComponentFitEllipseArea(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning A=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetEllipseArea();
		}

		/**
		* \brief Get component fit ellipse rot angle vs beam
		*/
		double GetComponentFitEllipseRotAngleVSBeam(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED				
					WARN_LOG("Component "<<componentId<<" does not exist, returning A=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetEllipseRotAngleVSBeam();
		}

		/**
		* \brief Get component fit ellipse pars
		*/
		int GetComponentFitEllipseParErrors(int componentId,double& x0_err,double& y0_err,double& bmaj_err,double& bmin_err,double& pa_err)
		{
			//Init values
			x0_err= 0;
			y0_err= 0;
			bmaj_err= 0;
			bmin_err= 0;
			pa_err= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetEllipseParErrors(x0_err,y0_err,bmaj_err,bmin_err,pa_err);
		}

		/**
		* \brief Get component WCS fit ellipse pars
		*/
		int GetComponentFitWCSEllipsePars(int componentId,double& x0_wcs,double& y0_wcs,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
		{
			//Init values
			x0_wcs= 0;
			y0_wcs= 0;
			bmaj_wcs= 0;
			bmin_wcs= 0;
			pa_wcs= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetWCSEllipsePars(x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs);
		}

		/**
		* \brief Get component WCS fit ellipse par errors
		*/
		int GetComponentFitWCSEllipseParErrors(int componentId,double& x0_wcs_err,double& y0_wcs_err,double& bmaj_wcs_err,double& bmin_wcs_err,double& pa_wcs_err)
		{
			//Init values
			x0_wcs_err= 0;
			y0_wcs_err= 0;
			bmaj_wcs_err= 0;
			bmin_wcs_err= 0;
			pa_wcs_err= 0;
			if(componentId<0 || componentId>=nComponents){	
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetWCSEllipseParErrors(x0_wcs_err,y0_wcs_err,bmaj_wcs_err,bmin_wcs_err,pa_wcs_err);
		}

		/**
		* \brief Get component WCS fit ellipse pars
		*/
		int GetComponentFitWCSDeconvolvedEllipsePars(int componentId,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
		{
			//Init values
			bmaj_wcs= 0;
			bmin_wcs= 0;
			pa_wcs= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetWCSDeconvolvedEllipsePars(bmaj_wcs,bmin_wcs,pa_wcs);
		}

		/**
		* \brief Get flux density
		*/
		double GetFluxDensity(){return fluxDensity;}

		/**
		* \brief Get flux density error
		*/
		double GetFluxDensityErr(){return fluxDensityErr;}

		/**
		* \brief Set flux density
		*/
		void SetFluxDensity(double value){fluxDensity=value;}
		/**
		* \brief Set flux density error
		*/
		void SetFluxDensityErr(double value){fluxDensityErr=value;}

		/**
		* \brief Set chi2 
		*/
		void SetChi2(double value){chi2=value;}

		/**
		* \brief Get chi2 
		*/
		double GetChi2(){return chi2;}
		/**
		* \brief Set NDF 
		*/
		void SetNDF(double value){ndof=value;}
		/**
		* \brief Get NDF 
		*/
		double GetNDF(){return ndof;}

		/**
		* \brief Set Status 
		*/
		void SetStatus(int value){status=value;}

		/**
		* \brief Get Status 
		*/
		int GetStatus(){return status;}

		/**
		* \brief Set Minimizer Status 
		*/
		void SetMinimizerStatus(int value){minimizer_status=value;}

		/**
		* \brief Get Minimizer Status 
		*/
		int GetMinimizerStatus(){return minimizer_status;}

		/**
		* \brief Set number of free pars 
		*/
		void SetNFreePars(int value){npars_free=value;}
		/**
		* \brief Get number of free parameters 
		*/
		int GetNFreePars(){return npars_free;}

		/**
		* \brief Set number of component fit parameters 
		*/
		void SetNComponentPars(int nc){npars_component=nc;}
		/**
		* \brief Get number of component fit parameters 
		*/
		int GetNComponentPars(){return npars_component;}


		
		/**
		* \brief Set total number of pars 
		*/
		void SetNPars(int value){npars=value;}
		/**
		* \brief Get number of free parameters 
		*/
		int GetNPars(){return npars;}

		/**
		* \brief Set number of fitted data 
		*/
		void SetNFitPoints(int value){nfit_points=value;}

		/**
		* \brief Get number of fitted data 
		*/
		int GetNFitPoints(){return nfit_points;}

		/**
		* \brief Get residual mean
		*/
		double GetResidualMean(){return residualMean;}
		/**
		* \brief Set residual mean
		*/
		void SetResidualMean(double value){residualMean=value;}
		/**
		* \brief Get residual rms
		*/
		double GetResidualRMS(){return residualRMS;}
		/**
		* \brief Set residual rms
		*/
		void SetResidualRMS(double value){residualRMS=value;}
		/**
		* \brief Get residual median
		*/
		double GetResidualMedian(){return residualMedian;}
		/**
		* \brief Set residual median
		*/
		void SetResidualMedian(double value){residualMedian=value;}
		/**
		* \brief Get residual mad
		*/
		double GetResidualMAD(){return residualMAD;}
		/**
		* \brief Set residual mad
		*/
		void SetResidualMAD(double value){residualMAD=value;}
		/**
		* \brief Get residual min
		*/
		double GetResidualMin(){return residualMin;}
		/**
		* \brief Set residual min
		*/
		void SetResidualMin(double value){residualMin=value;}		
		/**
		* \brief Get residual max
		*/
		double GetResidualMax(){return residualMax;}
		/**
		* \brief Set residual max
		*/
		void SetResidualMax(double value){residualMax=value;}
		/**
		* \brief Set fit covariance matrix
		*/
		int SetCovarianceMatrix(double* errMatrixValues,int ndim)
		{
			if(!errMatrixValues || ndim<=0) return -1;
			fitCovarianceMatrix.ResizeTo(ndim,ndim);
			fitCovarianceMatrix.Zero();
			try{
				for(int i=0;i<ndim;i++){
					for(int j=0;j<ndim;j++){
						int index= i*ndim+j;
						fitCovarianceMatrix(i,j)= errMatrixValues[index];
					}//end loop dim
				}//end loop dim
			}//close try block
			catch(...){	
				#ifdef LOGGING_ENABLED
					ERROR_LOG("C++ exception occurred while filling fit covariance matrix (hint: array size and given dim are different?");
				#endif
				fitCovarianceMatrix.ResizeTo(0,0);
				return 0;
			}
			return 0;
		}//close SetCovarianceMatrix()

		/**
		* \brief Set fit covariance matrix
		*/
		int SetCovarianceMatrix(TMatrixD& M)
		{
			int nCols= M.GetNcols();
			int nRows= M.GetNrows();
			if(nRows<=0 || nCols<=0) return -1;
			fitCovarianceMatrix.ResizeTo(nRows,nCols);
			fitCovarianceMatrix= M;
			return 0;
		}

		/**
		* \brief Remove component(s) from covariance matrix
		*/
		/*
		int RemoveComponentsFromCovarianceMatrix(std::vector<int> componentIds)
		{	
			//Check input
			if(componentIds.empty()) return 0;
			int nComponents= (fitCovarianceMatrix.GetNrows()-1)/npars_component;
			
			for(size_t i=0;i<componentIds.size();i++){
				if(componentIds[i]<0 || componentIds[i]>=nComponents){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Invalid component id ("<<componentIds[i]<<") to be removed (#"<<nComponents<<" present in cov matrix)!");
					#endif
					return -1;
				}
			}
		
			//Set components to be kept
			std::vector<int> keptComponentIds;
			for(int i=0;i<nComponents;i++){
				bool keepComponent= true;
				for(size_t k=0;k<componentIds.size();k++){
					if(i==componentIds[k]){
						keepComponent= false;
						break;
					}
				}
				if(keepComponent) keptComponentIds.push_back(i);
			}

			//Sort component ids to be kept
			std::sort(keptComponentIds.begin(),keptComponentIds.end());

			//Check if there are components left
			int nComponents_left= static_cast<int>(keptComponentIds.size());	
			int nDim_new= nComponents_left*npars_component + 1;
			if(nComponents_left<=0){
				fitCovarianceMatrix.ResizeTo(0,0);
				return 0;
			}

			//Copy cov matrix
			TMatrixD C= fitCovarianceMatrix;
			fitCovarianceMatrix.ResizeTo(nDim_new,nDim_new);
			for(int i=0;i<fitCovarianceMatrix.GetNrows();i++){
				
				for(int j=0;j<fitCovarianceMatrix.GetNcols();j++){
					int componentIndex= j/npars_component;
					int parIndex= j%npars_component;
					int componentId= keptComponentIds[componentIndex];
					int colId= componentId*npars_component + parIndex;
					int rowId= componentId*npars_component + j;
					//fitCovarianceMatrix(i,j)= C(colId,colId);
				}
			}

			for(int k=0;k<nComponents_left;k++){
				int componentId= keptComponentIds[k];
				int row_start= componentId*npars_component;
				for(int i=0;i<npars_component;i++){
					for(int j=0;j<npars_component;j++){
						int index= i*npars_component + j;
						
						//fitCovarianceMatrix(i,j)= errMatrixValues[index];
					}
				}

				for(int i=0;i<fitCovarianceMatrix.GetNrows();i++){
					for(int j=0;j<fitCovarianceMatrix.GetNcols();j++){
						int componentId= 
						int rowId= i*npars_component + j;
						
						//fitCovarianceMatrix(i,j)= errMatrixValues[index];
					}
				}
			}

			return 0;

		}//close RemoveComponentsFromCovarianceMatrix()
		*/
		

		/**
		* \brief Get fit covariance matrix
		*/	
		TMatrixD& GetCovarianceMatrix(){return fitCovarianceMatrix;}

		/**
		* \brief Print fit covariance matrix
		*/
		void PrintCovarianceMatrix(){
			fitCovarianceMatrix.Print();
		}

		/**
		* \brief Get flux density derivative matrix
		*/	
		TMatrixD& GetFluxDensityDerivMatrix(){return fluxDensityDerivMatrix;}

		/**
		* \brief Set flux density derivative matrix
		*/
		int SetFluxDensityDerivMatrix(TMatrixD& M)
		{	
			int nCols= M.GetNcols();
			int nRows= M.GetNrows();
			if(nRows<=0 || nCols<=0) return -1;
			fluxDensityDerivMatrix.ResizeTo(nRows,nCols);
			fluxDensityDerivMatrix= M;
			return 0;
		}


		/**
		* \brief Print flux density derivative matrix
		*/
		void PrintFluxDensityDerivMatrix(){
			fluxDensityDerivMatrix.Print();
		}

		/**
		* \brief Compute flux density derivative matrix
		*/
		int ComputeFluxDensityDerivMatrix(){
			//Check size
			//if(npars_free<=0 || pars.empty()) {
			if(npars<=0 || pars.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Cannot compute derivative matrix as no fit pars are stored and/or number of free pars is not initialized!");
				#endif
				return -1;
			}

			//Init matrix to 1 x Nfree_pars
			//fluxDensityDerivMatrix.ResizeTo(1,npars_free);
			fluxDensityDerivMatrix.ResizeTo(1,npars);
			fluxDensityDerivMatrix.Zero();

			//Fill matrix 
			//NB: consider the cases in which sigma/offset/theta are fixed
			int parCounter= 0;
			for(int i=0;i<nComponents;i++){
				//Retrieve fitted pars for this component
				double A= pars[i].GetParValue("A");
				//A*= 1.e+3;//convert to mJy
				A*= normFactor;//convert to mJy
				double sigmaX= pars[i].GetParValue("sigmaX");
				double sigmaY= pars[i].GetParValue("sigmaY");

				//Fill derivative wrt to amplitude
				//deriv= 2 pi sigmaX sigmaY
				double deriv_ampl= 2.*TMath::Pi()*sigmaX*sigmaY;
				fluxDensityDerivMatrix(0,parCounter)= deriv_ampl;
				parCounter++;

				//Fill derivative wrt to centroids
				//deriv= 0 (for both x & y)
				fluxDensityDerivMatrix(0,parCounter)= 0;
				fluxDensityDerivMatrix(0,parCounter+1)= 0;
				parCounter+= 2;

				//Fill derivative wrt to sigmas (if not fixed)
				//deriv_x= 2 pi A sigma_y
				//deriv_y= 2 pi A sigma_x
				if(!sigmaFixed){
					double deriv_sigmaX= 2.*TMath::Pi()*A*sigmaY;
					double deriv_sigmaY= 2.*TMath::Pi()*A*sigmaX;
					fluxDensityDerivMatrix(0,parCounter)= deriv_sigmaX;
					fluxDensityDerivMatrix(0,parCounter+1)= deriv_sigmaY;
					parCounter+= 2;
				}

				//Fill derivative wrt to theta (if not fixed)
				//deriv_theta= 0;
				if(!thetaFixed){
					fluxDensityDerivMatrix(0,parCounter)= 0;
					parCounter++;
				}

			}//end loop components	

			//Fill deriv wrt to offset (if not fixed)
			//deriv_offset= 0;
			if(!offsetFixed){
				fluxDensityDerivMatrix(0,parCounter)= 0;
			}
	
			return 0;
		}//close ComputeFluxDensityDerivMatrix()
		
		/**
		* \brief Compute flux density
		*/
		void ComputeFluxDensity(){
			//Summing up flux for all components
			fluxDensity= 0;
			for(int i=0;i<nComponents;i++){
				double componentFluxDensity= pars[i].GetFluxDensity();
				fluxDensity+= componentFluxDensity;
			}//end loop components
		}//close ComputeFluxDensity()

		/**
		* \brief Compute flux density error
		*/
		int ComputeFluxDensityError(){
			//Do some checks on covariance & deriv matrix
			fluxDensityErr= 0;
			int nRows= fitCovarianceMatrix.GetNrows();
			int nCols= fluxDensityDerivMatrix.GetNcols();
			if(nCols<=0 || nRows<=0 || nRows!=nCols){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Fit covariance and/or deriv matrix were not computed or have invalid dimensions!");
				#endif
				return -1;
			}
	
			//Compute fluxDensityVariance= D x CovMatrix x D^t  (D=deriv matrix)
			TMatrixD fluxDensityDerivMatrix_t= TMatrixD(TMatrixD::kTransposed,fluxDensityDerivMatrix);
			TMatrixD fluxDensityVarianceMatrix= fluxDensityDerivMatrix*fitCovarianceMatrix*fluxDensityDerivMatrix_t;
			double fluxDensityVariance= fluxDensityVarianceMatrix(0,0);
			if(fluxDensityVariance<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Flux density variance is negative (this should not occur, check for bugs or numerical roundoff errors!)");
				#endif
				return -1;
			}
			fluxDensityErr= sqrt(fluxDensityVariance);

			//Convert back to Jy
			//fluxDensityErr/= 1.e+3;
			fluxDensityErr/= normFactor;
			
			return 0;
		}//close ComputeFluxDensityError()


		/**
		* \brief Compute component ellipse pars
		*/
		int ComputeComponentEllipsePars()
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeEllipsePars()<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentEllipsePars()


		/**
		* \brief Compute component WCS ellipse pars
		*/
		int ComputeComponentWCSEllipsePars(WCS* wcs)
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeWCSEllipsePars(wcs)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute WCS ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentWCSEllipsePars()

		/**
		* \brief Compute component WCS ellipse pars
		*/
		int ComputeComponentWCSEllipseParsSimple(WCS* wcs)
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeWCSEllipseParsSimple(wcs)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute WCS ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentWCSEllipseParsSimple()

		
		/**
		* \brief Set component beam ellipse pars
		*/
		void SetComponentBeamEllipsePars(double bmaj,double bmin,double pa)
		{
			for(size_t i=0;i<pars.size();i++){
				pars[i].SetBeamEllipsePars(bmaj,bmin,pa);
			}
		}


		/**
		* \brief Compute component WCS ellipse pars
		*/
		int ComputeComponentWCSDeconvolvedEllipsePars()
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeWCSDeconvolvedEllipsePars()<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute WCS ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentWCSDeconvolvedEllipsePars()

		/**
		* \brief Print
		*/
		void Print(){
			cout<<"*** FIT RESULTS ***"<<endl;
			cout<<"nPars="<<npars<<", nParsFree="<<npars_free<<", fitStatus="<<status<<" (minimizer status="<<minimizer_status<<"), fitQuality="<<fitQuality<<endl;
			cout<<"Chi2="<<chi2<<", ndf="<<ndof<<", Chi2/NDF="<<chi2/ndof<<endl;
			cout<<"fluxDensity="<<fluxDensity<<" +- "<<fluxDensityErr<<endl;
			for(int i=0;i<nComponents;i++){
				cout<<"--> Component "<<i+1<<endl;
				cout<<"fluxDensity="<<GetComponentFluxDensity(i)<<" +- "<<GetComponentFluxDensityErr(i)<<endl;
				cout<<"A="<<pars[i].GetParValue("A")<<" +- "<<pars[i].GetParError("A")<<endl;
				cout<<"(x0,y0)=("<<	pars[i].GetParValue("x0")<<","<<pars[i].GetParValue("y0")<<") err("<<pars[i].GetParError("x0")<<","<<pars[i].GetParError("y0")<<")"<<endl;
				cout<<"(sigmaX,sigmaY)=("<<	pars[i].GetParValue("sigmaX")<<","<<pars[i].GetParValue("sigmaY")<<") err("<<pars[i].GetParError("sigmaX")<<","<<pars[i].GetParError("sigmaY")<<")"<<endl;
				cout<<"theta="<<pars[i].GetParValue("theta")<<" +- "<<pars[i].GetParError("theta")<<endl;
			}
			cout<<"****************"<<endl;
		}

		/**
		* \brief Set theta fixed
		*/
		void SetThetaFixed(bool choice){thetaFixed=choice;}
		/**
		* \brief Is theta fixed?
		*/
		bool IsThetaFixed(){return thetaFixed;}
		/**
		* \brief Set theta fixed
		*/
		void SetSigmaFixed(bool choice){sigmaFixed=choice;}
		/**
		* \brief Is sigma fixed?
		*/
		bool IsSigmaFixed(){return sigmaFixed;}
		/**
		* \brief Set offset fixed
		*/
		void SetOffsetFixed(bool choice){offsetFixed=choice;}
		/**
		* \brief Is offset fixed?
		*/
		bool IsOffsetFixed(){return offsetFixed;}

		/**
		* \brief Set norm factor
		*/
		void SetNormFactor(double value){normFactor=value;}
		/**
		* \brief Get norm factor
		*/
		double GetNormFactor(){return normFactor;}

		/**
		* \brief Get number of free parameters per component
		*/
		int GetFreeParsPerComponent(){
			int nFreeParsPerComponent= 3;//Ampl + centroids
			if(!sigmaFixed) nFreeParsPerComponent+= 2;//sigmas
			if(!thetaFixed) nFreeParsPerComponent++;
			return nFreeParsPerComponent;
		}

		/**
		* \brief Get component flux derivative matrix
		*/
		int GetComponentFluxDerivMatrix(TMatrixD& D,int componentId){
			//Keep only selected component
			//int nFreeParsPerComponent= GetFreeParsPerComponent();
			//int start_index= componentId*nFreeParsPerComponent;
			//int last_index= start_index + nFreeParsPerComponent-1;
			int start_index= componentId*npars_component;
			int last_index= start_index + npars_component-1;
			int nCols= fluxDensityDerivMatrix.GetNcols(); 
			if(nCols<=last_index){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Trying to access to an not-existing element (index="<<last_index<<") of derivative matrix (dim="<<nCols<<", npars_component="<<npars_component<<") (hint: derivative matrix not initialized)!");	
				#endif
				return -1;
			}
			
			//Init to flux derivative matrix
			D.ResizeTo(1,nCols);
			D.Zero();
			for(int i=start_index;i<=last_index;i++){
				D(0,i)= fluxDensityDerivMatrix(0,i);
			}

			//cout<<"*** DERIV MATRIX COMPONENT "<<componentId<<" ***"<<endl;
			//D.Print();

			return 0;

		}//close GetComponentFluxDerivMatrix()

		/**
		* \brief Set component flag
		*/
		int SetComponentFlag(int componentId,int flag)
		{		
			//Check component id	
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return 0;
			}

			//Get component fit ellipse pars
			pars[componentId].SetFlag(flag);

			return 0;
		}
	
		/**
		* \brief Get component flag
		*/
		int GetComponentFlag(int& flag,int componentId)
		{
			//Check component id
			flag= -1;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}
	
			//Retrieve flag
			flag= pars[componentId].GetFlag();

			return 0;
		}

		/**
		* \brief Set component flag
		*/
		int SetComponentType(int componentId,int type)
		{		
			//Check component id	
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return 0;
			}

			//Get component fit ellipse pars
			pars[componentId].SetType(type);

			return 0;
		}

		/**
		* \brief Get component type
		*/
		int GetComponentType(int& type,int componentId)
		{
			//Check component id
			type= -1;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}
	
			//Retrieve flag
			type= pars[componentId].GetType();

			return 0;
		}

		/**
		* \brief Set component image pix size
		*/
		void SetComponentImagePixSize(double value)
		{		
			for(size_t i=0;i<pars.size();i++){
				pars[i].SetImagePixSize(value);
			}
		}

		/**
		* \brief Get is good fit flag
		*/
		int GetFitQuality(){return fitQuality;}

		/**
		* \brief Set is good fit flag
		*/
		void SetFitQuality(int flag){fitQuality=flag;}

		/**
		* \brief Reset fit pars
		*/
		void Reset()
		{
			
		};

		/**
		* \brief Get fit pars vector
		*/
		int GetFitParVec(std::vector<std::vector<double>>& parVector);

	private:

		/**
		* \brief Initialize component pars
		*/
		void Init();

		
	private:

		int nComponents;
		double chi2;
		double ndof;
		int npars;
		int npars_free;
		int npars_component;
		int nfit_points;
		int status;
		int minimizer_status;
		double offset;
		double offset_err;
		double residualMean;
		double residualRMS;	
		double residualMedian;
		double residualMAD;
		double residualMin;	
		double residualMax;
		std::vector<SourceComponentPars> pars;
		TMatrixD fitCovarianceMatrix;
		TMatrixD fluxDensityDerivMatrix;
		bool thetaFixed;
		bool offsetFixed;
		bool sigmaFixed;
		double fluxDensity;
		double fluxDensityErr;
		int fitQuality;

		double normFactor;

	ClassDef(SourceFitPars,6)

};//close class


#ifdef __MAKECINT__
#pragma link C++ class SourceFitPars+;
#endif

}//close namespace

#endif
