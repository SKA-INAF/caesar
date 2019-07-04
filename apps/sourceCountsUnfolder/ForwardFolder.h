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
* @file ForwardFolder.h
* @class ForwardFolder
* @brief Forward folding of a count spectrum
* 
* @author S. Riggi
* @date 03/07/2019
*/
#ifndef _FORWARD_FOLDER_H_
#define _FORWARD_FOLDER_H_

#include <SpectrumUtils.h>

#include <TObject.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TMatrixD.h>
#include <TTree.h>


class ForwardFolder 
{
  public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
    ForwardFolder();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~ForwardFolder();
	
	public:
		/**
		* \brief Unfold a spectrum
		*/
		int RunUnfold(TH1D* spectrumHisto,TH2D* responseMatrix,SpectrumPars& initFitPars,bool computeUncertainties=true,int nRandomSamples=100,bool useFitRange=false,double fitMin=-1,double fitMax=-1);

		/**
		* \brief Return unfolded spectrum
		*/
		TH1D* GetUnfoldedSpectrum(){return fUnfoldedSpectrum;}
		/**
		* \brief Return true spectrum
		*/
		TH1D* GetTrueSpectrum(){return fTrueSpectrum;}
		/**
		* \brief Return rec spectrum
		*/
		TH1D* GetRecSpectrum(){return fRecSpectrum;}
		/**
		* \brief Return forward folded spectrum
		*/
		TH1D* GetForwardFoldedSpectrum(){return fForwardFoldedSpectrum;}

	private:
		/**
		* \brief Initialize class data
		*/
		int Init();
		/**
		* \brief Clear class data
		*/
		int Clear();
		/**
		* \brief Unfold spectrum (internal method)
		*/
		TH1D* UnfoldSpectrum(SpectrumPars& initFitPars,std::string runMode="FIT");
		/**
		* \brief Fit log-likelihood definition
		*/
		static void LogLikelihoodFcn(int& nPar, double* gin, double &f, double* par, int iflag);		
		/**
		* \brief Compute uncertainties on unfolded spectrum
		*/
		int ComputeUncertainties(SpectrumPars& initFitPars,int nRandomSamples=100);

	private:
		
		//Data
		static TF1* fTrueSpectrumModelFcn;
		TH1D* fRecSpectrum;
		static TH1D* fCurrentRecSpectrum;		
		static TH1D* fTrueSpectrum;
		static TH1D* fForwardFoldedSpectrum;
		static TH1D* fCurrentForwardFoldedSpectrum;
		TH1D* fUnfoldedSpectrum;
		static TH2D* fResponseMatrix;
		int fNTrueBins;
		double fLgEMin_true;
		double fLgEMax_true;
		std::vector<double> fLgEBins_true;

		int fNRecBins;
		double fLgEMin_rec;
		double fLgEMax_rec;
		std::vector<double> fLgEBins_rec;
		
		//Options
		static bool fUseFitRange;
		static double fLgEMin_fit;
		static double fLgEMax_fit;
		
		//Fit results
		int fFitStatus;
		static const int NMAXFITTEDPARS= 10;
		TMatrixD* fCovarianceMatrix;	
		double fFitPar[NMAXFITTEDPARS];
		double fFitParErr[NMAXFITTEDPARS];

		static SpectrumPars* fInitFitPars;

};//close class


#endif

