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
* @file SourceFitter.h
* @class SourceFitter
* @brief SourceFitter
*
* Class to fit a source image with a mixture of gaussian functions
* @author S. Riggi
* @date 01/09/2017
*/

#ifndef _SOURCE_FITTER_h
#define _SOURCE_FITTER_h 1

#include <SourceFitPars.h>

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

class TH2D;

namespace Caesar {

class Image;
class Source;

struct SourceFitOptions {

	public:

	
	//Default constructor
	SourceFitOptions() 
	{
    bmaj= 5;//pix
    bmin= 5;//pix
		bpa= 0;
		nMaxComponents= 3;	
		fixCentroidInPreFit= false;
		limitCentroidInFit= true;
		centroidLimit= 5;//pixels
		fixBkg= true;
		limitBkgInFit= true;
		useEstimatedBkgLevel= true;
		useBkgBoxEstimate= false;
		fixedBkgLevel= 0;
		fixAmplInPreFit= false;
		limitAmplInFit= true;
		amplLimit= 0.2;
		limitSigmaInFit= true;
		sigmaLimit= 0.2;
		fixSigmaInPreFit= true;
		fixSigma= false;
		fixThetaInPreFit= true;
		fixTheta= false;
		limitThetaInFit= true;
		thetaLimit= 5;//deg
		useFluxZCut= false;
		fluxZThrMin= 2.5;//in nsigmas
		peakZThrMin= 0;//in nsigmas
		peakMinKernelSize= 3;
		peakMaxKernelSize= 7;
		peakKernelMultiplicityThr= 1;
		peakShiftTolerance= 2;
		setBinErrorsToMapRMS= true;
		
		scaleMin= 3;
		scaleMax= 3;
		scaleStep= 1;
		minBlobSize= 5;
		blobMapThrFactor= 0;
		blobMapKernelFactor= 6; 
		useNestedAsComponents= false;
	
		fitImproveConvergence= true;
		fitNRetries= 1000;
		fitDoFinalMinimizerStep= true;
		fitFinalMinimizer= eHESS;
		chi2RegPar= 0.1;
		fitRetryWithLessComponents= true;

		//Supported minimizers & algo
		//  - Minuit: {Migrad,Simplex,Scan,Minimize}
		//  - Minuit2: {Migrad,Simplex,Combined,Scan,Fumili}
		//  - GSLMultiMin: {conjugatefr,conjugatepr,bfgs,bfgs2,steepestdescent} NB: No error calculation supported so do not use it
		//  - GSLMultiFit: {} NB: Requiring a different chi2 function definition, so do not use it
		fitMinimizer= "Minuit2";
		fitMinimizerAlgo= "Migrad";
		fitPrintLevel= 1;
		fitStrategy= 2;
		fitFcnTolerance= 1.e-2;//default tolerance used in ROOT
		fitMaxIters= 100000;
		fitParBoundIncreaseStepSize= 0.1;

		wcsType= eJ2000;
		fitScaleDataToMax= false;

		//Selection cuts
		useRedChi2Cut= true;
		fitRedChi2Cut= 5.;
		useFitEllipseCuts= false;
		fitEllipseEccentricityRatioMinCut= 0.5;
		fitEllipseEccentricityRatioMaxCut= 1.5;
		fitEllipseAreaRatioMinCut= 0.01;
		fitEllipseAreaRatioMaxCut= 10;
		fitEllipseRotAngleCut= 45;
		
		useSimpleWCSEllipseConversion= true;

	}//close constructor

	public:
		//- Blob start fit pars
		double bmaj;//in pixels
		double bmin;//in pixels
		double bpa;//in degrees

		//- Max number of components to be fitted & options
		int nMaxComponents;
		bool useNestedAsComponents;
	
		//- Number of matching peaks across kernels (-1=peak detected in all kernels, ...) and distance tolerance
		int peakKernelMultiplicityThr;
		int peakShiftTolerance;

		//- Dilation kernels to be used when finding peaks
		int peakMinKernelSize;
		int peakMaxKernelSize;

		//- Multiscale blob finder options
		double scaleMin;
		double scaleMax;
		double scaleStep;
		int minBlobSize;
		double blobMapThrFactor;
		int blobMapKernelFactor;

		//- Peak flux significance min threshold (in nsigmas wrt to avg bkg & rms)
		double peakZThrMin;	

		//- Centroid options
		bool limitCentroidInFit;
		double centroidLimit;//in pixels
		bool fixCentroidInPreFit;

		//- Bkg options
		bool fixBkg;
		bool limitBkgInFit;
		bool useEstimatedBkgLevel;
		bool useBkgBoxEstimate;
		double fixedBkgLevel;
		
		//- Amplitude fit par range (example +-20% around source peak)
		bool limitAmplInFit;
		double amplLimit;
		bool fixAmplInPreFit;

		//- Sigma fit par range
		bool fixSigmaInPreFit;
		bool fixSigma;	
		bool limitSigmaInFit;
		double sigmaLimit;

		//- Theta 
		bool fixThetaInPreFit;
		bool fixTheta;
		bool limitThetaInFit;
		double thetaLimit;//in deg

		//- Flux significance min threshold (in nsigmas above the bkg)
		bool useFluxZCut;
		double fluxZThrMin;

		//- Fit data bin error (if true set bin errors to pixel noise rms, otherwise to 1)
		bool setBinErrorsToMapRMS;

		//- Fit minimization options
		double fitFcnTolerance;
		long int fitMaxIters;
		bool fitImproveConvergence;
		long int fitNRetries;
		bool fitDoFinalMinimizerStep;
		int fitFinalMinimizer;
		double chi2RegPar;
		bool fitRetryWithLessComponents;

		//- Fit minimizer
		int fitPrintLevel;
		int fitStrategy;
		std::string fitMinimizer;
		std::string fitMinimizerAlgo;
		double fitParBoundIncreaseStepSize;		

		//- Fit ellipse pars
		int wcsType;
	
		//- Scale data to max
		bool fitScaleDataToMax;

		//- Selection cuts
		bool useRedChi2Cut;
		double fitRedChi2Cut;
		bool useFitEllipseCuts;
		double fitEllipseEccentricityRatioMinCut;
		double fitEllipseEccentricityRatioMaxCut;
		double fitEllipseAreaRatioMinCut;
		double fitEllipseAreaRatioMaxCut;
		double fitEllipseRotAngleCut;

		//- Use simple WCS ellipse conversion
		bool useSimpleWCSEllipseConversion;
		
};//close SourceFitOptions





//========================================
//==         SOURCE FITTER
//========================================

class SourceFitter : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFitter();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceFitter();

		/**
		* \brief Fit status enum flag
		*/
		/*
		enum FitStatusFlag {
			eFitUnknownStatus= 0,
			eFitAborted= 1,
			eFitNotConverged= 2,
			eFitConverged= 3,
			eFitConvergedWithWarns= 4
		};
		*/

		/**
		* \brief Source fit data
		*/
		struct SourceFitData {
			SourceFitData(){
				x= 0;
				y= 0;
				S= 0;
				Serr= 0;	
			}

			SourceFitData(double _x,double _y,double _S,double _Serr)
				: x(_x),y(_y),S(_S),Serr(_Serr)
			{}
			double x;
			double y;
			double S;
			double Serr;
		};


	public:

		/**
		* \brief Fit source
		*/
		int FitSource(Source* source,SourceFitOptions& fitOptions);
		
		/**
		* \brief Get fit pars
		*/
		SourceFitPars GetFitPars(){return m_sourceFitPars;}

		/**
		* \brief Get fit status
		*/
		int GetFitStatus(){return m_fitStatus;}

		/**
		* \brief 2D Gaussian mixture model used for the fit
		*/
		//static double Gaus2DMixtureFcn(double* x, double* p);
		double Gaus2DMixtureFcn(double* x, const double* p);

		/**
		* \brief 2D Gaussian model used for the fit
		*/
		//static double Gaus2DFcn(double* x, double* p);
		double Gaus2DFcn(double* x, const double* p);

	protected:
		/**
		* \brief Check if fit has parameters converged at bounds
		*/
		bool HasFitParsAtLimit(const ROOT::Fit::FitResult& fitRes);
		/**
		* \brief Returns list of parameters at limits (version with Minimizer)
		*/
		int GetParsAtLimits(std::vector<int>& parsAtLimits,ROOT::Math::Minimizer* fitter);
		
		/**
		* \brief Returns list of parameters at limits (version with TFitter)
		*/
		//int GetParsAtLimits(std::vector<int>& parsAtLimits,TFitter* fitter);
		
		/**
		* \brief Estimate fit components
		*/
		int EstimateFitComponents(std::vector<std::vector<double>>& fitPars_start,Source* aSource,SourceFitOptions& fitOptions);
		
		/**
		* \brief Perform fit using given initial fit components pars
		*/
		int DoChi2Fit(Source* aSource,SourceFitOptions& fitOptions,std::vector< std::vector<double> >& fitPars_start);

		/**
		* \brief Init data
		*/
		int InitData(Source* aSource,SourceFitOptions& fitOptions);
		/**
		* \brief Check fit options
		*/
		int CheckFitOptions(SourceFitOptions& fitOptions);

		/**
		* \brief Chi2 fit function 
		*/
		double Chi2Fcn(const double* par);
			
		/**
		* \brief Chi2 fit function (old MINUIT version)
		*/
		//static void Chi2Fcn(int& nPar,double* grad,double& fval,double* p,int iflag);
		

		/**
		* \brief Initialize fitter
		*/
		ROOT::Math::Minimizer* InitMinimizer(int nFitPars,SourceFitOptions& fitOptions);
	
	private:
	
		//static int m_NFitComponents;
		//static int m_fitStatus;
		//static SourceFitPars m_sourceFitPars;
		//static TH2D* m_fluxMapHisto;
		//static std::vector<SourceFitData> m_fitData;
		//static std::vector<SourceFitData> m_fitHaloData;
		//static double m_chi2RegPar;
		//static double m_bkgMean;
		//static double m_rmsMean;
		//static double m_sourceX0;
		//static double m_sourceY0;

		int m_NFitComponents;
		int m_fitStatus;
		SourceFitPars m_sourceFitPars;
		
		//TH2D* m_fluxMapHisto;
		std::vector<SourceFitData> m_fitData;
		std::vector<SourceFitData> m_fitHaloData;
		double m_chi2RegPar;
		double m_bkgMean;
		double m_rmsMean;
		double m_sourceX0;
		double m_sourceY0;
		
	ClassDef(SourceFitter,1)

};//close SourceFitter class

#ifdef __MAKECINT__
#pragma link C++ class SourceComponentPars+;
#pragma link C++ class SourceFitPars+;
#pragma link C++ class SourceFitter+;
#pragma link C++ struct SourceFitter::SourceFitData+;
#pragma link C++ enum FitStatusFlag;
#pragma link C++ struct SourceFitOptions+;
#endif

}//close namespace

#endif
