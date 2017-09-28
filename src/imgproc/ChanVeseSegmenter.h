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
* @file ChanVeseSegmenter.h
* @class ChanVeseSegmenter
* @brief Class implementing ChanVese segmentation algorithm
*
* @author R. Crandall, S. Riggi
* @date 15/06/2015
*/

#ifndef _CHANVESE_SEGMENTER_H
#define _CHANVESE_SEGMENTER_H 1


#include <Contour.h>

#include <TObject.h>
#include <TMatrixD.h>

#include <stdio.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <float.h>


namespace Caesar{

class Image;

class ChanVeseSegmenter : public TObject {

	public:
  	/**
		\brief Class constructor
		*/
    ChanVeseSegmenter();
		/**
		\brief Class destructor
		*/
    ~ChanVeseSegmenter();
		
		// Define structure containing parameters of
		struct CVsetup {
  		double dt; // time step 
  		double h;  // pixel spacing
  		double lambda1;
  		double lambda2;
  		double mu; // contour length weighting parameter
  		double nu; // region area weighting parameter
  		unsigned int p; // length weight exponent
			int niters;
		};

		struct CVdata {
			CVdata(){
				imgMatrix= 0;
				phi0= 0;		
				phi= 0;
				edges= 0;
				imgMatrixOut= 0;
			}
			~CVdata(){
				if(imgMatrix) imgMatrix->Delete();
				if(phi0) phi0->Delete();
				if(phi) phi->Delete();
				if(edges) edges->Delete();
				if(imgMatrixOut) imgMatrixOut->Delete();
			}
			TMatrixD* imgMatrix;
			TMatrixD* phi0;
			TMatrixD* phi;
			TMatrixD* edges;
			TMatrixD* imgMatrixOut;
		};


	public:
		
		/**
		* \brief Find the ChanVese segmentation of input image
		*/
		static Image* FindSegmentation(Image* img,Image* initSegmImg=0,bool returnContourImg=false,double dt=0.1,double h=1,double lambda1=1.0,double lambda2=2.0,double mu=0.5,double nu=0,double p=1,int nIterations=1000);

		

	private:

		/**
		* \brief Initialize algorithm
		*/
		static CVdata* Init(Image* img,Image* initSegmImg=0);


		//===============================================
		//==          CHAN-VESE INTERNAL METHODS
		//===============================================
		/**
		* \brief Main segmentation algorithm
		*/
		static void CVSegmentation(TMatrixD* img,TMatrixD* phi0,TMatrixD** phi,struct CVsetup* pCVinputs);
    
		/**
		* \brief Compute gray level averages in foreground and background regions defined by level set function phi
		*/
		static void GetRegionAverages(TMatrixD* img, TMatrixD* phi,double epsilon,double &c1,double &c2);

		/**
		* \brief Compute coefficients needed in Chan-Vese segmentation algorithm given current level set function
		*/
		static void GetChanVeseCoefficients(TMatrixD* phi,struct CVsetup* pCVinputs,
														 unsigned int i,
														 unsigned int j,
                             double L,
                             double& F1,
                             double& F2,
                             double& F3,
                             double& F4,
                             double& F,
                             double& deltaPhi);
                             
    /**
		* \brief Reinitialize a function to the signed distance function to its zero contour
		*/                         
		static void ReinitPhi(TMatrixD* phiIn,TMatrixD** psiOut,double dt,double h,unsigned int numIts);
		/**
		* \brief Compute zero crossings
		*/ 
		static void ZeroCrossings(TMatrixD* imageIn,TMatrixD** edges,double fg,double bg);

	
	private:

		ClassDef(ChanVeseSegmenter,1)

};//close class
		 

#ifdef __MAKECINT__
#pragma link C++ class ChanVeseSegmenter+;
#endif

}//close namespace 
                    
#endif


