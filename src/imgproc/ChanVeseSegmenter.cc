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

#include <ChanVeseSegmenter.h>
#include <Image.h>
#include <Source.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <iostream>

using namespace std;

ClassImp(Caesar::ChanVeseSegmenter)

namespace Caesar {

ChanVeseSegmenter::ChanVeseSegmenter() {

}//close constructor


ChanVeseSegmenter::~ChanVeseSegmenter() {

}//close destructor

void ChanVeseSegmenter::SetCircleLevelSet(TMatrixD* M)
{
	int nRows= M->GetNrows();
	int nCols= M->GetNcols();
	M->Zero();

	int rowCenter= nRows/2;
	int colCenter= nCols/2;
	//double R= std::min(nRows,nCols) * 3.0 / 8.0;
	double R= std::min(nRows,nCols) * 0.1;

	for (int i=0; i<nRows; ++i) {
    for (int j=0; j<nCols; ++j) {
			double x= i - rowCenter;
			double y= j - colCenter;
			double w = R - sqrt(x*x + y*y) ;
			M->operator()(i,j)= w;
    }//end loop cols
  }//end loop rows

}//close SetCircleLevelSet()

void ChanVeseSegmenter::SetCheckerBoardLevelSet(TMatrixD* M, double square_size)
{
	//Init level set
	int nRows= M->GetNrows();
	int nCols= M->GetNcols();
	M->Zero();

	//Fill level set
	for (int i=0; i<nRows; ++i) {
    for (int j=0; j<nCols; ++j) {
			double xx= sin(TMath::Pi()/square_size*i);
			double yy= sin(TMath::Pi()/square_size*j);
			double w= xx*yy;
			M->operator()(i,j)= w;
    }//end loop cols
  }//end loop rows

}//close SetCheckerBoardLevelSet()


ChanVeseSegmenter::CVdata* ChanVeseSegmenter::Init(Image* img,Image* initSegmImg){

	//## Normalize image
	double norm_min= 0;
	double norm_max= 255;
	Image* img_norm= img->GetNormalizedImage("LINEAR",norm_min,norm_max);
	if(!img_norm){
		ERROR_LOG("Failed to normalize input image!");
		return 0;
	}

	CVdata* pCVdata= new CVdata;
	
	//## Convert image to matrix
	(pCVdata->imgMatrix)= img_norm->GetMatrix();
	if(!pCVdata->imgMatrix) {
		ERROR_LOG("Failed to convert input norm image to matrix!");
		if(img_norm) img_norm->Delete();
		if(pCVdata){
			delete pCVdata;
			pCVdata= 0;
		}
		return 0;
	}
	img_norm->Delete();
	int nRows= (pCVdata->imgMatrix)->GetNrows();
	int nCols= (pCVdata->imgMatrix)->GetNcols();

	//## Init images
	pCVdata->phi0 = new TMatrixD(nRows,nCols);
	(pCVdata->phi0)->Zero();
  pCVdata->phi = new TMatrixD(nRows,nCols);
	(pCVdata->phi)->Zero();
  pCVdata->edges = new TMatrixD(nRows,nCols);
	(pCVdata->edges)->Zero();
	pCVdata->imgMatrixOut = new TMatrixD(nRows,nCols);
	(pCVdata->imgMatrixOut)->Zero();
  

	if(initSegmImg){
		INFO_LOG("Initial segmentation given, initializing level set from that...");
		for (int i=0; i<nRows; ++i) {
			long int iy= i;
    	for (int j=0; j<nCols; ++j) {
				long int ix= j;
				double w= initSegmImg->GetPixelValue(ix,iy);
				if(w>0) (pCVdata->phi0)->operator()(i,j)= 1;	
				else (pCVdata->phi0)->operator()(i,j)= -1;//bkg shall be set to negative (not to 0)
    	}//end loop cols
  	}//end loop rows
	}//close if
	else{

		//## Set level set to checkerboard
		INFO_LOG("Initializing level set to circle model...");
		//SetCheckerBoardLevelSet(pCVdata->phi0);
		//SetCircleLevelSet(pCVdata->phi0);

		
		//## Set up initial circular contour for a 256x256 image
		double rowCenter= nRows/2.0;
		double colCenter= nCols/2.0;
		double initContourRadius= 0.5;//1
		INFO_LOG("Initializing level set from dummy gaussian (center("<<rowCenter<<","<<colCenter<<", radiu="<<initContourRadius<<")");
	
		for (int i=0; i<nRows; ++i) {
    	for (int j=0; j<nCols; ++j) {
				double x= double(i) - rowCenter;
				double y= double(j) - colCenter;
				double w= 900.0/(900.0 + x*x + y*y ) - initContourRadius;
				(pCVdata->phi0)->operator()(i,j)= w;
    	}//end loop cols
  	}//end loop rows
		

	}//close else

	return pCVdata;
	
}//close Init()


Image* ChanVeseSegmenter::FindSegmentation(Image* img,Image* initSegmImg,bool returnContourImg,double dt,double h,double lambda1,double lambda2,double mu,double nu,double p,int nIterations)
{
	//Check input image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return nullptr;
	}

	//## Set up parameters
	double norm_min= 0;
	double norm_max= 255;
  struct CVsetup* pCVinputs = new struct CVsetup;
  pCVinputs->dt = dt;
  pCVinputs->h = h;
  pCVinputs->lambda1 = lambda1;
  pCVinputs->lambda2 = lambda2;
  pCVinputs->mu = mu;
  pCVinputs->nu = nu;
  pCVinputs->p = p;
	pCVinputs->niters = nIterations;
	INFO_LOG("Running Chan-Vese segmentation with pars {dt="<<dt<<", h="<<h<<", lambda1="<<lambda1<<", lambda2="<<lambda2<<" mu="<<mu<<" nu="<<nu<<" p="<<p<<", niters="<<nIterations<<"}");
	
	//## Initialize algo data
	ChanVeseSegmenter::CVdata* pCVdata= Init(img,initSegmImg);
	if(!pCVdata){
		ERROR_LOG("Failed to initialize!");
		if(pCVinputs){
			delete pCVinputs;
			pCVinputs= 0;
		}
		return 0;
	}
	int nRows= (pCVdata->imgMatrix)->GetNrows();
	int nCols= (pCVdata->imgMatrix)->GetNcols();

	//## Run segmentation
	INFO_LOG("Performing Chan-Vese segmentation...");
  CVSegmentation(pCVdata->imgMatrix,pCVdata->phi0, &(pCVdata->phi), pCVinputs);

	//## Get zero crossing pixels
  INFO_LOG("Getting zero crossings...");
  ZeroCrossings(pCVdata->phi,&(pCVdata->edges),norm_max,norm_min);
  *(pCVdata->imgMatrixOut)= *(pCVdata->imgMatrix);
  
	//## Fill segmented image (or contour image)
	INFO_LOG("Filling segmented image...");
	TString imgName= Form("%s_CVSegmented",img->GetName().c_str());
	Image* outputImg= img->GetCloned(std::string(imgName),true,true);
	outputImg->Reset();

	if(returnContourImg){
		for (int i=0;i<nRows; ++i) {
			long int iy= i;
			for (int j=0;j<nCols; ++j) {
				long int ix= j;
				double w= (pCVdata->edges)->operator()(i,j);
      	if (w == norm_max) outputImg->FillPixel(ix,iy,1.);
				else outputImg->FillPixel(ix,iy,0.);
    	}//end loop 
  	}//end loop

	}//close if
	else{		
		for (int i=0; i<nRows; ++i) {
			long int iy= i;
			for (int j=0;j<nCols;++j) {
				long int ix= j;
				double w= (pCVdata->phi)->operator()(i,j);	
				if (w>=0) outputImg->FillPixel(ix,iy,norm_max);
				else outputImg->FillPixel(ix,iy,0.);
  		}//end loop
		}//end loop
	}//close else

	//## Clear up
	if(pCVinputs){
		delete pCVinputs;
		pCVinputs= 0;
	}

	if(pCVdata){
		delete pCVdata;
		pCVdata= 0;
	}

	return outputImg;

}//close FindSegmentation()





//---------------------------------------------------------------------------//
// Main segmentation algorithm.  Segment a grayscale image into foreground and
// background regions, given an initial contour defined by the level set function
// phi.  Based on the algorithm described in the paper
// "Active Contours Without Edges" by Chan & Vese.
void ChanVeseSegmenter::CVSegmentation(TMatrixD* img,TMatrixD* phi0,TMatrixD** phi,struct CVsetup* pCVinputs) {
  
	double P_ij;
  double deltaPhi;
  double F1;
  double F2; 
  double F3;
  double F4;
  double F;
  double L;
  double c1;
  double c2;
  
  // Segmentation parameters
  double h = pCVinputs->h;
  double dt = pCVinputs->dt;
  double nu = pCVinputs->nu;
  double lambda1 = pCVinputs->lambda1;
  double lambda2 = pCVinputs->lambda2;
  unsigned int p = pCVinputs->p;
	int nIterations= pCVinputs->niters;
  
  // Variables to evaluate stopping condition
  bool fStop = false;
  TMatrixD* phiOld = new TMatrixD(img->GetNrows(),img->GetNcols());
  
  // Initialize phi
	if (*phi == 0){
    (*phi) = new TMatrixD(img->GetNrows(),img->GetNcols());
	}
 	else if ((*phi)->GetNrows() != phi0->GetNrows() || (*phi)->GetNcols() != phi0->GetNcols() ){
		(*phi)->Delete();
		(*phi) = new TMatrixD(phi0->GetNrows(),phi0->GetNcols());
	}
  
  (**phi) = *phi0;
  
  // Main loop
  for (unsigned int k = 0; k < 5 && fStop == false; ++k) {
    *phiOld= (**phi);
    
    // Compute region averages for current level set function
    // Main segmentation algorithm
    GetRegionAverages(img, *phi, h, c1, c2);

    // Inner loop...
    for (unsigned int l = 0; l < 5; ++l) {

      // Compute length of contour if p > 1
      if (1 == p) L = 1.0;
      else L = 1.0; // fix this!!
     
      // Loop through all interior image pixels
      for (int i = 1; i < img->GetNrows()-1; ++i) {
        for (int j = 1; j < img->GetNcols()-1; ++j) {
          // Compute coefficients needed in update
          GetChanVeseCoefficients(*phi,pCVinputs,i, j,L,F1,F2,F3,F4,F,deltaPhi);

					double phi_ij= (**phi)(i,j);
					double w= (*img)(i,j);
					P_ij = phi_ij - dt*deltaPhi*(nu + lambda1*pow(w-c1,2) - lambda2*pow(w-c2,2));

          // Update level set function
					double phi_iplus1_j= (**phi)(i+1,j);
					double phi_iminus1_j= (**phi)(i-1,j);
					double phi_i_jplus1= (**phi)(i,j+1);
					double phi_i_jminus1= (**phi)(i,j-1);
					(**phi)(i,j) = F1*phi_iplus1_j + F2*phi_iminus1_j + F3*phi_i_jplus1 + F4*phi_i_jminus1 + F*P_ij;
        }
      }
      
      // Update border values of phi by reflection
      for (int i = 0; i < img->GetNrows(); ++i) {
				(**phi)(i,0)= (**phi)(i,1);
				(**phi)(i,img->GetNcols()-1)= (**phi)(i,img->GetNcols()-2);
			}
      for (int j = 0; j < img->GetNcols(); ++j) {
				(**phi)(0,j)= (**phi)(1,j);
				(**phi)(img->GetNrows()-1,j)= (**phi)(img->GetNrows()-2,j);
      }
 
      // Reinitialize phi to the signed distance function to its zero contour
      ReinitPhi(*phi, phi, 0.1, h, nIterations);
    }
    
    // Check stopping condition
    //...
  }
}//close CVSegmentation()


//---------------------------------------------------------------------------//
// function GetRegionAverages
// Compute c1 and c2 as used in the Chan-Vese segmentation algorithm.
// c1 and c2 are given by 
//         c1 = integral(u0*H(phi))dxdy/integral(H(phi))dxdy
//         c2 = integral(u0*(1-H(phi))dxdy/integral(1-H(phi))dxdy
//
// If epsilon == 0, we define H as the usual Heaviside function. Then  c1 is 
// the average of the image pixels over the set where phi is >= 0, and c2 is
// the average over {phi < 0}.  
// If epsilon > 0, we use a smoothed version of the Heaviside function with
// parameter epsilon.
void ChanVeseSegmenter::GetRegionAverages(TMatrixD* img,TMatrixD* phi,double epsilon, double &c1, double &c2) {
  
	// Non-smoothed calculation
  if (0 == epsilon) {
    int n1 = 0;
    int n2 = 0;
    double Iplus = 0;
    double Iminus = 0;

    for (int i=0; i<img->GetNrows(); ++i) {
			for (int j=0; j<img->GetNcols();++j) {
				double w_phi= (*phi)(i,j);
				double w= (*img)(i,j);
        if (w_phi>= 0){
          ++n1;
          Iplus += w;
        }
        else {
          ++n2;
          Iminus += w;
        }
      }//end loop bins X
    }//end loop bins Y
    c1 = Iplus/double(n1);
    c2 = Iminus/double(n2);
  }//close if non smoothed
  // Smoothed calculation
  else {
    double num1 = 0;
    double den1 = 0;
    double num2 = 0;
    double den2 = 0;
    double H_phi;
    for (int i = 0; i < phi->GetNrows(); ++i){
      for (int j = 0; j < phi->GetNcols(); ++j) {
        // Compute value of H_eps(phi) where H_eps is a mollified Heavyside function
				double w_phi= (*phi)(i,j);
				double w= (*img)(i,j);
        H_phi = .5*(1+(2/TMath::Pi())*atan(w_phi/epsilon));
        num1 += w*H_phi;
        den1 += H_phi;
        num2 += w*(1-H_phi);
        den2 += 1-H_phi;
      }//end loop bins Y
    }//end loop bins X
    
    c1 = num1/den1;
    c2 = num2/den2;
  }//close else

}//close ChanVeseSegmenter::GetRegionAverages()

//---------------------------------------------------------------------------//
// function ReinitPhi
// Reinitialize phi to the signed distance function to its zero contour
void ChanVeseSegmenter::ReinitPhi(TMatrixD* phiIn,TMatrixD** psiOut,double dt,double h,unsigned int numIts) {

  if (*psiOut == 0){
    (*psiOut) = new TMatrixD(phiIn->GetNrows(),phiIn->GetNcols());
	}
  else if ((*psiOut)->GetNrows() != phiIn->GetNrows() || (*psiOut)->GetNcols() != phiIn->GetNcols() ){
		(*psiOut)->Delete();
		(*psiOut) = new TMatrixD(phiIn->GetNrows(),phiIn->GetNcols());
  }
  
  (**psiOut) = *phiIn;


  double a;
  double b;
  double c;
  double d;
  double x;
  double G;
  
  bool fStop = false;
  double Q;
  unsigned int M;
  TMatrixD* psiOld = new TMatrixD(phiIn->GetNrows(),phiIn->GetNcols());
  
  for (unsigned int k = 0; k < numIts && fStop == false; ++k) {
    *psiOld= (**psiOut);
    for (int i = 1; i < phiIn->GetNrows()-1; ++i) {
      for (int j = 1; j < phiIn->GetNcols()-1; ++j){
				double phiIn_ij= (*phiIn)(i,j);
				double phiIn_iminus1_j= (*phiIn)(i-1,j);
				double phiIn_iplus1_j= (*phiIn)(i+1,j);
				double phiIn_i_jminus1= (*phiIn)(i,j-1);
				double phiIn_i_jplus1= (*phiIn)(i,j+1);

        a = (phiIn_ij - phiIn_iminus1_j)/h;
        b = (phiIn_iplus1_j - phiIn_ij)/h;
        c = (phiIn_ij - phiIn_i_jminus1)/h;
        d = (phiIn_i_jplus1 - phiIn_ij)/h;
        
        if (phiIn_ij > 0)
          G = sqrt(max(max(a,0.0)*max(a,0.0),min(b,0.0)*min(b,0.0))
                 + max(max(c,0.0)*max(c,0.0),min(d,0.0)*min(d,0.0))) - 1.0;
        else if (phiIn_ij < 0)
          G = sqrt(max(min(a,0.0)*min(a,0.0),max(b,0.0)*max(b,0.0))
                 + max(min(c,0.0)*min(c,0.0),max(d,0.0)*max(d,0.0))) - 1.0;
        else
          G = 0;
        
        x = (phiIn_ij >= 0)?(1.0):(-1.0);
				(**psiOut)(i,j)= (**psiOut)(i,j) - dt*x*G;
      }
    }
    
    // Check stopping condition
    Q = 0.0;
    M = 0.0;
    for (int i = 0; i < phiIn->GetNrows(); ++i) {
      for (int j = 0; j < phiIn->GetNcols(); ++j) {
				double w_old= (*psiOld)(i,j);	
				double w_phi= (**psiOut)(i,j);
        if (abs(w_old) <= h) {
        	++M;
          Q += abs(w_old - w_phi);
        }
      }//end loop cols
    }//end loop rows
    if (M != 0)
      Q = Q/((double)M);
    else
      Q = 0.0;
    
    if (Q < dt*h*h){
      fStop = true;
      INFO_LOG("Stopping condition reached at " << k+1 << " iterations (Q="<<Q<<")");
    }
    else {
      //cout << "Iteration " << k << ", Q = " << Q << " > " << dt*h*h << endl;
    }
  }//end iteration loop

}//close ChanVeseSegmenter::ReinitPhi()

//---------------------------------------------------------------------------//
// function GetChanVeseCoefficients
// Compute coefficients needed in Chan-Vese segmentation algorithm given current 
// level set function
void ChanVeseSegmenter::GetChanVeseCoefficients(TMatrixD* phi,
                             struct CVsetup* pCVinputs,
                             unsigned int i,
                             unsigned int j,
                             double L,
                             double& F1,
                             double& F2,
                             double& F3,
                             double& F4,
                             double& F,
                             double& deltaPhi)
{
  // factor to avoid division by zero
  double eps = 0.000001;
  double h = pCVinputs->h;
  double dt = pCVinputs->dt;
  double mu = pCVinputs->mu;
  unsigned int p = pCVinputs->p;
  
	double phi_ij= (*phi)(i,j);	
	double phi_iminus1_j= (*phi)(i-1,j);
	double phi_iplus1_j= (*phi)(i+1,j);
	double phi_i_jminus1= (*phi)(i,j-1);
	double phi_i_jplus1= (*phi)(i,j+1);
	double phi_iminus1_jminus1= (*phi)(i-1,j-1);
	double phi_iplus1_jminus1= (*phi)(i+1,j-1);
	double phi_iminus1_jplus1= (*phi)(i-1,j+1);
	

  double C1 = 1/sqrt(eps + pow((phi_iplus1_j - phi_ij),2)
                         + pow((phi_i_jplus1 - phi_i_jminus1),2)/4.0);
  double C2 = 1/sqrt(eps + pow((phi_ij - phi_iminus1_j),2)
                         + pow((phi_iminus1_jplus1 - phi_iminus1_jminus1),2)/4.0);
  double C3 = 1/sqrt(eps + pow((phi_iplus1_j - phi_iminus1_j),2)/4.0
                         + pow((phi_i_jplus1- phi_ij),2));
  double C4 = 1/sqrt(eps + pow((phi_iplus1_jminus1 - phi_iminus1_jminus1),2)/4.0
                         + pow((phi_ij - phi_i_jminus1),2));

  deltaPhi = h/(TMath::Pi()*(h*h + (phi_ij)*(phi_ij)));
  
  double Factor = dt*deltaPhi*mu*(double(p)*pow(L,p-1));
  F = h/(h+Factor*(C1+C2+C3+C4));
  Factor = Factor/(h+Factor*(C1+C2+C3+C4));
  
  F1 = Factor*C1;
  F2 = Factor*C2;
  F3 = Factor*C3;
  F4 = Factor*C4;

}//close GetChanVeseCoefficients()

void ChanVeseSegmenter::ZeroCrossings(TMatrixD* imageIn,TMatrixD** edges,double fg,double bg) {

  // Allocate output image if necessary
  if (0 == (*edges)){
		(*edges) = new TMatrixD(imageIn->GetNrows(),imageIn->GetNcols());
  }
  else {
    if ((*edges)->GetNrows() != imageIn->GetNrows() || (*edges)->GetNcols() != imageIn->GetNcols()) {
			(*edges)->Delete();
			(*edges) = new TMatrixD(imageIn->GetNrows(),imageIn->GetNcols());
    } 
  }
  (**edges)= bg;

  for (int i = 0; i < imageIn->GetNrows(); ++i) {
    for (int j = 0; j < imageIn->GetNcols(); ++j) {
      // Currently only checking interior pixels to avoid bounds checking
      if (i > 0 && i < (imageIn->GetNrows()-1) && j > 0 && j < (imageIn->GetNcols()-1)) {

				double w_ij= (*imageIn)(i,j);
				double w_iminus1_jminus1= (*imageIn)(i-1,j-1);
				double w_iplus1_jplus1= (*imageIn)(i+1,j+1);
				double w_iminus1_j= (*imageIn)(i-1,j);
				double w_iplus1_j= (*imageIn)(i+1,j);
				double w_iplus1_jminus1= (*imageIn)(i+1,j-1);
				double w_iminus1_jplus1= (*imageIn)(i-1,j+1);
				double w_i_jminus1= (*imageIn)(i,j-1);
				double w_i_jplus1= (*imageIn)(i,j+1);

        if (0 == w_ij) {
          if (0 != w_iminus1_jminus1
           || 0 != w_iminus1_j
           || 0 != w_iminus1_jplus1
           || 0 != w_i_jminus1
           || 0 != w_i_jplus1
           || 0 != w_iplus1_jminus1
           || 0 != w_iplus1_j
           || 0 != w_iplus1_jplus1)
           {
             (**edges)(i,j) = fg;
           }
        }//close if
        else {
          if (abs(w_ij) < abs(w_iminus1_jminus1) && (w_ij>0) != (w_iminus1_jminus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iminus1_j) && (w_ij>0) != (w_iminus1_j>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iminus1_jplus1) && (w_ij>0) != (w_iminus1_jplus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_i_jminus1) && (w_ij>0) != (w_i_jminus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_i_jplus1) && (w_ij>0) != (w_i_jplus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iplus1_jminus1) && (w_ij>0) != (w_iplus1_jminus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iplus1_j) && (w_ij>0) != (w_iplus1_j>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iplus1_jplus1) && (w_ij>0) != (w_iplus1_jplus1>0))
             (**edges)(i,j) = fg;
        }
      }
    }
  }
}//close ZeroCrossing()


}//close namespace

