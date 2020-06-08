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
* @file StatsUtils.cc
* @class StatsUtils
* @brief Utility functions for statistical tasks
*
* Utility functions for statistical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <StatsUtils.h>

#include <TObject.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::StatsUtils)

ClassImp(Caesar::ClippedStats<int>)
ClassImp(Caesar::ClippedStats<long int>)
ClassImp(Caesar::ClippedStats<double>)
ClassImp(Caesar::ClippedStats<float>)

ClassImp(Caesar::StatMoments<int>)
ClassImp(Caesar::StatMoments<long int>)
ClassImp(Caesar::StatMoments<double>)
ClassImp(Caesar::StatMoments<float>)

ClassImp(Caesar::BoxStats<int>)
ClassImp(Caesar::BoxStats<long int>)
ClassImp(Caesar::BoxStats<double>)
ClassImp(Caesar::BoxStats<float>)

namespace Caesar {

StatsUtils::StatsUtils(){

}

StatsUtils::~StatsUtils(){

}

double StatsUtils::GetMahalanobisDistance(TMatrixD x, TMatrixD mean, TMatrixD Sigma,bool isInverted){
	double dist= 0;
	if(!isInverted) Sigma.Invert();
	TMatrixD diff= x-mean;
	TMatrixD diffT(TMatrixD::kTransposed,diff);
	TMatrixD M= diffT*Sigma*diff;
	dist= sqrt(M(0,0));
	return dist;
}//close GetMahalanobisDistance()

int StatsUtils::ComputePageRank(std::vector<double>& ranks,TMatrixD& M,double d,double tol){

	//M is the Google rank adjacency matrix
	int N= M.GetNcols();
	if(N<=0) return -1;
	
	ranks.clear();
	bool hasConverged= true;
	int maxIterToStop= 1000;
	
	//## Initialize rank with uniform prob (0.5)
	TVectorD v(N);
	TVectorD v_last(N);
	TMatrixD UnoMatrix(N,N);	
	UnoMatrix.Zero();
	for(int i=0;i<N;i++) {
		v(i)= 0.5;
		v_last(i)= 1.e+99;
		for(int j=0;j<N;j++) UnoMatrix(i,j)= 1;
	}
	
	double norm2= v.Norm2Sqr();
	v*= 1./sqrt(norm2);
	
	//## Compute Mhat= (d x M) + (1-d)/N*UnoMatrix(N,N)
	TMatrixD M_hat(N,N);
	M_hat.Zero();
	M_hat= d*M + ((1.-d)/(double)(N))*UnoMatrix;
	
	int iter= 0;
	int nIterToConverge= 0;
	double diff= (v-v_last).Norm2Sqr();
	
	while( diff>tol ){
		v_last= v;
    v = M_hat*v;
		double normFactor= sqrt(v.Norm2Sqr());
    v*= 1./normFactor;
		
		diff= sqrt((v-v_last).Norm2Sqr());	
		iter++;
		nIterToConverge= iter;

    if(iter >= maxIterToStop){
    	hasConverged= false;
			break;
		}
	}//end loop 
	
	if(hasConverged){
		for(int k=0;k<v.GetNoElements();k++) {
			ranks.push_back(v[k]);
		}
	}
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Page rank did not converge!");
		#endif
		return -1;
	}

	return 0;

}//close ComputePageRank()

int StatsUtils::GenerateFitParsAroundCovMatrix(std::vector<std::vector<double>>& fitPars_rand,const std::vector<double>& fitPars,const std::vector<std::vector<double>>& fitCovMatrix,int nsamples)
{
	//## Check input data
	//Clear data 
	fitPars_rand.clear();

	//Check input cov matrix
	if(fitCovMatrix.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cov matrix vector is empty!");
		#endif
		return -1;
	}

	int nRows= static_cast<int>(fitCovMatrix.size());
	int nCols= static_cast<int>(fitCovMatrix[0].size());
	for(int i=1;i<nRows;i++){
		int s= static_cast<int>(fitCovMatrix[i].size());
		if(s!=nCols){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Cov matrix vector is not a square matrix!");
			#endif
			return -1;
		}
	}

	//Create TMatrix	
	TMatrixD C(nRows,nCols);
	for(int i=0;i<nRows;i++){
		for(int j=0;j<nCols;j++){
			C(i,j)= fitCovMatrix[i][j];
		}
	}

	//Generate random fit pars
	return GenerateFitParsAroundCovMatrix(fitPars_rand,fitPars,C,nsamples);

}//close GenerateFitParsAroundCovMatrix()


int StatsUtils::GenerateFitParsAroundCovMatrix(std::vector<std::vector<double>>& fitPars_rand,const std::vector<double>& fitPars,TMatrixD& fitCovMatrix,int nsamples)
{
	//## Check inputs
	//Clear data
	fitPars_rand.clear();

	//Check if fit pars are given
	if(fitPars.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No pars given in input!");
		#endif
		return -1;
	}
	//Check cov matrix size
	int nRows= fitCovMatrix.GetNrows();
	int nCols= fitCovMatrix.GetNcols();
	if(nRows!=nCols || nRows<=0 || nCols<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given cov matrix is not a square matrix ("<<nRows<<"!="<<nCols<<") or matrix is not initialized!");
		#endif
		return -1;
	}
	int nData= static_cast<int>(fitPars.size());
	if(nRows!=nData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given cov matrix size is different from par size ("<<nData<<")!");
		#endif
		return -1;
	}

	
	//## Find a Cholesky decomposition of the covariance matrix for error propagation
	TDecompChol CCholDecomp(fitCovMatrix);
  CCholDecomp.Decompose();
	
	TMatrixD* CCholDecompTriang= new TMatrixD(nRows,nCols);
	CCholDecompTriang= (TMatrixD*)(&CCholDecomp.GetU());
	
	TMatrixD* CCholDecompTriangTransp= new TMatrixD(nRows,nCols);
	CCholDecompTriangTransp->Transpose(*CCholDecompTriang);


	//Generate rand pars
	for(int k=0;k<nsamples;k++)
	{
		//Initialize rand fit pars
		fitPars_rand.push_back( std::vector<double>(nRows,0) );

		//Randomize pars
		TMatrixD R(nRows,1);
		TMatrixD V(nRows,1);
		for(int i=0;i<nRows;i++) R(i,0)= gRandom->Gaus(0,1);		
		V= (*CCholDecompTriangTransp)*R;

		for(int i=0;i<nRows;i++) {
			double parValue= fitPars[i];
			double parValue_rand= parValue + V(i,0); 
			fitPars_rand[k][i]= parValue_rand;
		}

	}//end loop samples

	//## Delete data
	if(CCholDecompTriang) CCholDecompTriang->Delete();
	if(CCholDecompTriangTransp) CCholDecompTriangTransp->Delete();

	return 0;

}//close GenerateRandFitPars()

}//close namespace



