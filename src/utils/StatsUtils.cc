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
		cerr<<"Utils::ComputePageRank(): WARN: Page rank did not converge!"<<endl;
		return -1;
	}

	return 0;

}//close ComputePageRank()

}//close namespace



