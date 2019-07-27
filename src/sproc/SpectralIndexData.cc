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
* @file SpectralIndexData.cc
* @class SpectralIndexData
* @brief SpectralIndexData class
*
* Class representing spectral index data
* @author S. Riggi
* @date 26/03/2018
*/

#include <SpectralIndexData.h>
#include <Source.h>
#include <CodeUtils.h>


#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TPad.h>
#include <TExec.h>
#include <TStyle.h>
#include <TFitResult.h>

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

using namespace std;

ClassImp(Caesar::SpectralIndexData)

namespace Caesar {


SpectralIndexData::SpectralIndexData() 
	: TObject()
{
	//Initialize
	Init();

}//close costructor


SpectralIndexData::~SpectralIndexData()
{
	//...

}//close destructor

/*
SpectralIndexData::SpectralIndexData(const SpectralIndexData& sid)		
{
	// Copy constructor
  Init();
  ((SpectralIndexData&)sid).Copy(*this);
}


SpectralIndexData& SpectralIndexData::operator=(const SpectralIndexData& sid) { 
	// Operator =
  if (this != &sid) ((SpectralIndexData&)sid).Copy(*this);
  return *this;
}


void SpectralIndexData::Copy(TObject& obj) const 
{
	//Copy object
	TObject::Copy((SpectralIndexData&)obj);

	// Copy variables
  ((SpectralIndexData&)obj).hasSpectralIndex = hasSpectralIndex;
	((SpectralIndexData&)obj).isMultiSourceMatchIndex = isMultiSourceMatchIndex;
	((SpectralIndexData&)obj).spectralIndex = spectralIndex;
	((SpectralIndexData&)obj).spectralIndexErr = spectralIndexErr;
	((SpectralIndexData&)obj).isSpectralIndexFit = isSpectralIndexFit;
	((SpectralIndexData&)obj).spectralFitChi2 = spectralFitChi2;
	((SpectralIndexData&)obj).spectralFitNDF = spectralFitNDF;

}//close Copy()
*/

void SpectralIndexData::Init()
{
	hasSpectralIndex= false;
	spectralIndex= -999;
	spectralIndexErr= 0;
	isMultiSourceMatchIndex= false;
	spectralFitChi2= -999;
	spectralFitNDF= 0;
	isSpectralIndexFit= false;
}

}//close namespace
