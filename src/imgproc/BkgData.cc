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
* @file BkgData.cc
* @class BkgData
* @brief BkgData
*
* Class for storing bkg data
* @author S. Riggi
* @date 20/01/2015
*/

#include <BkgData.h>
#include <Image.h>

#include <TObject.h>
#include <TString.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
using namespace std;


ClassImp(Caesar::BkgSampleData)
ClassImp(Caesar::ImgBkgData)

namespace Caesar {


ImgBkgData::ImgBkgData() 
{
	//Initialize data
	Init();

}//close costructor

ImgBkgData::~ImgBkgData() 
{
	//Clear allocated data
	Clear();

}//close destructor


ImgBkgData::ImgBkgData(const ImgBkgData& data) {
  // Copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((ImgBkgData&)data).Copy(*this);
}

void ImgBkgData::Init(){
	
	BkgSamplings.clear();
	BkgMap= 0;
	NoiseMap= 0;
	gBkg= 0;
	gNoise= 0;

}//close Init()

void ImgBkgData::Copy(TObject &obj) const 
{	
	// Copy variables
  ((ImgBkgData&)obj).gBkg = gBkg;
	((ImgBkgData&)obj).gNoise = gNoise;
				
	//Copy sample bkg data
	//Delete first any existing collection
	((ImgBkgData&)obj).ClearSamplings();

	for(size_t i=0;i<BkgSamplings.size();i++){
		(((ImgBkgData&)obj).BkgSamplings).push_back(BkgSamplings[i]);
	}

	
	//Delete existing maps first
	DEBUG_LOG("Deleting existing bkg & noise maps...");
	if( ((ImgBkgData&)obj).BkgMap ){
		delete ((ImgBkgData&)obj).BkgMap;
		((ImgBkgData&)obj).BkgMap= 0;
	}
	if( ((ImgBkgData&)obj).NoiseMap ){
		delete ((ImgBkgData&)obj).NoiseMap;
		((ImgBkgData&)obj).NoiseMap= 0;
	}

	//Copy maps
	DEBUG_LOG("Copying bkg & noise maps...");
	if(BkgMap){
		((ImgBkgData&)obj).BkgMap= new Image;
		*((ImgBkgData&)obj).BkgMap = *BkgMap;
	}
	if(NoiseMap){
		((ImgBkgData&)obj).NoiseMap= new Image;
		*((ImgBkgData&)obj).NoiseMap = *NoiseMap;
	}
	
}//close Copy()


ImgBkgData& ImgBkgData::operator=(const ImgBkgData& data) { 
	// Operator =
  if (this != &data)  ((ImgBkgData&)data).Copy(*this);
  return *this;
}

void ImgBkgData::CopyBkgMap(Caesar::Image* aMap){
	if(!aMap) return;
	TString mapName= "bkgMap";
	if(BkgMap) {
		mapName= BkgMap->GetName();
		delete BkgMap;
		BkgMap= 0;
	}

	bool copyMetaData= true;
	bool resetStats= false;
	BkgMap= aMap->GetCloned(std::string(mapName),copyMetaData,resetStats);

}//close CopyBkgMap()
		
void ImgBkgData::CopyNoiseMap(Caesar::Image* aMap){
	if(!aMap) return;
	TString mapName= "noiseMap";
	if(NoiseMap) {
		mapName= NoiseMap->GetName();	
		delete NoiseMap;
		NoiseMap= 0;
	}

	bool copyMetaData= true;
	bool resetStats= false;
	NoiseMap= aMap->GetCloned(std::string(mapName),copyMetaData,resetStats);

}//close CopyNoiseMap()
		
}//close namespace

