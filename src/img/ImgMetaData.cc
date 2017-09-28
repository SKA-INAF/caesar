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
* @file ImgMetaData.cc
* @class ImgMetaData
* @brief ImgMetaData
*
* Image metadata class
* @author S. Riggi
* @date 20/01/2015
*/

#include <ImgMetaData.h>
#include <Consts.h>
#include <Logger.h>

//WCSTOOLS
#include <wcs.h>

ClassImp(Caesar::ImgMetaData)

namespace Caesar {

ImgMetaData::ImgMetaData()
{
	Init();
}

ImgMetaData::~ImgMetaData()
{

}//close destructor


void ImgMetaData::Init()
{
	//Initialize fields
	Nx= 0; Ny= 0;
	Cx= 0; Cy= 0;
	dX= 0; dY= 0;
	RotX= 0; RotY= 0;
	CoordTypeX= ""; CoordTypeY= "";
	BUnit= "";
	Bmaj= 0; Bmin= 0; Bpa= 0;
	Epoch= 0;
	m_wcsType= "";	

}//close Init()

void ImgMetaData::SetFITSCards(Caesar::FITSFileInfo& fits_info){
			
	//Initialize from FITS keywords
	Nx= (fits_info.header).Nx;
	Ny= (fits_info.header).Ny;
	Cx= (fits_info.header).Cx;
	Cy= (fits_info.header).Cy;
	Xc= (fits_info.header).Xc;
	Yc= (fits_info.header).Yc;
	CoordTypeX= (fits_info.header).CoordTypeX;
	CoordTypeY= (fits_info.header).CoordTypeY;
	BUnit= (fits_info.header).BUnit;
	Bmaj= (fits_info.header).Bmaj;
	Bmin= (fits_info.header).Bmin;
	Bpa= (fits_info.header).Bpa;
	dX= (fits_info.header).dX;
	dY= (fits_info.header).dY;
	RotX= (fits_info.header).RotX;
	RotY= (fits_info.header).RotY;
	Epoch= (fits_info.header).Epoch;
	
}//close SetFITSCards()


WorldCoor* ImgMetaData::GetWorldCoord(int coordSystem){

	//Compute the wcs from vars
	WorldCoor* wcs= wcskinit(Nx,Ny,(char*)CoordTypeX.c_str(),(char*)CoordTypeY.c_str(),Cx,Cy,Xc,Yc,NULL,dX,dY,RotY,(int)(Epoch),Epoch);
	std::string wcsType= std::string(getwcsout(wcs));
	
	//Convert wcs to desired type
	char* flag = (char*)("");
	if(coordSystem==eGALACTIC)
		flag = (char*)("GALACTIC");	
	else if(coordSystem==eJ2000)
		flag = (char*)("FK5");
	else if(coordSystem==eB1950)
		flag = (char*)("FK4");
	else if(coordSystem==-1 && m_wcsType!="")					
		flag = (char*)(m_wcsType.c_str());
			
	if(strcmp(flag,"")!=0) {
		wcsoutinit (wcs,flag);
		m_wcsType= std::string(flag);
	}
			
	wcsType= std::string(getwcsout(wcs));
	DEBUG_LOG("wcsType="<<wcsType);

	return wcs;
		
}//close GetWorldCoord()

}//close namespace
