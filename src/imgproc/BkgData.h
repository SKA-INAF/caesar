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
* @file BkgData.h
* @class BkgData
* @brief BkgData
*
* Class for storing bkg data
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _BKG_DATA_h
#define _BKG_DATA_h 1

#include <Image.h>
#include <CodeUtils.h>
#include <Logger.h>

#include <TObject.h>
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

using namespace std;


namespace Caesar {

class BkgSampleData : public TObject{

	public:

		/**
		* \brief Constructor
		*/
		BkgSampleData(){
			ix_min= 0; ix_max= 0;
			iy_min= 0; iy_max= 0;
			npix= 0;
			isReliable= true;
			bkgLevel= 0;
			bkgRMS= 0;
		}

		/**
		* \brief Destructor
		*/
		virtual ~BkgSampleData(){};

	public:
		/**
		* \brief Copy bkg data
		*/
		void CopyBkgData(BkgSampleData aBkgSample){
			bkgLevel= aBkgSample.bkgLevel;
			bkgRMS= aBkgSample.bkgRMS;
		}

		/**
		* \brief Log info
		*/
		void Log(std::string level="INFO"){
			LOG(level,GetPrintable());
		}

		/**
		* \brief Print info to stdout
		*/
		void Print(){
			cout<<"== BKG SAMPLE DATA NO. "<<id<<" =="<<endl;
			cout<<"N="<<npix<<" xrange("<<ix_min<<","<<ix_max<<") yrange("<<iy_min<<","<<iy_max<<")"<<endl;
			cout<<"bkgLevel="<<bkgLevel<<" bkgRMS="<<bkgRMS<<endl;
			cout<<"=================================="<<endl;
		}

		/**
		* \brief Get printable string
		*/
		std::string GetPrintable(){
			std::stringstream ss;
			ss<<"BkgSample no. "<<id<<": ";
			ss<<"N="<<npix<<", xrange("<<ix_min<<","<<ix_max<<"), yrange("<<iy_min<<","<<iy_max<<"), ";
			ss<<"bkgLevel="<<bkgLevel<<", bkgRMS="<<bkgRMS;
			return ss.str();
		}

	public:	
		int id;
		int ix_min;
		int iy_min;
		int ix_max;
		int iy_max;
		int npix;
		bool isReliable;
		double bkgLevel;
		double bkgRMS;	

	ClassDef(BkgSampleData,1)
};


class ImgBkgData : public TObject {

	public:

		/**
		* \brief Standard constructor
		*/
		ImgBkgData();

		/**
		* \brief Destructor
		*/
		virtual ~ImgBkgData();

		/**
		* \brief Copy constructor
		*/
		ImgBkgData(const ImgBkgData& data);

		/**
		* \brief Assignment Operator
		*/
		ImgBkgData& operator=(const ImgBkgData& data);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& data) const;

	public:

		/**
		* \brief Clear sampling data
		*/
		void ClearSamplings(){
			BkgSamplings.clear();
		}

		/**
		* \brief Clear bkg map
		*/
		void ClearBkgMap(){
			if(!BkgMap) return;
			delete BkgMap;
			BkgMap= 0;
		}

		/**
		* \brief Clear noise map
		*/
		void ClearNoiseMap(){
			if(!NoiseMap) return;
			delete NoiseMap;
			NoiseMap= 0;
		}

		/**
		* \brief Clear data
		*/
		void Clear(){
			ClearSamplings();
			ClearBkgMap();
			ClearNoiseMap();
		}

		/**
		* \brief Copy bkg map
		*/
		void CopyBkgMap(Caesar::Image* aMap);		

		/**
		* \brief Copy noise map
		*/
		void CopyNoiseMap(Caesar::Image* aMap);

		/**
		* \brief Has local bkg
		*/
		bool HasLocalBkg(){return (BkgMap && NoiseMap);}
		
	private:
		/**
		* \brief Initialize class data
		*/
		void Init();

	public:
		std::vector<BkgSampleData> BkgSamplings;
		Image* BkgMap;//the interpolated bkg map
		Image* NoiseMap;//the interpolated noise map
		double gBkg;
		double gNoise;

	ClassDef(ImgBkgData,1)

};//close class ImgBkgData


#ifdef __MAKECINT__
#pragma link C++ class BkgSampleData+;
#pragma link C++ class vector<BkgSampleData>+;
#pragma link C++ class ImgBkgData+;
#endif


}//close namespace 


#endif


