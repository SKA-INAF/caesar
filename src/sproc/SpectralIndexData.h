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
* @file SpectralIndexData.h
* @class SpectralIndexData
* @brief SpectralIndexData class
*
* Class representing spectral index data
* @author S. Riggi
* @date 26/03/2018
*/

#ifndef _SPECTRAL_INDEX_DATA_h
#define _SPECTRAL_INDEX_DATA_h 1


#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

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
#include <iostream>
#include <time.h>
#include <ctime>


namespace Caesar {


//======================================
//==      STRUCT: SPECTRAL INDEX DATA
//======================================
class SpectralIndexData : public TObject 
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SpectralIndexData();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SpectralIndexData();

		/**
		* \brief Copy constructor
		*/
		//SpectralIndexData(const SpectralIndexData& sid);
		
		/**
		* \brief Assignment Operator
		*/
		//SpectralIndexData& operator=(const SpectralIndexData& sid);
		/**
		* \brief Copy method
		*/
		//void Copy(TObject& obj) const;

	private:
		/**
		* \brief Init data
		*/
		void Init();
		
	public:
		bool hasSpectralIndex;
		bool isMultiSourceMatchIndex;
		double spectralIndex;
		double spectralIndexErr;
		bool isSpectralIndexFit;
		double spectralFitChi2;
		double spectralFitNDF;

	ClassDef(SpectralIndexData,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SpectralIndexData+;
#pragma link C++ class vector<SpectralIndexData>+;
#pragma link C++ class vector<SpectralIndexData*>+;
#endif

}//close namespace

#endif
