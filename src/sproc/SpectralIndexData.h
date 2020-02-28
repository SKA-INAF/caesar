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
//==   CLASS: SPECTRAL INDEX DATA
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

	public:
		void Reset(){Init();}

	private:
		/**
		* \brief Init data
		*/
		void Init();
		
	public:
		bool hasSpectralIndex;
		bool isMultiSourceMatchIndex;
		double norm;
		double normErr;
		double spectralIndex;
		double spectralIndexErr;
		bool isSpectralIndexFit;
		double spectralFitChi2;
		double spectralFitNDF;
		int spectralIndexFitStatus;

	ClassDef(SpectralIndexData,2)

};

#ifdef __MAKECINT__
#pragma link C++ class SpectralIndexData+;
#pragma link C++ class vector<SpectralIndexData>+;
#pragma link C++ class vector<SpectralIndexData*>+;
#endif

//============================================
//==   CLASS: PowerLaw SPECTRAL INDEX DATA
//============================================
class PLSpectralIndexData : public TObject  
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		PLSpectralIndexData();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~PLSpectralIndexData();

	private:
		/**
		* \brief Init data
		*/
		void Init();
		
	public:
		bool hasData;
		bool isMultiMatch;
		bool isFitted;
		double fitChi2;
		double fitNDF;	
		double norm;
		double normErr;
		double alpha;//spectral index
		double alphaErr;//spectral index err
		int status;//fit status
		
	ClassDef(PLSpectralIndexData,2)

};
#ifdef __MAKECINT__
#pragma link C++ class PLSpectralIndexData+;
#pragma link C++ class vector<PLSpectralIndexData>+;
#pragma link C++ class vector<PLSpectralIndexData*>+;
#endif


//============================================
//==   CLASS: Polynomial SPECTRAL INDEX DATA
//============================================
class PolSpectralIndexData : public TObject 
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		PolSpectralIndexData();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~PolSpectralIndexData();

	private:
		/**
		* \brief Init data
		*/
		void Init();
		
	public:
		bool hasData;
		bool isMultiMatch;
		double fitChi2;
		double fitNDF;
		std::vector<double> polPars;
		std::vector<double> polParErrors;
		int status;//fit status

	ClassDef(PolSpectralIndexData,2)

};
#ifdef __MAKECINT__
#pragma link C++ class PolSpectralIndexData+;
#pragma link C++ class vector<PolSpectralIndexData>+;
#pragma link C++ class vector<PolSpectralIndexData*>+;
#endif


//=====================================================
//==   CLASS: SYNCHROTRON + ABSORPTION SPECTRAL DATA
//=====================================================
class SASpectralIndexData : public TObject 
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SASpectralIndexData();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SASpectralIndexData();

	private:
		/**
		* \brief Init data
		*/
		void Init();
		
	public:
		bool hasData;
		bool isMultiMatch;
		double fitChi2;
		double fitNDF;
		double norm;
		double normErr;
		double alpha;
		double alphaErr;
		double nu_t;
		double nuErr_t;
		int status;//fit status
		
	ClassDef(SASpectralIndexData,2)

};
#ifdef __MAKECINT__
#pragma link C++ class SASpectralIndexData+;
#pragma link C++ class vector<SASpectralIndexData>+;
#pragma link C++ class vector<SASpectralIndexData*>+;
#endif




}//close namespace

#endif
