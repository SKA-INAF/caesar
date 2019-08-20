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
* @file AstroObject.h
* @class AstroObject
* @brief AstroObject class
*
* AstroObject class
* @author S. Riggi
* @date 06/08/2019
*/


#ifndef _ASTRO_OBJECT_h
#define _ASTRO_OBJECT_h 1

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//ROOT headers
#include <TObject.h>
#include <TEllipse.h>

//C++ headers
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
#include <time.h>
#include <ctime>

namespace Caesar {

//=================================
//==    ASTRO OBJECT CLASS
//=================================
class AstroObject : public TObject
{
	public:
		/**
		* \brief Standard constructor
		*/
		AstroObject();
		/**
		* \brief Destructor
		*/
		virtual ~AstroObject();

	public:
		
		/**
		* \brief Get fit ellipse (if ellipse info are available)
		*/
		TEllipse* GetFitEllipse();
		
	protected:
		/**
		* \brief Init fields
		*/
		void Init();

	public:

		long int index;
		std::string name;//object name
		std::string id_str;
		int id;//major object identifier
		int subid;//minor object identifier
		double x;//RA
		double y;//Dec	
		double xerr;
		double yerr;
		std::string refs;//catalog references
		bool confirmed;//whether the object is confirmed or not
		
		bool hasFrequencyInfo;
		double nu;
		double dnu;

		bool hasFluxInfo;
		double fluxDensity;//Jy/beam
		double fluxDensityErr;//Jy/beam
		double flux;//Jy
		double fluxErr;//Jy

		bool hasEllipseInfo;
		double bmaj;//arcsec
		double bmin;//arcsec
		double pa;//deg
		
		bool hasDeconvEllipseInfo;
		double bmaj_deconv;//arcsec
		double bmin_deconv;//arcsec
		double pa_deconv;//deg

	ClassDef(AstroObject,2)

};//close class AstroObject

#ifdef __MAKECINT__
#pragma link C++ class AstroObject+;
#pragma link C++ class AstroObject*+;
#pragma link C++ class vector<AstroObject>+;
#pragma link C++ class vector<AstroObject*>+;
#endif

}//close namespace


#endif
