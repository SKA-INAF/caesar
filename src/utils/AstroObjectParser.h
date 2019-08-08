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
* @file AstroObjectParser.h
* @class AstroObjectParser
* @brief AstroObjectParser class
*
* AstroObjectParser class
* @author S. Riggi
* @date 06/08/2019
*/


#ifndef _ASTRO_OBJECT_PARSER_h
#define _ASTRO_OBJECT_PARSER_h 1

#include <AstroObject.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//ROOT headers
#include <TObject.h>

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
//==    ASTRO OBJECT PARSER CLASS
//=================================
class AstroObjectParser : public TObject
{
	public:
		/**
		* \brief Standard constructor
		*/
		AstroObjectParser();
		/**
		* \brief Destructor
		*/
		virtual ~AstroObjectParser();

	public:
		/**
		* \brief Get objects from SIMBAD ascii catalog
		*/
		static int ParseSimbadData(std::vector<AstroObject*>& astroObjects,std::string filename,char delimiter='|');
		/**
		* \brief Get objects from MGPS ascii catalog
		*/
		static int ParseMGPSData(std::vector<AstroObject*>& astroObjects,std::string filename,char delimiter='|');
	

	protected:
		/**
		* \brief Parse SIMBAD ascii catalog line
		*/
		static int ParseSimbadObjectData(AstroObject& astroObject,std::string data,char delimiter='|');
		/**
		* \brief Parse MGPS ascii catalog line
		*/
		static int ParseMGPSObjectData(AstroObject& astroObject,std::string data,char delimiter='|');
		/**
		* \brief Read catalog ascii data
		*/
		static int Read(std::vector<std::string>& data,std::string filename);

	protected:

		static std::map<std::string,int> m_simbadObjIdMap;
		static std::map<std::string,int> m_simbadObjSubIdMap;


	ClassDef(AstroObjectParser,1)

};//close class AstroObject

#ifdef __MAKECINT__
#pragma link C++ class AstroObjectParser+;
#pragma link C++ class vector<AstroObjectParser>+;
#pragma link C++ class vector<AstroObjectParser*>+;
#endif

}//close namespace


#endif
