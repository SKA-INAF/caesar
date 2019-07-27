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
* @file SourceGroup.h
* @class SourceGroup
* @brief SourceGroup class
*
* Class representing a source group
* @author S. Riggi
* @date 26/03/2018
*/

#ifndef _SOURCE_GROUP_h
#define _SOURCE_GROUP_h 1

//#include <Source.h>
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

class Source;

//======================================
//==      CLASS: SOURCE GROUP
//======================================
class SourceGroup : public TObject 
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceGroup();

		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceGroup();

		/**
		* \brief Copy constructor
		*/
		//SourceGroup(const SourceGroup& sgroup);
		/**
		* \brief Assignment Operator
		*/
		//SourceGroup& operator=(const SourceGroup& sgroup);
		/**
		* \brief Copy method
		*/
		//void Copy(TObject& obj) const;
	
	public:
		/**
		* \brief Add source to group
		*/
		int AddSource(Source* aSource,int clone=false);
		/**
		* \brief Get size of group
		*/
		int GetNSources(){return static_cast<int>(m_sources.size());}


		/**
		* \brief Get sources
		*/
		std::vector<Source*>& GetSources(){return m_sources;}

		/**
		* \brief Get source i-th
		*/
		Source* GetSource(int index)
		{
			if(m_sources.empty()) return nullptr;
			if(index<0 || index>=(int)(m_sources.size())) return nullptr;
			return m_sources[index];
		}

		/**
		* \brief Get source names
		*/
		int GetSourceNames(std::vector<std::string>& snames);

		/**
		* \brief Get flux and its error. If Nsources>2 return sum
		*/
		int GetFlux(double& flux,double& fluxErr,bool& summed);

		/**
		* \brief Get frequency and relative width
		*/
		int GetFrequency(double& freq,double& dfreq);

	private:
		/**
		* \brief Init class data
		*/
		void Init();
		/**
		* \brief Delete class data
		*/
		void Clear();
		
		
	private:
		// - Source collection
		std::vector<Source*> m_sources;

		bool m_owned;

	ClassDef(SourceGroup,1)


};//close class SourceGroup

#ifdef __MAKECINT__
#pragma link C++ class SourceGroup+;
#pragma link C++ class vector<SourceGroup>+;
#pragma link C++ class vector<SourceGroup*>+;
#endif

}//close namespace

#endif
