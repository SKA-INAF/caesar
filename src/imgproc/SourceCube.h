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
* @file SourceCube.h
* @class SourceCube
* @brief SourceCube class
*
* Class representing an image source cube
* @author S. Riggi
* @date 26/03/2018
*/

#ifndef _SOURCE_CUBE_h
#define _SOURCE_CUBE_h 1

#include <Logger.h>

#include <TObject.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

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
//==      CLASS: SOURCE CUBE
//======================================
class SourceCube : public TNamed {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceCube();
		/** 
		\brief Parametric constructor
 		*/
		SourceCube(std::string name);		

		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceCube();

		/**
		* \brief Copy constructor
		*/
		SourceCube(const SourceCube& scube);
		
		/**
		* \brief Assignment Operator
		*/
		SourceCube& operator=(const SourceCube& scube);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& scube) const;
	
		
	public:

		/**
		* \brief Add source to cube
		*/
		int AddSource(Source* aSource){
			if(!aSource) return -1;
			m_sources.push_back(aSource);
			return 0;
		}

		/**
		* \brief Add component
		*/
		void AddComponent(){
			m_componentIndexes.push_back( std::vector<std::pair<size_t,size_t>>() );
		}

		/**
		* \brief Add index to existing component
		*/
		int AddIndexToComponent(int cubeComponentIndex,size_t sindex,size_t componentIndex){
			//Check if cube component id was allocated
			int nComponents= static_cast<int>(m_componentIndexes.size());
			if(nComponents<=0 || nComponents<cubeComponentIndex+1){
				ERROR_LOG("Component with index "<<cubeComponentIndex<<" was not allocated!");
				return -1;
			}
	
			//Fill component indexes
			m_componentIndexes[cubeComponentIndex].push_back(std::make_pair(sindex,componentIndex));
				
			return 0;
		}

		/**
		* \brief Draw source images
		*/
		int DoSourceImagePlot(bool useWCS=true,int coordSyst=0);

		/**
		* \brief Draw source SEDs
		*/
		int DoSourceSEDs();	

	private:
		
		/**
		* \brief Initialize class members
		*/
		void Init();

	protected:
	
		//- List of sources in cube
		std::vector<Source*> m_sources;

		//- List of component matches in cube
		std::vector<std::vector<std::pair<size_t,size_t>>> m_componentIndexes; 	

		//- Canvas with source plots	
		TCanvas* m_sourcePlot;

		//- Graph with source SED
		TGraphAsymmErrors* m_sourceSED;
		std::vector<TGraphAsymmErrors*> m_sourceComponentSED;

		//- Canvas with SED plot
		TCanvas* m_sourceSEDPlot;

	ClassDef(SourceCube,3)


};//close class SourceCube

#ifdef __MAKECINT__
#pragma link C++ class SourceCube+;
#pragma link C++ class vector<SourceCube>+;
#pragma link C++ class vector<SourceCube*>+;
#endif


}//close namespace

#endif


