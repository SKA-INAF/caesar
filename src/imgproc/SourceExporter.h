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
* @file SourceExporter.h
* @class SourceExporter
* @brief SourceExporter class
*
* Class to export an image source in different formats
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SOURCE_EXPORTER_h
#define _SOURCE_EXPORTER_h 1

#include <TObject.h>
#include <TMatrixD.h>

#include <Consts.h>

//WCSTOOLS (TO BE DEPRECATED)
//#include <wcs.h>

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
class WCS;

class SourceExporter : public TObject 
{
  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SourceExporter();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~SourceExporter();

		
	public:

		//=======================================
		//==      ASCII EXPORT
		//=======================================
		/**
		* \brief Write ascii file from source collection
		*/
		//static int WriteToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WorldCoor* wcs=0);
		static int WriteToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		
		/**
		* \brief Get source ascii string
		*/
		//static const std::vector<std::string> SourceToAscii(Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WorldCoor* wcs=0);
		static const std::vector<std::string> SourceToAscii(Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);

		/**
		* \brief Write ascii file from source component collection
		*/
		//static int WriteComponentsToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WorldCoor* wcs=0);
		static int WriteComponentsToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		
		/**
		* \brief Get source component ascii string
		*/
		//static const std::vector<std::string> SourceComponentsToAscii(Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WorldCoor* wcs=0);
		static const std::vector<std::string> SourceComponentsToAscii(Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
	
		//=======================================
		//==      DS9 EXPORT
		//=======================================
		/**
		* \brief Write DS9 regions from source collection
		*/
		//static int WriteToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS=false,int ds9WCSType=eJ2000,int ds9RegionFormat=ePolygonRegion,WorldCoor* wcs=0);
		static int WriteToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS=false,int ds9WCSType=eJ2000,int ds9RegionFormat=ePolygonRegion,WCS* wcs=0);

		/**
		* \brief Write DS9 regions for source fitted components
		*/
		//static int WriteComponentsToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS=false,int ds9WCSType=eJ2000,WorldCoor* wcs=0);
		static int WriteComponentsToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS=false,int ds9WCSType=eJ2000,WCS* wcs=0);


		/**
		* \brief Get DS9 region info
		*/
		//static const std::string SourceToDS9Region(Source* source,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WorldCoor* wcs=0,int coordSystem=-1);
		static const std::string SourceToDS9Region(Source* source,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1);

		/**
		* \brief Get DS9 ellipse info
		*/
		static const std::string SourceToDS9EllipseRegion(Source* source,bool dumpNestedSourceInfo=false);
	
		/**
		* \brief Get DS9 fitted ellipse info
		*/
		//static const std::string SourceToDS9FittedEllipseRegion(Source* source,bool useFWHM=true,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WorldCoor* wcs=0,int coordSystem=-1);
		static const std::string SourceToDS9FittedEllipseRegion(Source* source,bool useFWHM=true,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1);

	
	private:

		/**
		* \brief Get DS9 region color according to source type
		*/
		static std::string GetDS9RegionColor(Source* source);

	private:
	
		ClassDef(SourceExporter,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class SourceExporter+;
#endif

}//close namespace 


#endif
