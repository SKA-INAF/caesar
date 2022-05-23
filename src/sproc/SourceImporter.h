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
* @file SourceImporter.h
* @class SourceImporter
* @brief SourceImporter class
*
* Class to import image sources from different formats
* @author S. Riggi
* @date 17/05/2022
*/

#ifndef _SOURCE_IMPORTER_h
#define _SOURCE_IMPORTER_h 1


#include <Consts.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TTree.h>

#include <json/json.h>

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
class Image;

struct SourceInfoData 
{
	Source* source;
	int index;
	int parent_index;
	int parent_island_index;

	SourceInfoData(){
		source= nullptr;
		parent_index= -1;
		parent_island_index= -1;
	};
	SourceInfoData(Source* _source, int _index, int _pindex, int _pindex_island):
		source(_source), index(_index), parent_index(_pindex), parent_island_index(_pindex_island)
	{
		
	};
};

struct SourceInfoDataFinder 
{
	bool operator()(const SourceInfoData* lhs, const SourceInfoData* rhs) const { 
		bool equal= (
			lhs->index==rhs->index &&
			lhs->parent_index==rhs->parent_index &&
			lhs->parent_island_index==rhs->parent_island_index
		);
		return equal;
	}
};


class SourceImporter : public TObject 
{
  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SourceImporter();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~SourceImporter();

		
	public:

		//=======================================
		//==      JSON IMPORT
		//=======================================
		/**
		* \brief Import source collection from json
		*/
		static int ImportFromJson(std::string filename, std::vector<Source*>& sources, Image* img=0);

	private:
		//=======================================
		//==      JSON IMPORT
		//=======================================
		/**
		* \brief Read image from json fields
		*/
		static Image* ReadImageFromJson(Json::Value& json);
		/**
		* \brief Read source from json fields
		*/
		static Source* ReadSourceFromJson(Json::Value& json, Image* img=0);
		/**
		* \brief Read source island from json fields
		*/
		static Source* ReadSourceIslandFromJson(Json::Value& json, Image* img=0);


private:
	
		ClassDef(SourceImporter,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class SourceImporter+;
#endif

}//close namespace 


#endif

