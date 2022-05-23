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
* @file SourceImporter.cc
* @class SourceImporter
* @brief SourceImporter class
*
* Class to import image sources from different formats
* @author S. Riggi
* @date 17/05/2022
*/


#include <SourceImporter.h>
#include <Source.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <AstroUtils.h>
#include <WCSUtils.h>
#include <Consts.h>
#include <ImgMetaData.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TEllipse.h>

//WCSTOOLS (TO BE DEPRECATED)
//#include <wcs.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::SourceImporter)

namespace Caesar {

SourceImporter::SourceImporter()
{

}

SourceImporter::~SourceImporter()
{

}

//=================================================
//==        JSON IMPORTER
//=================================================
int SourceImporter::ImportFromJson(std::string filename, std::vector<Source*>& sources, Image* img)
{
	// - Init input collection
	sources.clear();

	// - Read json file
	Json::Value root;
	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading and parsing json file " << filename <<" ...");
	#endif
	try{
		std::ifstream fin(filename);
		fin >> root;
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse file "<< filename << "!");
		#endif
		return -1;
	}

	// - Read image
	if(!img){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading image listed in json file " << filename <<" ...");
		#endif
		img= ReadImageFromJson(root);
	}

	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image or failed to read image from json file!");
		#endif
		return -1;
	}

	// - Read sources
	auto root_slist= root["sources"];
	if(root_slist.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source list is empty in parsed json file, returning empty source collection!");
		#endif
		return 0;
	}

	
	for (Json::Value::ArrayIndex i=0; i!=root_slist.size(); i++)
	{
		Source* source= ReadSourceFromJson(root_slist[i], img);
		if(!source){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to read source no. "<<i+1<<" from json, skip to next...");
			#endif
			continue;
		}

		sources.push_back(source);

	}//end loop sources

	return 0;

}//close ImportFromJson()


Image* SourceImporter::ReadImageFromJson(Json::Value& root)
{
	//- Get image path from json
	std::string filename= "";
	try {
		filename= root["metadata"]["inputdata"]["filename"].asString();
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get image filename from json!");
		#endif
		return nullptr;
	}	

	//- Read image
	Image* img= new Image;
	if(img->ReadFITS(filename)<0)
	{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read image from file "<< filename << "!");
		#endif
		delete img;
		img= 0;
		return nullptr;
	}

	return img;

}//close ReadImageFromJson()

Source* SourceImporter::ReadSourceFromJson(Json::Value& json, Image* img)
{
	// - Retrieve some fields
	std::string sname= json["name"].asString();
	int sindex= json["index"].asInt();

	// - Check number of islands
	unsigned int island_size= json["islands"].size();
	unsigned int nislands= json["nislands"].asUInt();
	if(island_size==0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No islands for source "<<sname<<", this cannot occur, returning nullptr!");
		#endif
		return nullptr;
	}
	if(island_size!=nislands){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Mismatch between nislands field and island array size for source "<<sname<<", this cannot occur, returning nullptr!");
		#endif
		return nullptr;
	}

	// - Fill source data
	//   NB: If only one island is present and no parent, create standard source, otherwise create one source with nested 
	//Source* source= new Source(sname);
	Source* source= nullptr;
	auto root_ilist= json["islands"];

	if(nislands==1)
	{
		source= ReadSourceIslandFromJson(root_ilist[0], img);
	}
	else
	{
		//Create source
		source= new Source(sname);

		//Loop over islands
		std::vector<Source*> nested_sources;
		for (Json::Value::ArrayIndex i=0; i!=root_ilist.size(); i++)
		{
			Source* source_nested= ReadSourceIslandFromJson(root_ilist[i], img);
			if(!source_nested){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to read island no. "<<i+1<<" of source "<<sname<<" from json, skip to next...");
				#endif
				continue;
			}
			nested_sources.push_back(source_nested);

			//Add nested source
			source->AddNestedSource(source_nested);

			//Merge nested source into mother
			bool copyPixels= true;
			bool checkIfAdjacent= false;
			bool computeStatPars= false;
			bool computeMorphPars= false;
			bool sumMatchingPixels= false;
			source->MergeSource(
				source_nested, 
				copyPixels, checkIfAdjacent,
				computeStatPars, computeMorphPars, 
				sumMatchingPixels
			);

		}//end loop islands

		//Set composite source flag
		source->SetAreNestedComponentsOfCompositeSource(true);

		//Recompute merged source stats
		bool computeRobustStats= true;
		bool forceRecomputing= true;
		source->ComputeStats(computeRobustStats,forceRecomputing);

		//Recompute merged source pars
		source->ComputeMorphologyParams();
			
	}//close else

	return source;

}//close ReadSourceFromJson()

Source* SourceImporter::ReadSourceIslandFromJson(Json::Value& json, Image* img)
{
	// - Return if no image is given
	if(!img){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Input image is null, cannot retrieve pixels and create source, returning nullptr ...");
		#endif
		return nullptr;
	}

	// - Retrieve image metadata if available
	ImgMetaData* metadata= img->GetMetaData();
	float Xmin= img->GetXmin();
	float Ymin= img->GetYmin();
	float Xmax= img->GetXmax();
	float Ymax= img->GetYmax();

	// - Check vertices and/or pixel list field are present
	std::string sname= json["name"].asString();
	auto root_vertices= json["vertices"];
	auto root_pixels= json["pixels"];
	bool has_pixels= !root_pixels.empty();
	bool has_contour= !root_vertices.empty();
	if(!has_pixels && !has_contour){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source island "<<sname<<" has empty vertices & pixel fields, cannot create source, returning nullptr ...");
		#endif
		return nullptr;
	}

	// - Init source
	Source* source= new Source(sname);

	// - Retrieve pixels
	std::vector<Pixel*> pixels;
	Pixel* pixel= nullptr;
	bool hasBorderPixels= false;
	if(has_pixels)
	{
		for (Json::Value::ArrayIndex i=0; i!=root_pixels.size(); i++)
		{
			long int ix= root_pixels[i][0].asInt64();
			long int iy= root_pixels[i][1].asInt64();
			bool hasBin= img->HasBin(ix,iy);
			if(!hasBin){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Pixel "<<i+1<<" of source island "<<sname<<" does not exist in image, skipping it ...");
				#endif
				continue;
			}
			double x= img->GetX(ix);
			double y= img->GetX(iy);
			long int gbin= img->GetBin(ix,iy);
			double S= img->GetPixelValue(ix,iy);
			if( img->IsEdgeBin(ix,iy) ) {
				source->SetEdgeFlag(true);
				hasBorderPixels= true;
			}

			pixel= new Pixel;
			pixel->S= S;
			pixel->id= gbin;
			pixel->SetPhysCoords(x,y);
			pixel->SetCoords(ix,iy);

			source->AddPixel(pixel);
		}
	}
	else
	{//retrieve pixels from vertices

	}//close else

	
	// - Compute source stats
	source->ComputeStats();

	// - Compute source morph pars
	source->ComputeMorphologyParams();

	// - Adding image metadata to image (needed for WCS)
	source->SetImageMetaData(metadata, Xmin, Ymin);


	//- Set Bkg/noise estimators
	//  NB: pixel bkg & rms are not stored in json, so we need to reset them
	double nPixels= static_cast<double>(source->NPix);
	double bkgLevel= json["bkg"].asDouble();
	double bkgRMS= json["rms"].asDouble();
	double bkgSum= nPixels*bkgLevel;
	double bkgRMSSum= nPixels*bkgRMS;
	source->SetBkgSum(bkgSum);
	source->SetBkgRMSSum(bkgRMSSum);
	
	//- Source flags
	std::string morph_label= json["morph_label"].asString();
	int sourceMorphId= GetSourceMorphId(morph_label);
	source->MorphId= sourceMorphId;

	std::string sourceness_label= json["sourceness_label"].asString();
	int sourcenessId= GetSourcenessId(sourceness_label);
	source->SourcenessId= sourcenessId;

	bool at_edge= json["border"].asBool();
	source->SetEdgeFlag(at_edge);

	/*
	json["sourceness_score"]= -1;//not assessed by default
		
	//- User tags (empty list by default)
	json["tags"]= Json::arrayValue;
	*/

	return source;

}//close ReadSourceIslandFromJson()

}//close namespace
