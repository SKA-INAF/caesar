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
* @file DS9RegionParser.cc
* @class DS9RegionParser
* @brief Class to parse DS9 region files
*
* Class to parse DS9 region files
* @author S. Riggi
* @date 27/02/2019
*/

#include <DS9RegionParser.h>
#include <DS9Region.h>
#include <MathUtils.h>
#include <CodeUtils.h>
#include <Contour.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//ROOT headers
#include <TMath.h>


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
#include <deque>

using namespace std;

ClassImp(Caesar::DS9RegionParser)

namespace Caesar {

DS9RegionParser::DS9RegionParser()
{

}//close constructor

DS9RegionParser::~DS9RegionParser()
{

}//close destructor

int DS9RegionParser::Parse(std::vector<DS9Region*>& regions,std::string filename)
{
	//## Init data
	regions.clear();	

	//## Read region file and get region text lines
	std::vector<std::vector<std::string>> raw_data;
	int wcsType= eUNKNOWN_CS; 
	if(Read(raw_data,wcsType,filename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read region data from file "<<filename<<"!");
		#endif
		return -1;
	}

	//## Check lines
	if(raw_data.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty data list read from file "<<filename<<"!");
		#endif
		return -1;
	}

	//## Extract regions from parsed lines
	bool parseErr= false;

	for(size_t i=0;i<raw_data.size();i++){
		DS9Region* region= ParseRegion(raw_data[i]);
		if(!region) {
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to parse region from data line no. "<<i+1<<", stop parsing and return error!");
			#endif
			parseErr= true;
			break;
		}
	
		//Add WCS info and append region to list
		region->csType= wcsType;
		regions.push_back(region);

	}//end loop data

	//Clear data in case of errors
	if(parseErr){
		CodeUtils::DeletePtrCollection<DS9Region>(regions);
		return -1;
	}

	return 0;

}//close Parse()


DS9Region* DS9RegionParser::ParseRegion(const std::vector<std::string>& fields)
{
	//If parsing is correct there should be only one entry in fields
	if(fields.size()!=1){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty or invalid number of fields found (1 expected)!");
		#endif
		return nullptr;
	}

	//Detect type of region
	std::string data_str= fields[0];
	DS9Region* region= 0;

	if(CodeUtils::HasPatternInString(data_str,"polygon")){
		region= ParsePolygonRegion(data_str);
	}
	else if(CodeUtils::HasPatternInString(data_str,"box")){
		region= ParseBoxRegion(data_str);
	}
	else if(CodeUtils::HasPatternInString(data_str,"circle")){
		region= ParseCircleRegion(data_str);
	}
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Unknown/unsupported region type found!");
		#endif
		return nullptr;
	}


	return region;

}//close ParseRegion()


DS9PolygonRegion* DS9RegionParser::ParsePolygonRegion(const std::string& data_str)
{
	//Check string 
	if(data_str=="") return nullptr;

	//Parse string and get coord data
	bool hasComment= CodeUtils::HasPatternInString(data_str,"#");
	std::string delim_start= "polygon(";
	std::string delim_stop= ")";
	size_t first= data_str.find(delim_start);
	size_t last = data_str.find(delim_stop);
	if(first==std::string::npos){
		delim_start= "polygon";//try with another region format (without parenthesis)
		first= data_str.find(delim_start);
		if(first==std::string::npos){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to find polygon stirng pattern (this should not occur!)");
			#endif
			return nullptr;
		}
	}
	if(last==std::string::npos){
		if(hasComment) {
			delim_stop= "#";
			last = data_str.find(delim_stop);
		}
		else {
			last= data_str.size();
		} 
	}
	
	std::string parsed_data= data_str.substr(first+delim_start.size(),last-first-delim_start.size());

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("polygon parsed_data="<<parsed_data);
	#endif

	//Check if data are not given in real digit format (e.g. sexagesimal not supported)
	//bool hasPlusSigns= CodeUtils::HasPatternInString(parsed_data,"+");
	//bool hasMinusSigns= CodeUtils::HasPatternInString(parsed_data,"-");
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
	//if(hasPlusSigns || hasMinusSigns || hasColon){
	if(hasColon){
		#ifdef LOGGING_ENABLED		
			ERROR_LOG("Sexagesimal data format not supported, use degrees!");
		#endif
		return nullptr;
	}

	//Check if data are separated by ',' or by whitespaces
	bool hasCommaDelimiter= CodeUtils::HasPatternInString(parsed_data,",");
	std::vector<std::string> coords_str;

	if(hasCommaDelimiter){
		coords_str= CodeUtils::SplitStringOnPattern(parsed_data,',');
	}
	else{//splitting data on whitespaces
		coords_str= CodeUtils::SplitStringOnWhitespaces(parsed_data);
	}
	
	//Check coords
	if(coords_str.empty() || coords_str.size()%2==1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid number of coordinates (must be odd!)");
		#endif
		return nullptr;
	}

	//Extract coordinates and convert to double
	DS9PolygonRegion* region= new DS9PolygonRegion();
	for(size_t i=0;i<coords_str.size();i+=2){
		double x= atof(coords_str[i].c_str());	
		double y= atof(coords_str[i+1].c_str());
		(region->points).push_back(TVector2(x,y));
	}	

	return region;

}//close ParsePolygonRegion()

DS9BoxRegion* DS9RegionParser::ParseBoxRegion(const std::string& data_str)
{
	//Check string 
	if(data_str=="") return nullptr;

	//Parse string and get coord data
	bool hasComment= CodeUtils::HasPatternInString(data_str,"#");
	std::string delim_start= "box(";
	std::string delim_stop= ")";
	size_t first= data_str.find(delim_start);
	size_t last = data_str.find(delim_stop);
	if(first==std::string::npos){
		delim_start= "box";//try with another region format (without parenthesis)
		first= data_str.find(delim_start);
		if(first==std::string::npos){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to find polygon stirng pattern (this should not occur!)");
			#endif
			return nullptr;
		}
	}
	if(last==std::string::npos){
		if(hasComment) {
			delim_stop= "#";
			last = data_str.find(delim_stop);
		}
		else {
			last= data_str.size();
		} 
	}
	
	std::string parsed_data= data_str.substr(first+delim_start.size(),last-first-delim_start.size());

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("box parsed_data="<<parsed_data);
	#endif

	//Check if data are not given in real digit format (e.g. sexagesimal not supported)
	//bool hasPlusSigns= CodeUtils::HasPatternInString(parsed_data,"+");
	//bool hasMinusSigns= CodeUtils::HasPatternInString(parsed_data,"-");
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
	//if(hasPlusSigns || hasMinusSigns || hasColon){
	if(hasColon){
		#ifdef LOGGING_ENABLED		
			ERROR_LOG("Sexagesimal data format not supported, use degrees!");
		#endif
		return nullptr;
	}

	//Check if data are separated by ',' or by whitespaces
	bool hasCommaDelimiter= CodeUtils::HasPatternInString(parsed_data,",");
	std::vector<std::string> coords_str;

	if(hasCommaDelimiter){
		coords_str= CodeUtils::SplitStringOnPattern(parsed_data,',');
	}
	else{//splitting data on whitespaces
		coords_str= CodeUtils::SplitStringOnWhitespaces(parsed_data);
	}
	
	//Check coords
	if(coords_str.empty() || coords_str.size()!=5){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid number of coordinates (must be =5!)");
		#endif
		return nullptr;
	}

	//Extract coordinates and convert to double
	DS9BoxRegion* region= new DS9BoxRegion();
	double cx= atof(coords_str[0].c_str());	
	double cy= atof(coords_str[1].c_str());
	double width= atof(coords_str[2].c_str());	
	double height= atof(coords_str[3].c_str());
	double theta= atof(coords_str[4].c_str());
	region->cx= cx;
	region->cy= cy;
	region->width= width;
	region->height= height;
	region->theta= theta;
	region->ComputeBoxCoords();

	return region;

}//close ParseBoxRegion()


DS9CircleRegion* DS9RegionParser::ParseCircleRegion(const std::string& data_str)
{
	//Check string 
	if(data_str=="") return nullptr;

	//Parse string and get coord data
	bool hasComment= CodeUtils::HasPatternInString(data_str,"#");
	std::string delim_start= "circle(";
	std::string delim_stop= ")";
	size_t first= data_str.find(delim_start);
	size_t last = data_str.find(delim_stop);
	if(first==std::string::npos){
		delim_start= "circle";//try with another region format (without parenthesis)
		first= data_str.find(delim_start);
		if(first==std::string::npos){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to find polygon stirng pattern (this should not occur!)");
			#endif
			return nullptr;
		}
	}
	if(last==std::string::npos){
		if(hasComment) {
			delim_stop= "#";
			last = data_str.find(delim_stop);
		}
		else {
			last= data_str.size();
		} 
	}
	
	std::string parsed_data= data_str.substr(first+delim_start.size(),last-first-delim_start.size());

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("circle parsed_data="<<parsed_data);
	#endif

	//Check if data are not given in real digit format (e.g. sexagesimal not supported)
	//bool hasPlusSigns= CodeUtils::HasPatternInString(parsed_data,"+");
	//bool hasMinusSigns= CodeUtils::HasPatternInString(parsed_data,"-");
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
	//if(hasPlusSigns || hasMinusSigns || hasColon){
	if(hasColon){
		#ifdef LOGGING_ENABLED		
			ERROR_LOG("Sexagesimal data format not supported, use degrees!");
		#endif
		return nullptr;
	}

	//Check if data are separated by ',' or by whitespaces
	bool hasCommaDelimiter= CodeUtils::HasPatternInString(parsed_data,",");
	std::vector<std::string> coords_str;

	if(hasCommaDelimiter){
		coords_str= CodeUtils::SplitStringOnPattern(parsed_data,',');
	}
	else{//splitting data on whitespaces
		coords_str= CodeUtils::SplitStringOnWhitespaces(parsed_data);
	}
	
	//Check coords
	if(coords_str.empty() || coords_str.size()!=3){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid number of coordinates (must be =3!)");
		#endif
		return nullptr;
	}

	//Check radius units
	double convFactor= 1;//conversion factor to degrees
	bool hasUnits= CodeUtils::HasPatternInString(coords_str[2],"\"");
	std::string radiusStr= coords_str[2];
	if(hasUnits) {
		convFactor= 1./3600.;
		CodeUtils::RemovePatternInString(radiusStr,"\"");
		#ifdef LOGGING_ENABLED		
			DEBUG_LOG("DS9 circle region radius found in arcsec units, will convert to degree...");
		#endif
	}
	#ifdef LOGGING_ENABLED		
		DEBUG_LOG("DS9 circle region radius parsed: "<<radiusStr);
	#endif
		
	//Extract coordinates and convert to double
	DS9CircleRegion* region= new DS9CircleRegion();
	double cx= atof(coords_str[0].c_str());	
	double cy= atof(coords_str[1].c_str());
	double r= convFactor*atof(coords_str[2].c_str());	
	region->cx= cx;
	region->cy= cy;
	region->r= r;
	
	return region;

}//close ParseCircleRegion()


int DS9RegionParser::Read(std::vector<std::vector<std::string>>& data, int& wcsType,std::string filename)
{
	//Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filename<<" for reading!");
		#endif
		return -1;
  }

	//Parsing file
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading and parsing file: "<<filename);
	#endif

	std::string parsedline= "";	
	int line_counter= 0;
	std::map<std::string,int> wcsTypeMap;
	wcsTypeMap.insert(std::make_pair("image",eIMG_CS));
	wcsTypeMap.insert(std::make_pair("fk5",eJ2000));
	wcsTypeMap.insert(std::make_pair("fk4",eB1950));
	wcsTypeMap.insert(std::make_pair("galactic",eGALACTIC));
	
	std::vector<std::string> skipLinePatterns {"global","image","fk5","fk4","galactic"};

	wcsType= eUNKNOWN_CS;
	bool foundWCS= false;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;

		//Parse line fields
		istringstream ss(parsedline);
		std::string field= "";
		bool skipLine= false;
		std::vector<std::string> fields;

		while(ss >> field)
    {
			//Skip first line
			char first_char= *(field.c_str());
			if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
				skipLine= true;
				break;
			}

			//Search for WCS pattern
			for (auto const& it : wcsTypeMap){
				bool hasPattern= CodeUtils::HasPatternInString(field,it.first);
				if(hasPattern){
					wcsType= it.second;
					foundWCS= true;
					break;
				}
			}

			//Search for pattern and skip if present
			for(size_t i=0;i<skipLinePatterns.size();i++){
				bool hasPattern= CodeUtils::HasPatternInString(field,skipLinePatterns[i]);
				if(hasPattern) {
					skipLine= true;
					break;
				}
			}
			if(skipLine) break;
	
			//cout<<field<<"  ";
			fields.push_back(field);

		}//end loop line
		
		if(skipLine) continue;
		
		//Add region line to list
		data.push_back(fields);
		
		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	//Check if WCS type was found
	if(!foundWCS){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cannot find and parse WCS type in region file (check if given WCS is supported)!");
		#endif
		return -1;
	}

	return 0;

}//close Read()


}//close namespace
