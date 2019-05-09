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
	std::vector<std::vector<std::string>> raw_metadata;
	int wcsType= eUNKNOWN_CS; 
	if(Read(raw_data,raw_metadata,wcsType,filename)<0){
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
	if(!raw_metadata.empty() && raw_metadata.size()!=raw_data.size()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Region metadata read from "<<filename<<" but differs in size wrt region data!");
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
	
		//Add WCS info
		region->csType= wcsType;
		
		//Add metadata info
		if(!raw_metadata[i].empty()){	
			DS9RegionMetaData metadata;
			int status= ParseRegionMetaData(metadata,raw_metadata[i]);
			if(status==0){
				region->SetMetaData(metadata);
			}
			else{
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to parse region metadata from data line no. "<<i+1<<", will not set metadata for this region!");
				#endif
			}
		}
		else{
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Empty region raw metadata from data line no. "<<i+1<<", will not set metadata for this region!");
			#endif
		}

		//Append region to list
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
	//Check if empty fields
	if(fields.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty number of fields given!");
		#endif
		return nullptr;
	}

	/*
	//If parsing is correct there should be only one entry in fields
	if(fields.size()!=1){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty or invalid number of fields found (1 expected, "<<fields.size()<<" found)!");
		#endif
		return nullptr;
	}
	*/

	//Detect type of region
	//std::string data_str= fields[0];
	std::string data_str= CodeUtils::JoinVec(fields," ");

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
	else if(CodeUtils::HasPatternInString(data_str,"ellipse")){
		region= ParseEllipseRegion(data_str);
	}
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Unknown/unsupported region type found!");
		#endif
		return nullptr;
	}


	return region;

}//close ParseRegion()


int DS9RegionParser::ParseRegionMetaData(DS9RegionMetaData& metadata,const std::vector<std::string>& fields)
{
	//Define search string
	std::string name_delim_start= "text={";
	std::string name_delim_stop= "}";
	std::string tag_delim_start= "tag={";
	std::string tag_delim_stop= "}";
	std::vector<std::string> sourceTypes {"unknown-type","compact","point-like","extended","compact-extended"};
	std::vector<std::string> sourceFlags {"unknown-flag","real","candidate","fake"};
	std::vector<std::string> sourceFitQualities {"uq-fit","hq-fit","mq-fit","lq-fit","bad-fit"};

	//Loop over metadata and parse
	for(size_t i=0;i<fields.size();i++)
	{
		std::string data_str= fields[i];

		//Parse region name if present
		size_t first= data_str.find(name_delim_start);
		size_t last = data_str.find(name_delim_stop);
		if(first!=std::string::npos && last!=std::string::npos){
			std::string name= data_str.substr(first+name_delim_start.size(),last-first-name_delim_start.size());
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("regionName="<<name);
			#endif
			metadata.sourceName= name;
			metadata.hasSourceName= true;
		}	

		//Parse region tag
		first= data_str.find(tag_delim_start);
		last = data_str.find(tag_delim_stop);
		if(first!=std::string::npos && last!=std::string::npos){
			std::string tag= data_str.substr(first+tag_delim_start.size(),last-first-tag_delim_start.size());

			//Parse source type
			for(size_t j=0;j<sourceTypes.size();j++){
				bool hasPattern= CodeUtils::HasPatternInString(tag,sourceTypes[j]);
				if(!hasPattern) continue;
				int sourceType= GetSourceType(sourceTypes[j]); 
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("sourceType="<<sourceType);
				#endif
				metadata.sourceType= sourceType;	
				metadata.hasSourceType= true;
			}
		
			//Parse source flag
			for(size_t j=0;j<sourceFlags.size();j++){
				bool hasPattern= CodeUtils::HasPatternInString(tag,sourceFlags[j]);
				if(!hasPattern) continue;
				int sourceFlag= GetSourceFlag(sourceFlags[j]); 
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("sourceFlag="<<sourceFlag);
				#endif
				metadata.sourceFlag= sourceFlag;	
				metadata.hasSourceFlag= true;
			}

			//Parse source fit qualities
			for(size_t j=0;j<sourceFitQualities.size();j++){
				bool hasPattern= CodeUtils::HasPatternInString(tag,sourceFitQualities[j]);
				if(!hasPattern) continue;
				int sourceFitQuality= GetSourceFitQuality(sourceFitQualities[j]); 
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("sourceFitQuality="<<sourceFitQuality);
				#endif
				metadata.sourceFitQuality= sourceFitQuality;	
				metadata.hasSourceFitQuality= true;
			}
		}//close if has tag

	}//end loop metadata fields


	return 0;

}//close ParseRegionMetaData()

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
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
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
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
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
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
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
	double r= convFactor*atof(radiusStr.c_str());	
	region->cx= cx;
	region->cy= cy;
	region->r= r;
	
	return region;

}//close ParseCircleRegion()



DS9EllipseRegion* DS9RegionParser::ParseEllipseRegion(const std::string& data_str)
{
	//Check string 
	if(data_str=="") return nullptr;

	//Parse string and get coord data
	bool hasComment= CodeUtils::HasPatternInString(data_str,"#");
	std::string delim_start= "ellipse(";
	std::string delim_stop= ")";
	size_t first= data_str.find(delim_start);
	size_t last = data_str.find(delim_stop);
	if(first==std::string::npos){
		delim_start= "ellipse";//try with another region format (without parenthesis)
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
		DEBUG_LOG("ellipse parsed_data="<<parsed_data);
	#endif

	//Check if data are not given in real digit format (e.g. sexagesimal not supported)
	bool hasColon= CodeUtils::HasPatternInString(parsed_data,":");
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

	//Check radius units
	double convFactor= 1;//conversion factor to degrees
	bool hasUnits= CodeUtils::HasPatternInString(coords_str[2],"\"");
	std::string r1Str= coords_str[2];
	std::string r2Str= coords_str[3];
	if(hasUnits) {
		convFactor= 1./3600.;
		CodeUtils::RemovePatternInString(r1Str,"\"");
		CodeUtils::RemovePatternInString(r2Str,"\"");
		#ifdef LOGGING_ENABLED		
			DEBUG_LOG("DS9 ellipse region radius found in arcsec units, will convert to degree...");
		#endif
	}
	#ifdef LOGGING_ENABLED		
		DEBUG_LOG("DS9 ellipse region radius parsed: a="<<r1Str<<", b="<<r2Str);
	#endif
		
	//Extract coordinates and convert to double
	DS9EllipseRegion* region= new DS9EllipseRegion();
	double cx= atof(coords_str[0].c_str());	
	double cy= atof(coords_str[1].c_str());
	double r1= convFactor*atof(r1Str.c_str());
	double r2= convFactor*atof(r2Str.c_str());	
	double theta= atof(coords_str[4].c_str());
	region->cx= cx;
	region->cy= cy;
	region->a= r1;
	region->b= r2;
	region->theta= theta;
	
	return region;

}//close ParseEllipseRegion()


int DS9RegionParser::Read(std::vector<std::vector<std::string>>& data,std::vector<std::vector<std::string>>& metadata, int& wcsType,std::string filename)
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
	
	//std::vector<std::string> skipLinePatterns {"global","image","fk5","fk4","galactic"};
	std::vector<std::string> skipLinePatterns {"global","#"};
	std::vector<std::string> endDataReadPatterns {"#"};
	

	wcsType= eUNKNOWN_CS;
	bool foundWCS= false;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;


		//Check first character
		std::string field= "";
		istringstream parsedline_stream(parsedline);
		parsedline_stream>>field;
		char first_char= *(field.c_str());
		if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){//skip line
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("--> skip line (line="<<parsedline<<")");
			#endif
			continue;
		}

		//Search for pattern and skip if present
		bool skipLine= false;
		for(size_t i=0;i<skipLinePatterns.size();i++){
			bool hasPattern= CodeUtils::HasPatternInString(field,skipLinePatterns[i]);
			if(hasPattern) {
				skipLine= true;
				break;
			}
		}
		if(skipLine) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("--> skip line (line="<<parsedline<<")");
			#endif
			continue;
		}

		//Search for WCS line. If found skip line
		bool wcsLine= false;
		for (auto const& it : wcsTypeMap){
			bool hasPattern= CodeUtils::HasPatternInString(parsedline,it.first);
			if(hasPattern){
				wcsType= it.second;
				foundWCS= true;	
				wcsLine= true;	
			}
		}
		if(wcsLine) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("--> skip line (line="<<parsedline<<")");
			#endif
			continue;
		}

		//Parse data until '#'
		istringstream ss(parsedline);
		std::vector<std::string> fields;
		std::vector<std::string> fields_metadata;
		bool readMetadata= false;
		
		while(ss >> field)
    {			
			//Reached end data?
			char first_char= *(field.c_str());
			if(first_char=='#'){
				readMetadata= true;
				break;
			}
			if(first_char=='\n') {
				break;
			}

			//Save data
			fields.push_back(field);

		}//end loop data
	
		
		for(size_t k=0;k<fields.size();k++){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("--> fields["<<k<<"]="<<fields[k]);
			#endif
		}
		
		//Parse metadata
		if(readMetadata){	
			while(ss >> field)
    	{			
				//Reached end line
				char first_char= *(field.c_str());
				if(first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
					break;
				}

				//Save data
				fields_metadata.push_back(field);

			}//end loop data
		}//close if

		//Append data
		data.push_back(fields);
		metadata.push_back(fields_metadata);

		/*
		//Parse line fields and read data
		istringstream ss(parsedline);
		std::string field= "";
		//bool skipLine= false;
		bool dataEndReached= false;
		std::vector<std::string> fields;
		std::vector<std::string> fields_metadata;

		while(ss >> field)
    {
			
			//Skip first line
			//char first_char= *(field.c_str());
			//if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
			//	skipLine= true;
			//	break;
			//}
	
			//Search for WCS pattern
			//for (auto const& it : wcsTypeMap){
			//	bool hasPattern= CodeUtils::HasPatternInString(field,it.first);
			//	if(hasPattern){
			//		wcsType= it.second;
			//		foundWCS= true;
			//		break;
			//	}
			//}
			

			//Search for pattern and skip if present
			//for(size_t i=0;i<skipLinePatterns.size();i++){
			for(size_t i=0;i<endDataReadPatterns.size();i++){
				//bool hasPattern= CodeUtils::HasPatternInString(field,skipLinePatterns[i]);
				bool hasPattern= CodeUtils::HasPatternInString(field,endDataReadPatterns[i]);
				if(hasPattern) {
					//skipLine= true;
					dataEndReached= true;
					break;
				}
			}
			cout<<field<<"  (skip? "<<skipLine<<", dataEndReached? "<<dataEndReached<<") "<<endl;

			//if(skipLine) break;
			//fields.push_back(field);

			if(dataEndReached) {
				if(!skipLine) fields_metadata.push_back(field);
			}
			else{
				fields.push_back(field);
			}			

		}//end loop line
		
		//if(skipLine) continue;
		
		//Add region line to list
		if(!fields.empty()) data.push_back(fields);
		if(!fields_metadata.empty()) metadata.push_back(fields_metadata);
		*/

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
