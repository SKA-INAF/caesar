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
* @file SourceExporter.cc
* @class SourceExporter
* @brief SourceExporter class
*
* Class to export an image source in different formats
* @author S. Riggi
* @date 20/01/2015
*/


#include <SourceExporter.h>
#include <Source.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <AstroUtils.h>

#include <TObject.h>
#include <TEllipse.h>

#include <wcs.h>

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

ClassImp(Caesar::SourceExporter)

namespace Caesar {

SourceExporter::SourceExporter()
{

}

SourceExporter::~SourceExporter()
{

}


//=================================================
//==        DS9 EXPORTER
//=================================================
int SourceExporter::WriteToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS,int ds9WCSType,int ds9RegionFormat,WorldCoor* wcs)
{
	//## Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//## Saving DS9 file region
	std::string ds9WCSTypeHeader= "image";
	if(convertDS9RegionsToWCS) ds9WCSTypeHeader= AstroUtils::GetDS9WCSTypeHeader(ds9WCSType);

	DEBUG_LOG("Saving DS9 region header...");
	fprintf(fout,"global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"%s\n",ds9WCSTypeHeader.c_str());

	DEBUG_LOG("Saving "<<sources.size()<<" sources to file...");

	for(size_t k=0;k<sources.size();k++){
		int source_type= sources[k]->Type;
		bool isAtEdge= sources[k]->IsAtEdge();

		//If WCS is not computed, compute it
		if(convertDS9RegionsToWCS && !wcs){
			wcs= sources[k]->GetWCS(ds9WCSType);
			if(!wcs) WARN_LOG("Failed to compute WCS from source no "<<k<<"!");
		}
	
		//Get DS9 regions
		DEBUG_LOG("Dumping DS9 region info for source no. "<<k<<" ...");
		std::string regionInfo= "";
		if(ds9RegionFormat==ePolygonRegion) {
			regionInfo= SourceToDS9Region(sources[k],true,convertDS9RegionsToWCS,wcs,ds9WCSType);
		}
		else if(ds9RegionFormat==eEllipseRegion) {
			regionInfo= SourceToDS9EllipseRegion(sources[k],true);
		}
		else {
			WARN_LOG("Invalid DS9RegionType given ("<<ds9RegionFormat<<")");
			return -1;
		}

		//Write source region to file
		fprintf(fout,"%s\n",regionInfo.c_str());
	  	
	}//end loop sources
		
	DEBUG_LOG("Closing DS9 file region...");
	fclose(fout);

	return 0;

}//close WriteToDS9()


int SourceExporter::WriteComponentsToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS,int ds9WCSType,WorldCoor* wcs)
{
	//## Open file
	FILE* fout_fit= fopen(filename.c_str(),"w");

	std::string ds9WCSTypeHeader= "image";
	if(convertDS9RegionsToWCS) ds9WCSTypeHeader= AstroUtils::GetDS9WCSTypeHeader(ds9WCSType);

	//## Saving DS9 file region
	DEBUG_LOG("Saving DS9 region header for fitted source catalog...");
	fprintf(fout_fit,"global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout_fit,"%s\n",ds9WCSTypeHeader.c_str());

	DEBUG_LOG("Saving "<<sources.size()<<" sources to file...");
	bool useFWHM= true;

	for(size_t k=0;k<sources.size();k++){
		DEBUG_LOG("Dumping DS9 region fitting info for source no. "<<k<<" ...");

		//If WCS is not computed, compute it
		if(convertDS9RegionsToWCS && !wcs){
			wcs= sources[k]->GetWCS(ds9WCSType);
			if(!wcs) WARN_LOG("Failed to compute WCS from source no "<<k<<"!");
		}

		//Get DS9 regions for fitted components
		std::string regionInfo= SourceToDS9FittedEllipseRegion(sources[k],useFWHM,true,convertDS9RegionsToWCS,wcs,ds9WCSType);

		fprintf(fout_fit,"%s\n",regionInfo.c_str());
	}//end loop sources
		
	DEBUG_LOG("Closing DS9 file region for fitted sources...");
	fclose(fout_fit);

	return 0;

}//close WriteComponentsToDS9()

std::string SourceExporter::GetDS9RegionColor(Source* source)
{
	std::string colorStr= "white";
	if(source->Type==Source::eExtended) colorStr= "green";
	else if(source->Type==Source::eCompactPlusExtended) colorStr= "magenta";
	else if(source->Type==Source::ePointLike) colorStr= "red";
	else if(source->Type==Source::eCompact) colorStr= "blue";
	else colorStr= "white";
			
	return colorStr;
		
}//close GetDS9RegionColor()

const std::string SourceExporter::SourceToDS9Region(Source* source,bool dumpNestedSourceInfo,bool convertToWCS,WorldCoor* wcs,int coordSystem)
{
	//Check source
	if(!source){
		WARN_LOG("Null input source ptr given!");
		return std::string("");
	}

	//Check if has pixels
	//NB: DS9 crashes miserably when given a polygon region with one point 
	if(source->NPix<=1) return std::string("");

	//Convert contours to WCS?
	std::stringstream sstream;
	std::string regionText= source->GetName();
	std::string regionColor= GetDS9RegionColor(source);
	std::vector<std::string> regionTags {source->GetDS9RegionTag()};
	std::string region= "";
	bool useImageCoords= true;
	if(convertToWCS) useImageCoords= false;
	
	if(convertToWCS){
		int pixOffset= 1;
		std::vector<Contour*> contours_wcs= source->GetWCSContours(wcs,coordSystem,pixOffset);		
		if(contours_wcs.empty()){
			WARN_LOG("Failed to convert contours in WCS, region will be empty!");
		}
		else{
			//Loop over WCS contours
			for(size_t i=0; i<contours_wcs.size(); i++){ 
				region= AstroUtils::ContourToDS9Region(contours_wcs[i],regionText,regionColor,regionTags,useImageCoords);
			}

			//Delete contours at the end
			CodeUtils::DeletePtrCollection<Contour>(contours_wcs);
		}

	}//close if
	else{
		std::vector<Contour*> contours= source->GetContours();
		for(size_t i=0; i<contours.size(); i++){ 
			region= AstroUtils::ContourToDS9Region(contours[i],regionText,regionColor,regionTags,useImageCoords);
		}
	}//close else

	sstream<<region;

	//###### FILL NESTED SOURCE REGIONS ###########
	//Fill nested source regions
	bool hasNestedSources= source->HasNestedSources();
	if(dumpNestedSourceInfo && hasNestedSources){
		std::vector<Source*> nestedSources= source->GetNestedSources();

		sstream<<endl;
		for(size_t k=0;k<nestedSources.size();k++)
		{
			std::string regionText_nested= nestedSources[k]->GetName();
			std::string regionColor_nested= nestedSources[k]->GetDS9RegionColor();
			std::vector<std::string> regionTags_nested {nestedSources[k]->GetDS9RegionTag()};
			std::string region_nested= "";
			if(convertToWCS){
				int pixOffset= 1;
				std::vector<Contour*> contours_wcs= nestedSources[k]->GetWCSContours(wcs,coordSystem,pixOffset);
				if(contours_wcs.empty()){
					WARN_LOG("Failed to convert contours in WCS, region will be empty!");
				}
				else{
					//Loop over WCS contours
					for(size_t i=0;i<contours_wcs.size(); i++){ 
						region_nested= AstroUtils::ContourToDS9Region(contours_wcs[i],regionText_nested,regionColor_nested,regionTags_nested,useImageCoords);
					}

					//Delete contours at the end
					CodeUtils::DeletePtrCollection<Contour>(contours_wcs);
				}
			}//close if
			else{
				std::vector<Contour*> nestedContours= nestedSources[k]->GetContours();
				for(size_t i=0;i<nestedContours.size(); i++){ 
					region_nested= AstroUtils::ContourToDS9Region(nestedContours[i],regionText_nested,regionColor_nested,regionTags_nested,useImageCoords);
				}
			}//close else

			sstream<<region_nested;
			if(k!=nestedSources.size()-1) sstream<<endl;

		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close SourceToDS9Region()


const std::string SourceExporter::SourceToDS9EllipseRegion(Source* source,bool dumpNestedSourceInfo)
{
	//Check source
	if(!source){
		WARN_LOG("Null input source ptr given!");
		return std::string("");
	}

	//Get contours
	std::vector<Contour*> contours= source->GetContours();
			
	//Ellipse x y radius radius angle
	std::stringstream sstream;
	sstream<<"ellipse ";
	for(size_t i=0; i<contours.size(); i++){ 
		if(!contours[i]->HasEllipseFit) continue;
		double EllX= contours[i]->EllipseCenter.X();
		double EllY= contours[i]->EllipseCenter.Y();
		double EllMajAxis= contours[i]->EllipseMajAxis;
		double EllMinAxis= contours[i]->EllipseMinAxis;
		double EllRotAxis= contours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
		sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
	}
	sstream<<"# text={S"<<source->Id<<"}";

	//Nested sources
	bool hasNestedSources= source->HasNestedSources();
	if(dumpNestedSourceInfo && hasNestedSources){		
		std::vector<Source*> nestedSources= source->GetNestedSources();
	
		sstream<<endl;
		for(size_t k=0;k<nestedSources.size();k++){
			std::vector<Contour*> nestedContours= nestedSources[k]->GetContours();

			sstream<<"ellipse ";
			
			for(unsigned int i=0; i<nestedContours.size(); i++){ 
				if(!nestedContours[i]->HasEllipseFit) continue;
				double EllX= nestedContours[i]->EllipseCenter.X();
				double EllY= nestedContours[i]->EllipseCenter.Y();
				double EllMajAxis= nestedContours[i]->EllipseMajAxis;
				double EllMinAxis= nestedContours[i]->EllipseMinAxis;
				double EllRotAxis= nestedContours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
				sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
			}//end loop contours
				
			sstream<<"# text={"<<source->GetName()<<"_Nest"<<k<<"}";
			if(k!=nestedSources.size()-1) sstream<<endl;
		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close SourceToDS9EllipseRegion()


const std::string SourceExporter::SourceToDS9FittedEllipseRegion(Source* source,bool useFWHM,bool dumpNestedSourceInfo,bool convertToWCS,WorldCoor* wcs,int coordSystem)
{
	//Check source
	if(!source){
		WARN_LOG("Null input source ptr given!");
		return std::string("");
	}

	//Check WCS & metadata
	ImgMetaData* metadata= source->GetImageMetaData();
	if(convertToWCS && !wcs && !metadata){
		WARN_LOG("Requested to convert ellipse to WCS but no wcs was provided and no metadata to compute it are available!");
		return std::string("");
	}

	bool useImageCoords= true;
	int pixOffset= 0;
	if(convertToWCS) {
		useImageCoords= false;
		pixOffset= 1;
	}

	//Check if source has fit info
	std::stringstream sstream;
	bool hasFitInfo= source->HasFitInfo();
	if(hasFitInfo){
		//Get fit ellipses
		std::vector<TEllipse*> ellipses;
		
		if(source->GetFitEllipses(ellipses,useFWHM,convertToWCS,wcs,coordSystem,pixOffset)<0){
			ERROR_LOG("Failed to get WorldCoord system from metadata!");
			return std::string("");
		}
	
		//Loop over fit ellipses and convert to DS9 regions
		for(size_t i=0;i<ellipses.size();i++){
			if(!ellipses[i]) continue;

			//Get encoded string region
			std::string regionText(Form("%s_fitcomp%d",source->GetName(),(int)(i+1)));
			std::string regionColor= "red";
			std::vector<std::string> regionTags {"point-like","fitted component"};
			std::string region= AstroUtils::EllipseToDS9Region(ellipses[i],regionText,regionColor,regionTags,useImageCoords);
			sstream<<region;

			if(i!=ellipses.size()-1) sstream<<endl;
		}//end loop ellipses

		//Delete ellipses
		CodeUtils::DeletePtrCollection<TEllipse>(ellipses);

	}//close if has fit info

	//Loop over nested components and get fit ellipse regions
	bool hasNestedSources= source->HasNestedSources();
	if(dumpNestedSourceInfo && hasNestedSources){	
		std::vector<Source*> nestedSources= source->GetNestedSources();
	
		for(size_t k=0;k<nestedSources.size();k++){	
			std::string nestedRegionStr= SourceToDS9FittedEllipseRegion(nestedSources[k],useFWHM,false,convertToWCS,wcs,coordSystem);
			if(nestedRegionStr!="") sstream<<nestedRegionStr<<endl;
		}//end loop nested sources
	}//close if

	return sstream.str();
	
}//close SourceToDS9FittedEllipseRegion()


}//close namespace

