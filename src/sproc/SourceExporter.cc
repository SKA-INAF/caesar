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
#include <WCSUtils.h>
#include <Consts.h>

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

ClassImp(Caesar::SourceExporter)

namespace Caesar {

SourceExporter::SourceExporter()
{

}

SourceExporter::~SourceExporter()
{

}

//=================================================
//==        ASCII EXPORTER
//=================================================
int SourceExporter::WriteToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo,bool convertBrightnessToFlux,char delimiter)
{
	//Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Write header
	std::stringstream ss_coldescr;

	std::stringstream ss;
	ss<<"#======== CAESAR RESULTS ======\n";
	ss<<"# This ascii file contains source islands (or blobs) found by CAESAR.\n";
	ss<<"# Source parameters are described in the header below.\n";
	ss<<"# NB: Fitted components, if available, are reported in a different ascii file. \n";
	ss<<"#     We include here the total fitted flux density summed over all components.\n";
	ss<<"#\n";
	ss<<"# ---- HEADER ----\n";
	ss<<"# name - Source name assigned by finder\n";
	ss<<"# iauName - Source name in IAU notation\n";
	ss<<"# npix - Number of pixels in source\n";
	ss<<"# nComponents - Number of fitted components (=0 if fit not performed or failed)\n";
	ss<<"# nNestedSources - Number of nested sources detected\n";
	ss<<"# x - Source centroid in image coordinates along x axis \n";
	ss<<"# y - Source centroid in image coordinates along y axis \n";
	ss<<"# x_w - Source centroid in image coordinates along x axis, weighted by pixel fluxes \n";
	ss<<"# y_w - Source centroid in image coordinates along y axis, weighted by pixel fluxes \n";
	ss<<"# x_wcs - Source centroid in selected WCS coordinate (deg) along x axis\n";
	ss<<"# y_wcs - Source centroid in selected WCS coordinate (deg) along y axis\n";
	ss<<"# x_w_wcs - Source centroid in selected selected WCS coordinate (deg) along x axis, weighted by pixel fluxes\n";
	ss<<"# y_w_wcs - Source centroid in selected selected WCS coordinate (deg) along y axis, weighted by pixel fluxes\n";
	ss<<"# xmin - Source minimum pixel image coordinate along x axis\n";
	ss<<"# xmax - Source maximum pixel image coordinate along x axis\n";
	ss<<"# ymin - Source minimum pixel image coordinate along y axis\n";
	ss<<"# ymax - Source maximum pixel image coordinate along y axis\n";
	ss<<"# xmin_wcs - Source minimum pixel WCS coordinate (deg) along x axis\n";
	ss<<"# xmax_wcs - Source maximum pixel WCS coordinate (deg) along x axis\n";
	ss<<"# ymin_wcs - Source minimum pixel WCS coordinate (deg) along y axis\n";
	ss<<"# ymax_wcs - Source maximum pixel WCS coordinate (deg) along y axis\n";
	ss<<"# nu - Spectral axis value present in image header. If frequency it is given in GHz units.\n";
	ss<<"# Stot - Sum of source pixel brightness in Jy/beam units.\n";
	ss<<"# Smax - Maximum source pixel flux in Jy/beam units.\n";
	if(convertBrightnessToFlux){
		ss<<"# S - Fitted source flux density in Jy units\n";
		ss<<"# SErr - Fitted source flux density error in Jy units\n";
	}
	else{
		ss<<"# S - Fitted source brightness in Jy/beam units (divide by beam area to obtain flux density).\n";
		ss<<"# SErr - Fitted source brightness error in Jy/beam units (divide by beam area to obtain flux density).\n";
	}
	ss<<"# beamArea - Number of pixels in beam. Used to convert flux parameters from Jy/beam to Jy/pixel (e.g. Jy/pixel=Jy/beam/beamarea).\n";
	ss<<"# bkgSum - Background estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# rmsSum - Noise (rms) estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# sourceType - Source type tag (eUnknownType=0,eCompact=1,ePointLike=2,eExtended=3,eCompactPlusExtended=4)\n";
	ss<<"# sourceFlag - Source flag (eReal=1,eCandidate=2,eFake=3)\n";
	ss<<"# isGoodSource - Bool flag indicating if source was tagged as good (true) or bad (false) in the finding process\n";
	ss<<"# sourceNestedLevel - Source nested level (0=mother source,1=nested source,...)\n";

	ss_coldescr<<"# ";
	ss_coldescr<<"name"<<delimiter;
	ss_coldescr<<"iauName"<<delimiter;
	ss_coldescr<<"npix"<<delimiter;
	ss_coldescr<<"nComponents"<<delimiter;
	ss_coldescr<<"nNestedSources"<<delimiter;
	ss_coldescr<<"x"<<delimiter;
	ss_coldescr<<"y"<<delimiter;
	ss_coldescr<<"x_w"<<delimiter;
	ss_coldescr<<"y_w"<<delimiter;
	ss_coldescr<<"x_wcs(deg)"<<delimiter;
	ss_coldescr<<"y_wcs(deg)"<<delimiter;
	ss_coldescr<<"x_w_wcs(deg)"<<delimiter;
	ss_coldescr<<"y_w_wcs(deg)"<<delimiter;
	ss_coldescr<<"xmin"<<delimiter;
	ss_coldescr<<"xmax"<<delimiter;
	ss_coldescr<<"ymin"<<delimiter;
	ss_coldescr<<"ymax"<<delimiter;
	ss_coldescr<<"xmin_wcs(deg)"<<delimiter;
	ss_coldescr<<"xmax_wcs(deg)"<<delimiter;
	ss_coldescr<<"ymin_wcs(deg)"<<delimiter;
	ss_coldescr<<"ymax_wcs(deg)"<<delimiter;
	ss_coldescr<<"nu(GHz)"<<delimiter;
	ss_coldescr<<"Stot(Jy/beam)"<<delimiter;
	ss_coldescr<<"Smax(Jy/beam)"<<delimiter;
	if(convertBrightnessToFlux){
		ss_coldescr<<"S(Jy)"<<delimiter;
		ss_coldescr<<"Serr(Jy)"<<delimiter;
	}
	else{
		ss_coldescr<<"S(Jy/beam)"<<delimiter;
		ss_coldescr<<"Serr(Jy/beam)"<<delimiter;
	}
	ss_coldescr<<"beamArea"<<delimiter;
	ss_coldescr<<"bkgSum(Jy/beam)"<<delimiter;
	ss_coldescr<<"rmsSum(Jy/beam)"<<delimiter;
	ss_coldescr<<"sourceType"<<delimiter;
	ss_coldescr<<"sourceFlag"<<delimiter;
	ss_coldescr<<"isGoodSource"<<delimiter;
	ss_coldescr<<"sourceNestedLevel"<<delimiter;

	if(writeAdditionalSourceInfo)
	{
		// - Spectral index data
		ss<<"# hasSpectralIndexData - Flag indicating if spectral index data are available\n";
		ss<<"# isMultiSourceMatchIndex - Flag indicating if spectral index data was computed by summing up component fluxes\n";
		ss<<"# spectralIndex - Spectral index value\n";
		ss<<"# spectralIndexErr - Spectral index value error\n";
		ss<<"# isSpectralIndexFit - Flag indicating if spectral index was fitted from a SED model\n";
		ss<<"# spectralFitChi2 - Spectral index fit chi2\n";
		ss<<"# spectralFitNDF - Spectral index fit ndf\n";

		// - Astro object ids
		//ss<<"# objLocationId - Object location id (0=UNKNOWN,1=GALACTIC,2=EXTRAGALACTIC)\n";
		ss<<"# objClassStrId - Object class string id\n";
		ss<<"# objClassId - Object class id\n";
		ss<<"# objClassSubId - Object class subid\n";
		ss<<"# objConfirmed - Object confirmed\n";

		// - Astro crossmatched object names
		ss<<"# hasMatchedObject - Bool flag indicating if source has matched astro object\n";
		ss<<"# matchedObjNames - Cross-matched astronomical object names\n";

		ss_coldescr<<"hasSpectralIndexData"<<delimiter;
		ss_coldescr<<"isMultiSourceMatchIndex"<<delimiter;
		ss_coldescr<<"spectralIndex"<<delimiter;
		ss_coldescr<<"spectralIndexErr"<<delimiter;
		ss_coldescr<<"isSpectralIndexFit"<<delimiter;
		ss_coldescr<<"spectralFitChi2"<<delimiter;
		ss_coldescr<<"spectralFitNDF"<<delimiter;
		//ss_coldescr<<"objLocationId"<<delimiter;
		ss_coldescr<<"objClassStrId"<<delimiter;
		ss_coldescr<<"objClassId"<<delimiter;
		ss_coldescr<<"objClassSubId"<<delimiter;
		ss_coldescr<<"objConfirmed"<<delimiter;
		ss_coldescr<<"hasMatchedObject"<<delimiter;
		ss_coldescr<<"matchedObjNames";
	
	}//close if

	ss<<"# ---------------\n";
	ss<<"#\n";
	ss<<ss_coldescr.str()<<endl;

	fprintf(fout,"%s",ss.str().c_str());

	//Loop sources
	bool deleteWCS= false;

	for(size_t k=0;k<sources.size();k++){
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			#ifdef LOGGING_ENABLED
			 INFO_LOG("WCS not given, retrieving it from the first source ...");
			#endif
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No metadata are available to retrieve WCS!");
				#endif
				return -1;
			}
			wcs= metadata->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WCS from metadata!");
				#endif
				return -1;
			}
			deleteWCS= true;
		}

		//Print source info to ascii
		std::vector<std::string> slist= SourceToAscii(sources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightnessToFlux,delimiter);
		for(size_t j=0;j<slist.size();j++){
			fprintf(fout,"%s\n",slist[j].c_str());
		}
		
	}//end loop sources

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	//Close file
	fclose(fout);

	return 0;

}//close WriteToAscii()


const std::vector<std::string> SourceExporter::SourceToAscii(Source* source,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo,bool convertBrightnessToFlux,char delimiter)
{
	//Init string list
	std::vector<std::string> sourceStrList;

	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return sourceStrList;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWCS(wcsType);
			if(wcs) deleteWCS= true;
			else {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get WCS from metadata!");
				#endif
			}
		}
		else {
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
		}		
	}
	

	//Compute IAU name
	WCS* wcs_j2000= metadata->GetWCS(eJ2000);
	if(!wcs_j2000){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get J2000 WCS from metadata!");
		#endif
	}
	bool useWeightedPos= false;
	std::string iauName= source->GetName();
	//if(wcs) iauName= source->GetIAUName(useWeightedPos,wcs,wcsType);
	if(wcs_j2000) iauName= source->GetIAUName(useWeightedPos,wcs_j2000,eJ2000);

	//Compute WCS centroid
	double X0_wcs= 0;
	double Y0_wcs= 0;
	double X0_weighted_wcs= 0;
	double Y0_weighted_wcs= 0;
	if(wcs){
		source->GetWCSPos(X0_wcs,Y0_wcs,wcs,wcsType);
		source->GetWCSWeightedPos(X0_weighted_wcs,Y0_weighted_wcs,wcs,wcsType);
	}

	//Get spectral axis info
	double Nu= -999;
	double dNu= -999;
	std::string units= "";
	source->GetSpectralAxisInfo(Nu,dNu,units);
	CodeUtils::StripBlankSpaces(units);
	if(metadata && units=="Hz") {//convert to GHz
		Nu/= 1.e+9;
		dNu/= 1.e+9;
	}
	
	
	double beamArea= source->GetBeamFluxIntegral();
	double fluxDensity= 0;
	double fluxDensityErr= 0;
	double flux= 0;
	double fluxErr= 0;
	if(source->HasFitInfo()){
		source->GetFluxDensity(fluxDensity);
		source->GetFluxDensityErr(fluxDensityErr);
		flux= fluxDensity/beamArea;
		fluxErr= fluxDensityErr/beamArea;
	}
	
	//Get WCS bounding box
	float xmin, xmax, ymin, ymax;
	source->GetSourceRange(xmin,xmax,ymin,ymax);
	
	double xmin_wcs, xmax_wcs, ymin_wcs, ymax_wcs;
	source->GetWCSSourceRange(xmin_wcs,xmax_wcs,ymin_wcs,ymax_wcs,wcs,wcsType);
		
	//## Set fields
	std::stringstream ss;
	
	//- Source name
	ss<<source->GetName()<<delimiter;
	ss<<iauName<<delimiter;

	
	//- Number of pixels
	ss<<source->NPix<<delimiter;

	//- Number of sub-components
	ss<<source->GetNSelFitComponents()<<delimiter;
	ss<<source->GetNestedSourceNumber()<<delimiter;

	//- Pixel centroids
	ss<<source->X0<<"\t"<<source->Y0<<delimiter;
	ss<<source->GetSx()<<"\t"<<source->GetSy()<<delimiter;
	ss<<setprecision(8)<<X0_wcs<<delimiter<<Y0_wcs<<delimiter<<X0_weighted_wcs<<delimiter<<Y0_weighted_wcs<<delimiter;

	//- Bounding box
	ss<<xmin<<delimiter<<xmax<<delimiter<<ymin<<delimiter<<ymax<<delimiter;
	ss<<xmin_wcs<<delimiter<<xmax_wcs<<delimiter<<ymin_wcs<<delimiter<<ymax_wcs<<delimiter;
	
	//- Frequency
	ss<<Nu<<delimiter;

	//- Flux 
	ss<<source->GetS()<<delimiter;
	ss<<source->GetSmax()<<delimiter;

	
	if(convertBrightnessToFlux){
		ss<<flux<<delimiter;
		ss<<fluxErr<<delimiter;
	}
	else{
		ss<<fluxDensity<<delimiter;
		ss<<fluxDensityErr<<delimiter;
	}
	ss<<beamArea<<delimiter;

	//- Bkg/noise estimators
	ss<<source->GetBkgSum()<<delimiter;
	ss<<source->GetBkgRMSSum()<<delimiter;

	//- Source flags
	//ss<<source->Type<<"\t"<<source->Flag<<"\t"<<source->IsGoodSource()<<"\t";
	ss<<source->Type<<delimiter<<source->Flag<<delimiter<<source->IsGoodSource()<<delimiter;
	ss<<source->GetDepthLevel();

	//## Additional source info
	if(writeAdditionalSourceInfo)
	{	
		ss<<delimiter;

		//- Spectral index data
		SpectralIndexData sid= source->GetSpectralIndexData();	
		bool hasSpectralIndexData= ((source->HasSpectralIndexData()) && (sid.hasSpectralIndex));
			
		ss<<hasSpectralIndexData<<delimiter;
		ss<<sid.isMultiSourceMatchIndex<<delimiter;
		ss<<sid.spectralIndex<<delimiter;
		ss<<sid.spectralIndexErr<<delimiter;
		ss<<sid.isSpectralIndexFit<<delimiter;
		ss<<sid.spectralFitChi2<<delimiter;
		ss<<sid.spectralFitNDF<<delimiter;
	
		//- Astro object IDs
		//ss<<source->ObjLocationId<<delimiter;
		ss<<"\""<<source->ObjClassStrId<<"\""<<delimiter;
		ss<<source->ObjClassId<<delimiter;
		ss<<source->ObjClassSubId<<delimiter;
		ss<<source->ObjConfirmed<<delimiter;

		//- Astro object cross matches
		bool hasAstroObjs= source->HasAstroObjects();
		std::vector<AstroObject> astroObjs= source->GetAstroObjects();
		
		ss<<hasAstroObjs<<delimiter;

		if(hasAstroObjs && !astroObjs.empty()){
			std::stringstream sstream;
			if(astroObjs.size()==1){
				sstream<<astroObjs[0].name;
			}
			else{
				for(size_t j=0;j<astroObjs.size()-1;j++){
					sstream<<astroObjs[j].name<<";";
				}
				sstream<<astroObjs[astroObjs.size()-1].name;
			}
			
			//ss<<sstream.str();
			ss<<"\""<<sstream.str()<<"\"";
	
		}//close if
		else{
			//ss<<"XXX";
			ss<<"\"XXX\"";
		}
		
	}//close if writeAdditionalSourceInfo
	
	//Add string to list
	sourceStrList.push_back(ss.str());

	//Store nested sources
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++)
		{
			std::vector<std::string> sourceStrList_nested= SourceToAscii(nestedSources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightnessToFlux,delimiter);
			if(!sourceStrList_nested.empty()){
				sourceStrList.insert(sourceStrList.end(),sourceStrList_nested.begin(),sourceStrList_nested.end());
			}
		}
	}//close if has nested sources 
		

	//Delete WCS
	if(deleteWCS) {
		WCSUtils::DeleteWCS(&wcs);
		WCSUtils::DeleteWCS(&wcs_j2000);
	}

	return sourceStrList;

}//close SourceToAscii()



int SourceExporter::WriteComponentsToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo,bool convertBrightnessToFlux,char delimiter)
{
	//Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Write header
	std::stringstream ss_coldescr;

	std::stringstream ss;
	ss<<"#======== CAESAR RESULTS ======\n";
	ss<<"# This ascii file contains source components fitted to islands/blobs detected by CAESAR.\n";
	ss<<"# Source parameters are described in the header below.\n";
	ss<<"#\n";
	ss<<"# ---- HEADER ----\n";
	ss<<"# name - Island source name assigned by finder\n";
	ss<<"# npix - Number of pixels in source\n";
	ss<<"# componentId - Fitted component id\n";
	ss<<"# iauName - Fitted component name in IAU notation\n";
	ss<<"# x - Fitted component centroid in image coordinates along x axis \n";
	ss<<"# y - Fitted component centroid in image coordinates along y axis \n";
	ss<<"# x_err - Fitted component centroid error in image coordinates along x axis \n";
	ss<<"# y_err - Fitted component centroid error in image coordinates along y axis \n";
	ss<<"# x_wcs - Fitted component centroid in selected WCS coordinate (deg) along x axis\n";
	ss<<"# y_wcs - Fitted component centroid in selected WCS coordinate (deg) along y axis\n";
	ss<<"# x_wcs_err - Fitted component centroid error in selected WCS coordinate (deg) along x axis\n";
	ss<<"# y_wcs_err - Fitted component centroid error in selected WCS coordinate (deg) along y axis\n";
	ss<<"# nu - Spectral axis value present in image header. If frequency it is given in GHz units.\n";
	ss<<"# Speak - Fitted component peak brightness in Jy/beam units.\n";
	ss<<"# Speak_err - Fitted component peak brightness error in Jy/beam units.\n";
	if(convertBrightnessToFlux){
		ss<<"# S - Fitted component flux density in Jy units\n";
		ss<<"# S_err - Fitted component flux density error in Jy units\n";
		ss<<"# S_island - Fitted island flux density in Jy units\n";
		ss<<"# S_island_err - Fitted island flux density error in Jy units\n";
	}
	else{
		ss<<"# S - Fitted component brightness in Jy/beam units (divide by beam area to obtain flux density).\n";
		ss<<"# S_err - Fitted component brightness error in Jy/beam units (divide by beam area to obtain flux density).\n";
		ss<<"# S_island - Fitted island brightness in Jy/beam units (divide by beam area to obtain flux density).\n";
		ss<<"# S_island_err - Fitted island brightness error in Jy/beam units (divide by beam area to obtain flux density).\n";
	}
	ss<<"# beamArea - Number of pixels in beam. Used to convert flux parameters from Jy/beam to Jy/pixel (e.g. Jy/pixel=Jy/beam/beamarea).\n";
	ss<<"# bmaj - Fitted component ellipse major axis in image coordinates\n";
	ss<<"# bmin - Fitted component ellipse major axis in image coordinates\n";
	ss<<"# pa - Fitted component ellipse position angles in deg (measured counterclock-wise from North)\n";
	ss<<"# bmaj_err - Fitted component ellipse major axis error in image coordinates\n";
	ss<<"# bmin_err - Fitted component ellipse major axis error in image coordinates\n";
	ss<<"# pa_err - Fitted component ellipse position angles error in deg\n";
	ss<<"# bmaj_wcs - Fitted component ellipse major axis in selected WCS coordinates (arcsec)\n";
	ss<<"# bmin_wcs - Fitted component ellipse major axis in selected WCS coordinates (arcsec)\n";
	ss<<"# pa_wcs - Fitted component ellipse position angles in deg (measured counterclock-wise from North)\n";
	ss<<"# bmaj_wcs_err - Fitted component ellipse major axis error in selected WCS coordinates (arcsec)\n";
	ss<<"# bmin_wcs_err - Fitted component ellipse major axis error in selected WCS coordinates (arcsec)\n";
	ss<<"# pa_wcs_err - Fitted component ellipse position angles error in deg\n";
	ss<<"# bmaj_beam - Beam ellipse major axis (arcsec)\n";
	ss<<"# bmin_beam - Beam ellipse minor axis (arcsec)\n";
	ss<<"# pa_beam - Beam ellipse position angles in deg (measured counterclock-wise from North)\n";
	ss<<"# bmaj_deconv_wcs - Fitted component ellipse major axis in selected WCS coordinates, deconvolved by beam (arcsec)\n";
	ss<<"# bmin_deconv_wcs - Fitted component ellipse major axis in selected WCS coordinates, deconvolved by beam (arcsec)\n";
	ss<<"# pa_deconv_wcs - Fitted component ellipse position angles in deg, deconvolved by beam (measured counterclock-wise from North)\n";
	ss<<"# fitBeamEllipseEccentricityRatio - Ratio between eccentricities of fitted and beam ellipses\n";
	ss<<"# fitBeamEllipseAreaRatio - Ratio between areas of fitted ellipse and beam ellipse \n";
	ss<<"# fitBeamEllipseRotAngle - Rotation angle in degrees (range 0-180) between fit ellipse and beam ellipse \n";
	ss<<"# bkgSum - Background estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# rmsSum - Noise (rms) estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# chi2 - Fit chisquare.\n";
	ss<<"# ndf - Fit number of degrees of freedom.\n";
	ss<<"# fitQuality - Fit quality flag (eBadFit=0,eLQFit=1,eMQFit=2,eHQFit=3)\n";
	ss<<"# fitComponentFlag- Fitted component flag (eReal=1,eCandidate=2,eFake=3)\n";
	ss<<"# fitComponentType- Fitted component type (eUnknown=0,eCompact=1,ePoint-Like=2,eExtended=3)\n";

	ss_coldescr<<"# ";
	ss_coldescr<<"name"<<delimiter;
	ss_coldescr<<"npix"<<delimiter;
	ss_coldescr<<"componentId"<<delimiter;
	ss_coldescr<<"iauName"<<delimiter;
	ss_coldescr<<"x"<<delimiter;
	ss_coldescr<<"y"<<delimiter;
	ss_coldescr<<"x_err"<<delimiter;
	ss_coldescr<<"y_err"<<delimiter;
	ss_coldescr<<"x_wcs(deg)"<<delimiter;
	ss_coldescr<<"y_wcs(deg)"<<delimiter;
	ss_coldescr<<"x_wcs_err(deg)"<<delimiter;
	ss_coldescr<<"y_wcs_err(deg)"<<delimiter;
	ss_coldescr<<"nu(GHz)"<<delimiter;
	ss_coldescr<<"Speak(Jy/beam)"<<delimiter;
	ss_coldescr<<"Speak_err(Jy/beam)"<<delimiter;
	if(convertBrightnessToFlux){
		ss_coldescr<<"S(Jy)"<<delimiter;
		ss_coldescr<<"S_err(Jy)"<<delimiter;
		ss_coldescr<<"S_island(Jy)"<<delimiter;
		ss_coldescr<<"S_island_err(Jy)"<<delimiter;
	}
	else{
		ss_coldescr<<"S(Jy/beam)"<<delimiter;
		ss_coldescr<<"S_err(Jy/beam)"<<delimiter;
		ss_coldescr<<"S_island(Jy/beam)"<<delimiter;
		ss_coldescr<<"S_island_err(Jy/beam)"<<delimiter;
	}
	ss_coldescr<<"beamArea"<<delimiter;
	ss_coldescr<<"bmaj"<<delimiter;
	ss_coldescr<<"bmin"<<delimiter;
	ss_coldescr<<"pa"<<delimiter;
	ss_coldescr<<"bmaj_err"<<delimiter;
	ss_coldescr<<"bmin_err"<<delimiter;
	ss_coldescr<<"pa_err"<<delimiter;
	ss_coldescr<<"bmaj_wcs(arcsec)"<<delimiter;
	ss_coldescr<<"bmin_wcs(arcsec)"<<delimiter;
	ss_coldescr<<"pa_wcs(deg)"<<delimiter;
	ss_coldescr<<"bmaj_wcs_err(arcsec)"<<delimiter;
	ss_coldescr<<"bmin_wcs_err(arcsec)"<<delimiter;
	ss_coldescr<<"pa_wcs_err(deg)"<<delimiter;
	ss_coldescr<<"bmaj_beam(arcsec)"<<delimiter;
	ss_coldescr<<"bmin_beam(arcsec)"<<delimiter;
	ss_coldescr<<"pa_beam(deg)"<<delimiter;
	ss_coldescr<<"bmaj_deconv_wcs(arcsec)"<<delimiter;
	ss_coldescr<<"bmin_deconv_wcs(arcsec)"<<delimiter;
	ss_coldescr<<"pa_deconv_wcs(deg)"<<delimiter;
	ss_coldescr<<"fitBeamEllipseEccentricityRatio"<<delimiter;
	ss_coldescr<<"fitBeamEllipseAreaRatio"<<delimiter;
	ss_coldescr<<"fitBeamEllipseRotAngle"<<delimiter;
	ss_coldescr<<"bkgSum(Jy/beam)"<<delimiter;
	ss_coldescr<<"rmsSum(Jy/beam)"<<delimiter;
	ss_coldescr<<"chi2"<<delimiter;
	ss_coldescr<<"ndf"<<delimiter;
	ss_coldescr<<"fitQuality"<<delimiter;
	ss_coldescr<<"fitComponentFlag"<<delimiter;
	ss_coldescr<<"fitComponentType"<<delimiter;

	if(writeAdditionalSourceInfo)
	{
		// - Spectral index data
		ss<<"# hasSpectralIndexData - Flag indicating if spectral index data are available\n";
		ss<<"# isMultiSourceMatchIndex - Flag indicating if spectral index data was computed by summing up component fluxes\n";
		ss<<"# spectralIndex - Spectral index value\n";
		ss<<"# spectralIndexErr - Spectral index value error\n";
		ss<<"# isSpectralIndexFit - Flag indicating if spectral index was fitted from a SED model\n";
		ss<<"# spectralFitChi2 - Spectral index fit chi2\n";
		ss<<"# spectralFitNDF - Spectral index fit ndf\n";

		// - Astro object ids
		//ss<<"# objLocationId - Object location id (0=UNKNOWN,1=GALACTIC,2=EXTRAGALACTIC)\n";
		ss<<"# objClassStrId - Object class string id\n";
		ss<<"# objClassId - Object class id\n";
		ss<<"# objClassSubId - Object class subid\n";
		ss<<"# objClassConfirmed - Object confirmed\n";

		// - Astro crossmatched object names
		ss<<"# hasMatchedObject - Bool flag indicating if source component has matched astro object\n";
		ss<<"# matchedObjNames - Cross-matched astronomical object names\n";

		ss_coldescr<<"hasSpectralIndexData"<<delimiter;
		ss_coldescr<<"isMultiSourceMatchIndex"<<delimiter;
		ss_coldescr<<"spectralIndex"<<delimiter;
		ss_coldescr<<"spectralIndexErr"<<delimiter;
		ss_coldescr<<"isSpectralIndexFit"<<delimiter;
		ss_coldescr<<"spectralFitChi2"<<delimiter;
		ss_coldescr<<"spectralFitNDF"<<delimiter;
		ss_coldescr<<"objClassStrId"<<delimiter;
		ss_coldescr<<"objClassId"<<delimiter;
		ss_coldescr<<"objClassSubId"<<delimiter;
		ss_coldescr<<"objClassConfirmed"<<delimiter;
		ss_coldescr<<"hasMatchedObject"<<delimiter;
		ss_coldescr<<"matchedObjNames";

	}//close if

	ss<<"# ---------------\n";
	ss<<"#\n";
	ss<<ss_coldescr.str()<<endl;
	fprintf(fout,"%s",ss.str().c_str());

	//Loop sources
	bool deleteWCS= false;

	for(size_t k=0;k<sources.size();k++){
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No metadata are available to retrieve WCS!");
				#endif
				return -1;
			}
			wcs= metadata->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WCS from metadata!");
				#endif
				return -1;
			}
			deleteWCS= true;
		}

		//Print source info to ascii
		std::vector<std::string> slist= SourceComponentsToAscii(sources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightnessToFlux,delimiter);
		for(size_t j=0;j<slist.size();j++){
			fprintf(fout,"%s\n",slist[j].c_str());
		}
		
	}//end loop sources

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	//Close file
	fclose(fout);

	return 0;

}//close WriteComponentsToAscii()


const std::vector<std::string> SourceExporter::SourceComponentsToAscii(Source* source,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo,bool convertBrightnessToFlux,char delimiter)
{
	//Init vector
	std::vector<std::string> fitComponentStrList;
	
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return fitComponentStrList;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWCS(wcsType);
			if(wcs) deleteWCS= true;
			else {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get WCS from metadata!");
				#endif
			}
		}
		else {	
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
		}		
	}

	//Compute WCS J2000
	WCS* wcs_j2000= metadata->GetWCS(eJ2000);
	if(!wcs_j2000){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get J2000 WCS from metadata!");
		#endif
	}

	//Check if fit info are available
	bool hasFitInfo= source->HasFitInfo();
	if(hasFitInfo){
		//Retrieve fit pars
		SourceFitPars fitPars= source->GetFitPars();
		int nComponents= fitPars.GetNComponents();
	
		//Get spectral axis info
		double Nu= -999;
		double dNu= -999;
		std::string units= "";
		source->GetSpectralAxisInfo(Nu,dNu,units);
		CodeUtils::StripBlankSpaces(units);
		if(metadata && units=="Hz") {//convert to GHz
			Nu/= 1.e+9;
			dNu/= 1.e+9;
		}

		//Get total flux density
		double beamArea= source->GetBeamFluxIntegral();
		double fluxDensityTot= 0;
		double fluxDensityTot_err= 0;
		source->GetFluxDensity(fluxDensityTot);
		source->GetFluxDensityErr(fluxDensityTot_err);
		double fluxTot= fluxDensityTot/beamArea;
		double fluxTot_err= fluxDensityTot_err/beamArea;
		
		//Retrieve component spectral index data
		std::vector<SpectralIndexData> compSpectralIndexData= source->GetComponentSpectralIndexData();
		
		//Retrieve component astro objects
		std::vector<std::vector<AstroObject>> compAstroObjects= source->GetComponentAstroObjects();
	
		//Loop over fit components
		for(int k=0;k<nComponents;k++)
		{
			//Skip component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			//Get component fit pars
			//- Amplitude
			double A= fitPars.GetParValue(k,"A");
			double A_err= fitPars.GetParError(k,"A");

			//- Flux density
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			double fluxDensity_err= fitPars.GetComponentFluxDensityErr(k);
			double flux= fluxDensity/beamArea;
			double flux_err= fluxDensity_err/beamArea;

			//- Fit ellipse pars
			double x0= 0;
			double y0= 0;
			double bmaj= 0;
			double bmin= 0;
			double pa= 0;
			if(fitPars.GetComponentFitEllipsePars(k,x0,y0,bmaj,bmin,pa)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}
		
			double x0_err= 0;
			double y0_err= 0;
			double bmaj_err= 0;
			double bmin_err= 0;
			double pa_err= 0;
			if(fitPars.GetComponentFitEllipseParErrors(k,x0_err,y0_err,bmaj_err,bmin_err,pa_err)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Fit ellipse eccentricity, area & rot angle vs beam
			double E= fitPars.GetComponentFitEllipseEccentricity(k);
			double Area= fitPars.GetComponentFitEllipseArea(k);
			double RotAngle= fitPars.GetComponentFitEllipseRotAngleVSBeam(k);

			//- WCS fit ellipse pars
			double x0_wcs= 0;
			double y0_wcs= 0;
			double bmaj_wcs= 0;
			double bmin_wcs= 0;
			double pa_wcs= 0;
			if(fitPars.GetComponentFitWCSEllipsePars(k,x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve WCS ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}
			
			double x0_wcs_err= 0;
			double y0_wcs_err= 0;
			double bmaj_wcs_err= 0;
			double bmin_wcs_err= 0;
			double pa_wcs_err= 0;
			if(fitPars.GetComponentFitWCSEllipseParErrors(k,x0_wcs_err,y0_wcs_err,bmaj_wcs_err,bmin_wcs_err,pa_wcs_err)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve WCS ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Beam ellipse pars
			double bmaj_beam= 0;
			double bmin_beam= 0;
			double pa_beam= 0;
			if(fitPars.GetComponentBeamEllipsePars(k,bmaj_beam,bmin_beam,pa_beam)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve beam ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Beam eccentricity & area
			double E_beam= fitPars.GetComponentBeamEllipseEccentricity(k);
			double Area_beam= fitPars.GetComponentBeamEllipseArea(k);

			//- Ratio between fit ellipse pars & beam ellipse pars
			double EccentricityRatio= 0;
			if(E_beam!=0) EccentricityRatio= E/E_beam;
			
			double AreaRatio= 0;
			if(Area_beam!=0) AreaRatio= Area/Area_beam;

			//- WCS beam-deconvolved ellipse pars
			double bmaj_deconv_wcs= 0;
			double bmin_deconv_wcs= 0;
			double pa_deconv_wcs= 0;
			if(fitPars.GetComponentFitWCSDeconvolvedEllipsePars(k,bmaj_deconv_wcs,bmin_deconv_wcs,pa_deconv_wcs)<0){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Failed to retrieve WCS beam-deconvolved ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//Compute IAU name
			std::stringstream ssname;
			ssname<<source->GetName()<<"_fitcomp"<<k+1;
			std::string iau_default= ssname.str();		
			std::string iau= iau_default;
			if(wcs_j2000){
				std::string wcspos_str= "";
				AstroUtils::PixelToWCSStrCoords(wcspos_str,wcs_j2000,x0,y0);	
				int status= AstroUtils::GetIAUCoords(iau,wcspos_str);
				if(status<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute IAU name for component no. "<<k+1<<", assume default name "<<iau_default<<" ...");
					#endif
					iau= iau_default;		
				}
			}

			
			//Get component flag
			int componentFlag= -1;
			if(fitPars.GetComponentFlag(componentFlag,k)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve flag for component no. "<<k+1<<"!");
				#endif
			}

			//Get component type
			int componentType= -1;
			if(fitPars.GetComponentType(componentType,k)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve type for component no. "<<k+1<<"!");
				#endif
			}
	
			//## Fill source component
			std::stringstream ss;

			//- Mother source name & npix
			ss<<source->GetName()<<delimiter;
			ss<<source->NPix<<delimiter;

			//- Component id
			ss<<k+1<<delimiter;

			//- Component IAU name
			ss<<iau<<delimiter;
				
			//- Position in pixel coordinates
			ss<<x0<<delimiter<<y0<<delimiter;
			ss<<x0_err<<delimiter<<y0_err<<delimiter;

			//- Position in WCS coordinates
			ss<<setprecision(8)<<x0_wcs<<delimiter<<y0_wcs<<delimiter;
			ss<<setprecision(8)<<x0_wcs_err<<delimiter<<y0_wcs_err<<delimiter;

			//- Frequency
			ss<<Nu<<delimiter;
		
			//- Flux amplitude & density
			ss<<A<<delimiter<<A_err<<delimiter;
			if(convertBrightnessToFlux){
				ss<<flux<<delimiter<<flux_err<<delimiter;
				ss<<fluxTot<<delimiter<<fluxTot_err<<delimiter;
			}
			else{
				ss<<fluxDensity<<delimiter<<fluxDensity_err<<delimiter;
				ss<<fluxDensityTot<<delimiter<<fluxDensityTot_err<<delimiter;
			}
			ss<<beamArea<<delimiter;
	
			//- Ellipse pars in pixel coordinates
			ss<<bmaj<<delimiter<<bmin<<delimiter<<pa<<delimiter;
			ss<<bmaj_err<<delimiter<<bmin_err<<delimiter<<pa_err<<delimiter;
				
			//- Ellipse pars in WCS coordinates
			ss<<bmaj_wcs<<delimiter<<bmin_wcs<<delimiter<<pa_wcs<<delimiter;
			ss<<bmaj_wcs_err<<delimiter<<bmin_wcs_err<<delimiter<<pa_wcs_err<<delimiter;
		
			//Beam ellipse pars
			ss<<bmaj_beam<<delimiter<<bmin_beam<<delimiter<<pa_beam<<delimiter;

			//- Beam-deconvolved ellipse pars in WCS coordinates
			ss<<bmaj_deconv_wcs<<delimiter<<bmin_deconv_wcs<<delimiter<<pa_deconv_wcs<<delimiter;
			
			//- Fit ellipse vs beam ellipse pars
			ss<<EccentricityRatio<<delimiter<<AreaRatio<<delimiter<<RotAngle<<delimiter;

			//Bkg/noise estimators
			ss<<source->GetBkgSum()<<delimiter;
			ss<<source->GetBkgRMSSum()<<delimiter;

			//Fit chi2/ndf
			ss<<fitPars.GetChi2()<<delimiter<<fitPars.GetNDF()<<delimiter;

			//Source component flags
			ss<<fitPars.GetFitQuality()<<delimiter<<componentFlag<<delimiter<<componentType;

			//## Additional source info
			if(writeAdditionalSourceInfo)
			{	
				ss<<delimiter;

				//- Spectral index data
				if(compSpectralIndexData.size()==nComponents)
				{
					bool hasComponentSpectralIndexData= ((source->HasComponentSpectralIndexData()) && (compSpectralIndexData[k].hasSpectralIndex));
			
					ss<<hasComponentSpectralIndexData<<delimiter;
					ss<<compSpectralIndexData[k].isMultiSourceMatchIndex<<delimiter;
					ss<<compSpectralIndexData[k].spectralIndex<<delimiter;
					ss<<compSpectralIndexData[k].spectralIndexErr<<delimiter;
					ss<<compSpectralIndexData[k].isSpectralIndexFit<<delimiter;
					ss<<compSpectralIndexData[k].spectralFitChi2<<delimiter;
					ss<<compSpectralIndexData[k].spectralFitNDF<<delimiter;
				}
				else{
					SpectralIndexData dummySid;
					ss<<false<<delimiter;
					ss<<dummySid.isMultiSourceMatchIndex<<delimiter;
					ss<<dummySid.spectralIndex<<delimiter;
					ss<<dummySid.spectralIndexErr<<delimiter;
					ss<<dummySid.isSpectralIndexFit<<delimiter;
					ss<<dummySid.spectralFitChi2<<delimiter;
					ss<<dummySid.spectralFitNDF<<delimiter;
				}

				//- Astro object IDs
				/*
				if(source->componentObjLocationIds.size()==nComponents){
					ss<<source->componentObjLocationIds[k]<<delimiter;
				}
				else{
					ss<<eUNKNOWN_OBJECT_LOCATION<<delimiter;
				}
				*/
				if(source->componentObjClassStrIds.size()==nComponents){
					ss<<"\""<<source->componentObjClassStrIds[k]<<"\""<<delimiter;
				}
				else{
					ss<<"\"0000\""<<delimiter;
				}
				if(source->componentObjClassIds.size()==nComponents){
					ss<<source->componentObjClassIds[k]<<delimiter;
				}
				else{
					ss<<eUNKNOWN_OBJECT<<delimiter;
				}
				if(source->componentObjClassSubIds.size()==nComponents){
					ss<<source->componentObjClassSubIds[k]<<delimiter;
				}
				else{
					ss<<eUNKNOWN_OBJECT<<delimiter;
				}
				if(source->componentObjConfirmed.size()==nComponents){
					ss<<source->componentObjConfirmed[k]<<delimiter;
				}
				else{
					ss<<eUNKNOWN_OBJECT<<delimiter;
				}
				
							
				//- Astro object cross matches
				bool hasComponentAstroObjs= source->HasComponentAstroObjects();
		
				ss<<hasComponentAstroObjs<<delimiter;

				if(hasComponentAstroObjs && !compAstroObjects.empty() && !compAstroObjects[k].empty()){
					std::stringstream sstream;
					if(compAstroObjects[k].size()==1){
						sstream<<compAstroObjects[k][0].name;
					}
					else{
						for(size_t j=0;j<compAstroObjects[k].size()-1;j++){
							sstream<<compAstroObjects[k][j].name<<";";
						}
						sstream<<compAstroObjects[k][compAstroObjects[k].size()-1].name;
					}
			
					//ss<<sstream.str();
					ss<<"\""<<sstream.str()<<"\"";

				}//close if
				else{
					//ss<<"XXX";
					ss<<"\"XXX\"";
				}
		
			}//close if writeAdditionalSourceInfo
	

			//Store component string
			fitComponentStrList.push_back(ss.str());

		}//end loop components
	}//close if has fit

	
	//Store nested components
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++)
		{
			std::vector<std::string> fitComponentStrList_nested= SourceComponentsToAscii(nestedSources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightnessToFlux,delimiter);
			if(!fitComponentStrList_nested.empty()){
				fitComponentStrList.insert(fitComponentStrList.end(),fitComponentStrList_nested.begin(),fitComponentStrList_nested.end());
			}
		}
	}//close if has nested sources 
	
	
	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	return fitComponentStrList;

}//close SourceComponentsToAscii()





//=================================================
//==        JSON EXPORTER
//=================================================
int SourceExporter::WriteToJson(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WCS* wcs)
{
	//Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	
	//Loop sources
	bool deleteWCS= false;
	Json::Value root;
	std::vector<Json::Value> jsonValues;

	for(size_t k=0;k<sources.size();k++){
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No metadata are available to retrieve WCS!");
				#endif
				return -1;
			}
			wcs= metadata->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WCS from metadata!");
				#endif
				return -1;
			}
			deleteWCS= true;
		}

		//Get source in json format
		SourceToJson(jsonValues,sources[k],dumpNestedSourceInfo,wcsType,wcs);
		
	}//end loop sources

	for(size_t j=0;j<jsonValues.size();j++){
		root["sources"].append(jsonValues[j]);
	}

	//Get string from json
	std::string jsonString= "";
	bool minimizeJson= false;
	if(CodeUtils::JsonToString(jsonString,root,minimizeJson)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to convert json object to string!");
		#endif
		return -1;
	}


	//Write to file
	fprintf(fout,"%s",jsonString.c_str());

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	//Close file
	fclose(fout);

	return 0;

}//close WriteToJson()

int SourceExporter::SourceToJson(std::vector<Json::Value>& jsonValues,Source* source,bool dumpNestedSourceInfo,int wcsType,WCS* wcs)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return -1;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWCS(wcsType);
			if(wcs) deleteWCS= true;
			else {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get WCS from metadata!");
				#endif
			}
		}
		else {
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
		}		
	}
	
	//Compute IAU name
	std::string name= source->GetName();
	bool useWeightedPos= false;
	std::string iauName= name;
	if(wcs) iauName= source->GetIAUName(useWeightedPos,wcs,wcsType);

	//Compute WCS centroid
	double X0_wcs= 0;
	double Y0_wcs= 0;
	double X0_weighted_wcs= 0;
	double Y0_weighted_wcs= 0;
	if(wcs){
		source->GetWCSPos(X0_wcs,Y0_wcs,wcs,wcsType);
		source->GetWCSWeightedPos(X0_weighted_wcs,Y0_weighted_wcs,wcs,wcsType);
	}

	//Get spectral axis info
	double Nu= -999;
	double dNu= -999;
	std::string units= "";
	source->GetSpectralAxisInfo(Nu,dNu,units);
	CodeUtils::StripBlankSpaces(units);
	if(metadata && units=="Hz") {//convert to GHz
		Nu/= 1.e+9;
		dNu/= 1.e+9;
	}
	
	double fluxDensity= 0;
	double fluxDensityErr= 0;
	if(source->HasFitInfo()){
		source->GetFluxDensity(fluxDensity);
		source->GetFluxDensityErr(fluxDensityErr);
	}

	//Get WCS bounding box
	float xmin, xmax, ymin, ymax;
	source->GetSourceRange(xmin,xmax,ymin,ymax);
	
	double xmin_wcs, xmax_wcs, ymin_wcs, ymax_wcs;
	source->GetWCSSourceRange(xmin_wcs,xmax_wcs,ymin_wcs,ymax_wcs,wcs,wcsType);
		
	//## Set fields
	Json::Value json;

	//- Source name
	json["name"]= name;
	json["iauName"]= iauName;
	
	//- Number of pixels
	json["nPix"]= source->NPix;
	
	//- Number of sub-components
	json["nComponents"]= source->GetNSelFitComponents();
	json["nNestedSources"]= source->GetNestedSourceNumber();
	
	//- Pixel centroids
	json["X0"]= source->X0;
	json["Y0"]= source->Y0;
	json["X0w"]= source->GetSx();
	json["Y0w"]= source->GetSy();
	json["X0_wcs"]= X0_wcs;
	json["Y0_wcs"]= Y0_wcs;
	json["X0w_wcs"]= X0_weighted_wcs;
	json["Y0w_wcs"]= Y0_weighted_wcs;

	//- Bounding box
	json["Xmin"]= xmin;
	json["Xmax"]= xmax;
	json["Ymin"]= ymin;
	json["Ymax"]= ymax;

	json["Xmin_wcs"]= xmin_wcs;
	json["Xmax_wcs"]= xmax_wcs;
	json["Ymin_wcs"]= ymin_wcs;
	json["Ymax_wcs"]= ymax_wcs;

	//- Frequency
	json["nu"]= Nu;
	
	//- Flux 
	json["S"]= source->GetS();
	json["Smax"]= source->GetSmax();
	json["fluxDensity"]= fluxDensity;
	json["fluxDensityErr"]= fluxDensityErr;
	json["beamArea"]= source->GetBeamFluxIntegral();

	//Bkg/noise estimators
	json["bkgSum"]= source->GetBkgSum();
	json["bkgRMSSum"]= source->GetBkgRMSSum();
	
	//- Source flags
	json["type"]= source->Type;
	json["flag"]= source->Flag;
	json["isGoodSource"]= source->IsGoodSource();
	json["depthLevel"]= source->GetDepthLevel();
	
	jsonValues.push_back(json);
	
	//Store nested sources
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++)
		{
			SourceToJson(jsonValues,nestedSources[k],dumpNestedSourceInfo,wcsType,wcs);
		}
	}//close if has nested sources 
		
	//Append object
	//root["sources"].append(json);


	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	
	return 0;

}//close SourceToJson()

int SourceExporter::WriteComponentsToJson(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WCS* wcs)
{
	//Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Loop sources
	bool deleteWCS= false;
	Json::Value root;
	std::vector<Json::Value> jsonValues;

	for(size_t k=0;k<sources.size();k++){
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No metadata are available to retrieve WCS!");
				#endif
				return -1;
			}
			//wcs= metadata->GetWorldCoord(wcsType);
			wcs= metadata->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WCS from metadata!");
				#endif
				return -1;
			}
			deleteWCS= true;
		}

		//Print source info to ascii
		SourceComponentsToJson(jsonValues,sources[k],dumpNestedSourceInfo,wcsType,wcs);
		
		
	}//end loop sources

	for(size_t j=0;j<jsonValues.size();j++){
		root["components"].append(jsonValues[j]);
	}

	//Get string from json
	std::string jsonString= "";
	bool minimizeJson= false;
	if(CodeUtils::JsonToString(jsonString,root,minimizeJson)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to convert json object to string!");
		#endif
		return -1;
	}


	//Write to file
	fprintf(fout,"%s",jsonString.c_str());

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	//Close file
	fclose(fout);

	return 0;

}//close WriteComponentsToJson)

int SourceExporter::SourceComponentsToJson(std::vector<Json::Value>& jsonValues,Source* source,bool dumpNestedSourceInfo,int wcsType,WCS* wcs)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return -1;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWCS(wcsType);
			if(wcs) deleteWCS= true;
			else {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get WCS from metadata!");
				#endif
			}
		}
		else {	
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
		}		
	}


	//Check if fit info are available
	bool hasFitInfo= source->HasFitInfo();
	if(hasFitInfo){
		//Retrieve fit pars
		SourceFitPars fitPars= source->GetFitPars();
		int nComponents= fitPars.GetNComponents();
	
		//Get spectral axis info
		double Nu= -999;
		double dNu= -999;
		std::string units= "";
		source->GetSpectralAxisInfo(Nu,dNu,units);
		CodeUtils::StripBlankSpaces(units);
		if(metadata && units=="Hz") {//convert to GHz
			Nu/= 1.e+9;
			dNu/= 1.e+9;
		}

		//Get total flux density
		double fluxDensityTot= 0;
		double fluxDensityTot_err= 0;
		source->GetFluxDensity(fluxDensityTot);
		source->GetFluxDensityErr(fluxDensityTot_err);
	
	
		//Loop over fit components
		for(int k=0;k<nComponents;k++)
		{
			//Skip component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			//Get component fit pars
			//- Amplitude
			double A= fitPars.GetParValue(k,"A");
			double A_err= fitPars.GetParError(k,"A");

			//- Flux density
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			double fluxDensity_err= fitPars.GetComponentFluxDensityErr(k);

			//- Fit ellipse pars
			double x0= 0;
			double y0= 0;
			double bmaj= 0;
			double bmin= 0;
			double pa= 0;
			if(fitPars.GetComponentFitEllipsePars(k,x0,y0,bmaj,bmin,pa)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}
		
			double x0_err= 0;
			double y0_err= 0;
			double bmaj_err= 0;
			double bmin_err= 0;
			double pa_err= 0;
			if(fitPars.GetComponentFitEllipseParErrors(k,x0_err,y0_err,bmaj_err,bmin_err,pa_err)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Fit ellipse eccentricity, area & rot angle vs beam
			double E= fitPars.GetComponentFitEllipseEccentricity(k);
			double Area= fitPars.GetComponentFitEllipseArea(k);
			double RotAngle= fitPars.GetComponentFitEllipseRotAngleVSBeam(k);

			//- WCS fit ellipse pars
			double x0_wcs= 0;
			double y0_wcs= 0;
			double bmaj_wcs= 0;
			double bmin_wcs= 0;
			double pa_wcs= 0;
			if(fitPars.GetComponentFitWCSEllipsePars(k,x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve WCS ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}
			
			double x0_wcs_err= 0;
			double y0_wcs_err= 0;
			double bmaj_wcs_err= 0;
			double bmin_wcs_err= 0;
			double pa_wcs_err= 0;
			if(fitPars.GetComponentFitWCSEllipseParErrors(k,x0_wcs_err,y0_wcs_err,bmaj_wcs_err,bmin_wcs_err,pa_wcs_err)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve WCS ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Beam ellipse pars
			double bmaj_beam= 0;
			double bmin_beam= 0;
			double pa_beam= 0;
			if(fitPars.GetComponentBeamEllipsePars(k,bmaj_beam,bmin_beam,pa_beam)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve beam ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Beam eccentricity & area
			double E_beam= fitPars.GetComponentBeamEllipseEccentricity(k);
			double Area_beam= fitPars.GetComponentBeamEllipseArea(k);

			//- Ratio between fit ellipse pars & beam ellipse pars
			double EccentricityRatio= 0;
			if(E_beam!=0) EccentricityRatio= E/E_beam;
			
			double AreaRatio= 0;
			if(Area_beam!=0) AreaRatio= Area/Area_beam;

			//- WCS beam-deconvolved ellipse pars
			double bmaj_deconv_wcs= 0;
			double bmin_deconv_wcs= 0;
			double pa_deconv_wcs= 0;
			if(fitPars.GetComponentFitWCSDeconvolvedEllipsePars(k,bmaj_deconv_wcs,bmin_deconv_wcs,pa_deconv_wcs)<0){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Failed to retrieve WCS beam-deconvolved ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//Compute IAU name
			std::stringstream ssname;
			ssname<<source->GetName()<<"_fitcomp"<<k+1;
			std::string iau_default= ssname.str();		
			std::string iau= iau_default;
			if(wcs){
				std::string wcspos_str= "";
				AstroUtils::PixelToWCSStrCoords(wcspos_str,wcs,x0,y0);	
				int status= AstroUtils::GetIAUCoords(iau,wcspos_str);
				if(status<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute IAU name for component no. "<<k+1<<", assume default name "<<iau_default<<" ...");
					#endif
					iau= iau_default;		
				}
			}

			
			//Get component flag
			int componentFlag= -1;
			if(fitPars.GetComponentFlag(componentFlag,k)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve flag for component no. "<<k+1<<"!");
				#endif
			}

			//Get component type
			int componentType= -1;
			if(fitPars.GetComponentType(componentType,k)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve type for component no. "<<k+1<<"!");
				#endif
			}
	
			//## Fill source component
			Json::Value json;

			//- Mother source name & npix
			json["name"]= source->GetName();
			json["nPix"]= source->NPix;
			
			//- Component id
			json["componentId"]= k+1;
			
			//- Component IAU name
			json["iauName"]= iau;
				
			//- Position in pixel coordinates
			json["x0"]= x0;
			json["y0"]= y0;
			json["x0Err"]= x0_err;
			json["y0Err"]= y0_err;

			//- Position in WCS coordinates
			json["x0_wcs"]= x0_wcs;
			json["y0_wcs"]= y0_wcs;
			json["x0Err_wcs"]= x0_wcs_err;
			json["y0Err_wcs"]= y0_wcs_err;


			//- Frequency
			json["nu"]= Nu;
		
			//- Flux amplitude & density
			json["A"]= A;
			json["Aerr"]= A_err;
			json["fluxDensity"]= fluxDensity;
			json["fluxDensityErr"]= fluxDensity_err;
			json["fluxDensityTot"]= fluxDensityTot;
			json["fluxDensityTotErr"]= fluxDensityTot_err;
			json["beamArea"]= source->GetBeamFluxIntegral();
	
			//- Ellipse pars in pixel coordinates
			json["bmaj"]= bmaj;
			json["bmin"]= bmin;
			json["pa"]= pa;
			json["bmajErr"]= bmaj_err;
			json["bminErr"]= bmin_err;
			json["paErr"]= pa_err;
					
			//- Ellipse pars in WCS coordinates
			json["bmaj_wcs"]= bmaj_wcs;
			json["bmin_wcs"]= bmin_wcs;
			json["pa_wcs"]= pa_wcs;
 			json["bmajErr_wcs"]= bmaj_wcs_err;
			json["bminErr_wcs"]= bmin_wcs_err;
			json["paErr_wcs"]= pa_wcs_err;
 		
			//Beam ellipse pars
			json["bmaj_beam"]= bmaj_beam;
			json["bmin_beam"]= bmin_beam;
			json["pa_beam"]= pa_beam;

			//- Beam-deconvolved ellipse pars in WCS coordinates
			json["bmaj_deconv_wcs"]= bmaj_deconv_wcs;
			json["bmin_deconv_wcs"]= bmin_deconv_wcs;
			json["pa_deconv_wcs"]= pa_deconv_wcs;
 
			//- Fit ellipse vs beam ellipse pars
			json["eccentricityRatio"]= EccentricityRatio;
			json["areaRatio"]= AreaRatio;
			json["rotAngle"]= RotAngle;
	
			//Bkg/noise estimators	
			json["bkgSum"]= source->GetBkgSum();
			json["bkgRMSSum"]= source->GetBkgRMSSum();

			//Fit chi2/ndf
			json["chi2"]= fitPars.GetChi2();
			json["ndf"]= fitPars.GetNDF();

			//Source component flags
			json["fitQuality"]= fitPars.GetFitQuality();
			json["flag"]= componentFlag;
			json["type"]= componentType;

			//Add json value to collection
			jsonValues.push_back(json);

		}//end loop components
	}//close if has fit

	
	//Store nested components
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++)
		{
			SourceComponentsToJson(jsonValues,nestedSources[k],dumpNestedSourceInfo,wcsType,wcs);
		}
	}//close if has nested sources 
	
	
	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	return 0;

}//close SourceComponentsToJson()



//=================================================
//==        ROOT EXPORTER
//=================================================
int SourceExporter::WriteToROOT(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo)
{
	//Open output file
	TFile* fout= new TFile(filename.c_str(),"RECREATE");
	if(!fout || !fout->IsOpen()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to create out ROOT file "<<filename<<"!");
		#endif
		return -1;
	}

	//Create TTree
	SourceTreeData sourceTreeData;
	
	TTree* dataTree= new TTree("data","data");
	dataTree->Branch("name",&sourceTreeData.name);//source name
	dataTree->Branch("iau",&sourceTreeData.iau);//source iau name
	dataTree->Branch("nPix",&sourceTreeData.nPix);//number of pixels
	dataTree->Branch("nFitComponents",&sourceTreeData.nFitComponents);//number of fit components
	dataTree->Branch("nSelFitComponents",&sourceTreeData.nSelFitComponents);//number of fit components
	dataTree->Branch("nNestedSources",&sourceTreeData.nNestedSources);//number of nested sources
	dataTree->Branch("X0",&sourceTreeData.X0);//pos X0 in image coords
	dataTree->Branch("Y0",&sourceTreeData.Y0);//pos Y0 in image coords
	dataTree->Branch("X0w",&sourceTreeData.X0w);//signal-weighted pos X0 in image coords
	dataTree->Branch("Y0w",&sourceTreeData.Y0w);//signal-weighted pos Y0 in image coords
	dataTree->Branch("X0_wcs",&sourceTreeData.X0_wcs);//pos X0 in WCS coords
	dataTree->Branch("Y0_wcs",&sourceTreeData.Y0_wcs);//pos Y0 in WCS coords
	dataTree->Branch("X0w_wcs",&sourceTreeData.X0w_wcs);//signal-weighted pos X0 in WCS coords
	dataTree->Branch("Y0w_wcs",&sourceTreeData.Y0w_wcs);//signal-weighted pos Y0 in WCS coords
	dataTree->Branch("Xmin",&sourceTreeData.Xmin);//bounding box Xmin
	dataTree->Branch("Xmax",&sourceTreeData.Xmax);//bounding box Xmax
	dataTree->Branch("Ymin",&sourceTreeData.Ymin);//bounding box Ymin
	dataTree->Branch("Ymax",&sourceTreeData.Ymax);//bounding box Ymax
	dataTree->Branch("Xmin_wcs",&sourceTreeData.Xmin_wcs);//bounding box Xmin in WCS coords
	dataTree->Branch("Xmax_wcs",&sourceTreeData.Xmax_wcs);//bounding box Xmax in WCS coords
	dataTree->Branch("Ymin_wcs",&sourceTreeData.Ymin_wcs);//bounding box Ymin in WCS coords
	dataTree->Branch("Ymax_wcs",&sourceTreeData.Ymax_wcs);//bounding box Ymax in WCS coords
	dataTree->Branch("Nu",&sourceTreeData.Nu);//frequency
	dataTree->Branch("S",&sourceTreeData.S);//flux sum over pixels
	dataTree->Branch("Smax",&sourceTreeData.Smax);//flux max
	dataTree->Branch("fittedFlux",&sourceTreeData.fittedFlux);//flux density
	dataTree->Branch("fittedFluxErr",&sourceTreeData.fittedFluxErr);//flux density err
	dataTree->Branch("beamArea",&sourceTreeData.beamArea);//beamArea
	dataTree->Branch("bkgSum",&sourceTreeData.bkgSum);//bkg sum over pixels
	dataTree->Branch("rmsSum",&sourceTreeData.rmsSum);//noise rms sum over pixels
	dataTree->Branch("type",&sourceTreeData.type);//source type
	dataTree->Branch("flag",&sourceTreeData.flag);//source flag
	dataTree->Branch("good",&sourceTreeData.good);//source isGoodFlag
	dataTree->Branch("depthLevel",&sourceTreeData.depthLevel);//source depth level
	dataTree->Branch("chi2",&sourceTreeData.chi2);
	dataTree->Branch("ndf",&sourceTreeData.ndf);
	dataTree->Branch("fitQuality",&sourceTreeData.fitQuality);
	dataTree->Branch("residualMean",&sourceTreeData.residualMean);
	dataTree->Branch("residualRMS",&sourceTreeData.residualRMS);
	dataTree->Branch("residualMedian",&sourceTreeData.residualMedian);
	dataTree->Branch("residualMAD",&sourceTreeData.residualMAD);
	dataTree->Branch("residualMin",&sourceTreeData.residualMin);
	dataTree->Branch("residualMax",&sourceTreeData.residualMax);

	if(writeAdditionalSourceInfo)
	{
		//Spectral index
		dataTree->Branch("hasSpectralIndexData",&sourceTreeData.hasSpectralIndexData);
		dataTree->Branch("isMultiSourceMatchIndex",&sourceTreeData.isMultiSourceMatchIndex);
		dataTree->Branch("spectralIndex",&sourceTreeData.spectralIndex);
		dataTree->Branch("spectralIndexErr",&sourceTreeData.spectralIndexErr);
		dataTree->Branch("isSpectralIndexFit",&sourceTreeData.isSpectralIndexFit);
		dataTree->Branch("spectralFitChi2",&sourceTreeData.spectralFitChi2);
		dataTree->Branch("spectralFitNDF",&sourceTreeData.spectralFitNDF);

		//Astro obj match
		dataTree->Branch("objLocationId",&sourceTreeData.objLocationId);
		dataTree->Branch("objClassId",&sourceTreeData.objClassId);
		dataTree->Branch("objClassSubId",&sourceTreeData.objClassSubId);
		dataTree->Branch("objConfirmed",&sourceTreeData.objConfirmed);
		dataTree->Branch("hasAstroObjs",&sourceTreeData.hasAstroObjs);
		dataTree->Branch("objName",&sourceTreeData.objName);
	}

	//Loop sources
	bool deleteWCS= false;
	
	for(size_t k=0;k<sources.size();k++)
	{
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No metadata are available to retrieve WCS!");
				#endif
				return -1;
			}
			wcs= metadata->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WCS from metadata!");
				#endif
				return -1;
			}
			deleteWCS= true;
		}

		//Get source data and fill TTree
		FillSourceTTree(dataTree,sourceTreeData,sources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo);
		
	}//end loop sources

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	//Write TTree to file
	fout->cd();
	dataTree->Write();
	fout->Close();

	return 0;

}//close WriteToROOT()

int SourceExporter::FillSourceTTree(TTree* dataTree,SourceTreeData& sourceTreeData,Source* source,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return -1;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWCS(wcsType);
			if(wcs) deleteWCS= true;
			else {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get WCS from metadata!");
				#endif
			}
		}
		else {
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
		}		
	}

	//## Fill source data struct
	//Get source name	
	sourceTreeData.name= source->GetName();

	//Compute IAU name
	bool useWeightedPos= false;
	sourceTreeData.iau= source->GetName();
	if(wcs) sourceTreeData.iau= source->GetIAUName(useWeightedPos,wcs,wcsType);

	//Get number of pixels
	sourceTreeData.nPix= source->NPix;

	//GEt number of nested sources & fit components
	sourceTreeData.nFitComponents= source->GetNFitComponents();
	sourceTreeData.nSelFitComponents= source->GetNSelFitComponents();
	sourceTreeData.nNestedSources= source->GetNestedSourceNumber();

	//Get source centroid
	sourceTreeData.X0= source->X0;
	sourceTreeData.Y0= source->Y0;
	sourceTreeData.X0w= source->GetSx();
	sourceTreeData.Y0w= source->GetSy();

	//Compute WCS centroid
	sourceTreeData.X0_wcs= 0;
	sourceTreeData.Y0_wcs= 0;
	sourceTreeData.X0w_wcs= 0;
	sourceTreeData.Y0w_wcs= 0;
	if(wcs){
		source->GetWCSPos(sourceTreeData.X0_wcs,sourceTreeData.Y0_wcs,wcs,wcsType);
		source->GetWCSWeightedPos(sourceTreeData.X0w_wcs,sourceTreeData.Y0w_wcs,wcs,wcsType);
	}

	//Get spectral axis info
	sourceTreeData.Nu= -999;
	double dNu= -999;
	std::string units= "";
	source->GetSpectralAxisInfo(sourceTreeData.Nu,dNu,units);
	CodeUtils::StripBlankSpaces(units);
	if(metadata && units=="Hz") {//convert to GHz
		sourceTreeData.Nu/= 1.e+9;
		dNu/= 1.e+9;
	}
	
	sourceTreeData.fittedFlux= 0;
	sourceTreeData.fittedFluxErr= 0;
	if(source->HasFitInfo()){
		source->GetFluxDensity(sourceTreeData.fittedFlux);
		source->GetFluxDensityErr(sourceTreeData.fittedFluxErr);
	}

	//Get WCS bounding box
	float xmin, xmax, ymin, ymax;
	source->GetSourceRange(xmin,xmax,ymin,ymax);
	sourceTreeData.Xmin= xmin;
	sourceTreeData.Xmax= xmax;
	sourceTreeData.Ymin= ymin;
	sourceTreeData.Ymax= ymax;
	
	double xmin_wcs, xmax_wcs, ymin_wcs, ymax_wcs;
	source->GetWCSSourceRange(xmin_wcs,xmax_wcs,ymin_wcs,ymax_wcs,wcs,wcsType);
	sourceTreeData.Xmin_wcs= xmin_wcs;	
	sourceTreeData.Xmax_wcs= xmax_wcs;	
	sourceTreeData.Ymin_wcs= ymin_wcs;	
	sourceTreeData.Ymax_wcs= ymax_wcs;	
	
	//Get flux 
	sourceTreeData.S= source->GetS();
	sourceTreeData.Smax= source->GetSmax();
	sourceTreeData.beamArea= source->GetBeamFluxIntegral();

	//Bkg/noise estimators
	sourceTreeData.bkgSum= source->GetBkgSum();
	sourceTreeData.rmsSum= source->GetBkgRMSSum();

	//Flags	
	sourceTreeData.type= source->Type;
	sourceTreeData.flag= source->Flag;
	sourceTreeData.good= static_cast<int>(source->IsGoodSource());
	sourceTreeData.depthLevel= source->GetDepthLevel();

	//Fit chi2 & residuals
	int fitQuality= eUnknownFitQuality;
	double chi2= -999;
	double ndf= -999;
	double residualMean= -999;
	double residualRMS= -999;
	double residualMedian= -999;
	double residualMAD= -999;
	double residualMin= -999;
	double residualMax= -999;

	if(source->HasFitInfo()){
		SourceFitPars fitPars= source->GetFitPars();

		chi2= fitPars.GetChi2();
		ndf= fitPars.GetNDF();
		fitQuality= fitPars.GetFitQuality();
		residualMean= fitPars.GetResidualMean();
		residualRMS= fitPars.GetResidualRMS();
		residualMedian= fitPars.GetResidualMedian();
		residualMAD= fitPars.GetResidualMAD();
		residualMin= fitPars.GetResidualMin();
		residualMax= fitPars.GetResidualMax();
	}
		
	sourceTreeData.chi2= chi2;
	sourceTreeData.ndf= ndf;
	sourceTreeData.fitQuality= fitQuality;
	sourceTreeData.residualMean= residualMean;
	sourceTreeData.residualRMS= residualRMS;
	sourceTreeData.residualMedian= residualMedian;
	sourceTreeData.residualMAD= residualMAD;
	sourceTreeData.residualMin= residualMin;
	sourceTreeData.residualMax= residualMax;
	
	//- Spectral index data
	SpectralIndexData sid= source->GetSpectralIndexData();	
	bool hasSpectralIndexData= ((source->HasSpectralIndexData()) && (sid.hasSpectralIndex));
			
	
	sourceTreeData.hasSpectralIndexData= hasSpectralIndexData;
	sourceTreeData.isMultiSourceMatchIndex = sid.isMultiSourceMatchIndex;
	sourceTreeData.spectralIndex = sid.spectralIndex;
	sourceTreeData.spectralIndexErr = sid.spectralIndexErr;
	sourceTreeData.isSpectralIndexFit = sid.isSpectralIndexFit;
	sourceTreeData.spectralFitChi2 = sid.spectralFitChi2;
	sourceTreeData.spectralFitNDF = sid.spectralFitNDF;
	
	//- Astro object IDs
	sourceTreeData.objLocationId = source->ObjLocationId;
	sourceTreeData.objClassId = source->ObjClassId;
	sourceTreeData.objClassSubId = source->ObjClassSubId;
	sourceTreeData.objConfirmed = source->ObjConfirmed;

	//- Astro object cross matches
	bool hasAstroObjs= source->HasAstroObjects();
	std::vector<AstroObject> astroObjs= source->GetAstroObjects();
		
	sourceTreeData.hasAstroObjs= hasAstroObjs;

	if(hasAstroObjs && !astroObjs.empty()){
		std::stringstream sstream;
		if(astroObjs.size()==1){
			sstream<<astroObjs[0].name;
		}
		else{
			for(size_t j=0;j<astroObjs.size()-1;j++){
				sstream<<astroObjs[j].name<<";";
			}
			sstream<<astroObjs[astroObjs.size()-1].name;
		}
			
		sourceTreeData.objName= sstream.str();

	}//close if
	else{
		sourceTreeData.objName= "XXX";
	}

	//## Fill Tree
	dataTree->Fill();

	//## Store nested sources
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++)
		{			
			FillSourceTTree(dataTree,sourceTreeData,nestedSources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo);		
		}
	}//close if has nested sources 
		

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	return 0;

}//close FillSourceTTree()





int SourceExporter::WriteComponentsToROOT(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo)
{
	//Open output file
	TFile* fout= new TFile(filename.c_str(),"RECREATE");
	if(!fout || !fout->IsOpen()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to create out ROOT file "<<filename<<"!");
		#endif
		return -1;
	}

	//Create TTree
	SourceComponentTreeData sourceTreeData;
	
	TTree* dataTree= new TTree("data","data");
	dataTree->Branch("name",&sourceTreeData.name);
	dataTree->Branch("nPix",&sourceTreeData.nPix);
	dataTree->Branch("componentId",&sourceTreeData.componentId);
	dataTree->Branch("iau",&sourceTreeData.iau);	
	dataTree->Branch("nFitComponents",&sourceTreeData.nFitComponents);
	dataTree->Branch("nSelFitComponents",&sourceTreeData.nSelFitComponents);
	dataTree->Branch("X0",&sourceTreeData.X0);
	dataTree->Branch("Y0",&sourceTreeData.Y0);
	dataTree->Branch("X0_err",&sourceTreeData.X0_err);
	dataTree->Branch("Y0_err",&sourceTreeData.Y0_err);
	dataTree->Branch("X0_wcs",&sourceTreeData.X0_wcs);
	dataTree->Branch("Y0_wcs",&sourceTreeData.Y0_wcs);
	dataTree->Branch("X0_err_wcs",&sourceTreeData.X0_err_wcs);
	dataTree->Branch("Y0_err_wcs",&sourceTreeData.Y0_err_wcs);	
	dataTree->Branch("Nu",&sourceTreeData.Nu);
	dataTree->Branch("S",&sourceTreeData.S);
	dataTree->Branch("Smax",&sourceTreeData.Smax);
	dataTree->Branch("A",&sourceTreeData.A);
	dataTree->Branch("A_err",&sourceTreeData.A_err);
	dataTree->Branch("fittedFlux",&sourceTreeData.fittedFlux);
	dataTree->Branch("fittedFlux_err",&sourceTreeData.fittedFlux_err);
	dataTree->Branch("fittedIslandFlux",&sourceTreeData.fittedIslandFlux);
	dataTree->Branch("fittedIslandFlux_err",&sourceTreeData.fittedIslandFlux_err);
	dataTree->Branch("beamArea",&sourceTreeData.beamArea);
	dataTree->Branch("Bmaj",&sourceTreeData.Bmaj);
	dataTree->Branch("Bmin",&sourceTreeData.Bmin);
	dataTree->Branch("Pa",&sourceTreeData.Pa);
	dataTree->Branch("Bmaj_err",&sourceTreeData.Bmaj_err);
	dataTree->Branch("Bmin_err",&sourceTreeData.Bmin_err);
	dataTree->Branch("Pa_err",&sourceTreeData.Pa_err);
	dataTree->Branch("Bmaj_wcs",&sourceTreeData.Bmaj_wcs);
	dataTree->Branch("Bmin_wcs",&sourceTreeData.Bmin_wcs);
	dataTree->Branch("Pa_wcs",&sourceTreeData.Pa_wcs);
	dataTree->Branch("Bmaj_err_wcs",&sourceTreeData.Bmaj_err_wcs);
	dataTree->Branch("Bmin_err_wcs",&sourceTreeData.Bmin_err_wcs);
	dataTree->Branch("Pa_err_wcs",&sourceTreeData.Pa_err_wcs);
	dataTree->Branch("Bmaj_beam",&sourceTreeData.Bmaj_beam);
	dataTree->Branch("Bmin_beam",&sourceTreeData.Bmin_beam);
	dataTree->Branch("Pa_beam",&sourceTreeData.Pa_beam);
	dataTree->Branch("Bmaj_deconv_wcs",&sourceTreeData.Bmaj_deconv_wcs);
	dataTree->Branch("Bmin_deconv_wcs",&sourceTreeData.Bmin_deconv_wcs);
	dataTree->Branch("Pa_deconv_wcs",&sourceTreeData.Pa_deconv_wcs);
	dataTree->Branch("eccentricityRatio",&sourceTreeData.eccentricityRatio);
	dataTree->Branch("areaRatio",&sourceTreeData.areaRatio);
	dataTree->Branch("fitVSBeamRotAngle",&sourceTreeData.fitVSBeamRotAngle);
	dataTree->Branch("bkgSum",&sourceTreeData.bkgSum);
	dataTree->Branch("rmsSum",&sourceTreeData.rmsSum);
	dataTree->Branch("chi2",&sourceTreeData.chi2);
	dataTree->Branch("ndf",&sourceTreeData.ndf);
	dataTree->Branch("fitQuality",&sourceTreeData.fitQuality);
	dataTree->Branch("type",&sourceTreeData.type);
	dataTree->Branch("flag",&sourceTreeData.flag);

	if(writeAdditionalSourceInfo)
	{
		//Spectral index
		dataTree->Branch("hasSpectralIndexData",&sourceTreeData.hasSpectralIndexData);
		dataTree->Branch("isMultiSourceMatchIndex",&sourceTreeData.isMultiSourceMatchIndex);
		dataTree->Branch("spectralIndex",&sourceTreeData.spectralIndex);
		dataTree->Branch("spectralIndexErr",&sourceTreeData.spectralIndexErr);
		dataTree->Branch("isSpectralIndexFit",&sourceTreeData.isSpectralIndexFit);
		dataTree->Branch("spectralFitChi2",&sourceTreeData.spectralFitChi2);
		dataTree->Branch("spectralFitNDF",&sourceTreeData.spectralFitNDF);

		//Astro obj match
		dataTree->Branch("objLocationId",&sourceTreeData.objLocationId);
		dataTree->Branch("objClassId",&sourceTreeData.objClassId);
		dataTree->Branch("objClassSubId",&sourceTreeData.objClassSubId);
		dataTree->Branch("objConfirmed",&sourceTreeData.objConfirmed);
		dataTree->Branch("hasAstroObjs",&sourceTreeData.hasAstroObjs);
		dataTree->Branch("objName",&sourceTreeData.objName);
	}
	
	//Loop sources
	bool deleteWCS= false;
	
	for(size_t k=0;k<sources.size();k++)
	{
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No metadata are available to retrieve WCS!");
				#endif
				return -1;
			}
			wcs= metadata->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WCS from metadata!");
				#endif
				return -1;
			}
			deleteWCS= true;
		}

		//Get source data and fill TTree
		FillSourceComponentTree(dataTree,sourceTreeData,sources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo);
		
	}//end loop sources

	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	//Write TTree to file
	fout->cd();
	dataTree->Write();
	fout->Close();

	return 0;

}//close WriteComponentsToROOT()


int SourceExporter::FillSourceComponentTree(TTree* dataTree,SourceComponentTreeData& sourceData,Source* source,bool dumpNestedSourceInfo,int wcsType,WCS* wcs,bool writeAdditionalSourceInfo)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return -1;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWCS(wcsType);
			if(wcs) deleteWCS= true;
			else {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get WCS from metadata!");
				#endif
			}
		}
		else {	
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
		}		
	}

	//Get island info
	double S= source->GetS();
	double Smax= source->GetSmax();
	sourceData.S= S;
	sourceData.Smax= Smax;

	//Check if fit info are available
	bool hasFitInfo= source->HasFitInfo();

	if(hasFitInfo){
		//Retrieve fit pars
		SourceFitPars fitPars= source->GetFitPars();
		int nComponents= fitPars.GetNComponents();
		int nSelComponents= fitPars.GetNSelComponents();
		sourceData.nFitComponents= nComponents;
		sourceData.nSelFitComponents= nSelComponents;

		//Get spectral axis info
		sourceData.Nu= -999;
		double dNu= -999;
		std::string units= "";
		source->GetSpectralAxisInfo(sourceData.Nu,dNu,units);
		CodeUtils::StripBlankSpaces(units);
		if(metadata && units=="Hz") {//convert to GHz
			sourceData.Nu/= 1.e+9;
			dNu/= 1.e+9;
		}

		//Get total flux density
		sourceData.fittedIslandFlux= 0;
		sourceData.fittedIslandFlux_err= 0;
		source->GetFluxDensity(sourceData.fittedIslandFlux);
		source->GetFluxDensityErr(sourceData.fittedIslandFlux_err);
	
		//Get source name
		sourceData.name= source->GetName();

		//Get source nPix
		sourceData.nPix= source->NPix;

		//- Bkg pars
		sourceData.bkgSum= source->GetBkgSum();
		sourceData.rmsSum= source->GetBkgRMSSum();

		//Retrieve component spectral index data
		std::vector<SpectralIndexData> compSpectralIndexData= source->GetComponentSpectralIndexData();
		
		//Retrieve component astro objects
		std::vector<std::vector<AstroObject>> compAstroObjects= source->GetComponentAstroObjects();
	
		//Loop over fit components
		for(int k=0;k<nComponents;k++)
		{
			//Skip component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			//Get component fit pars
			//- Component id
			sourceData.componentId= k+1;

			//- Amplitude
			sourceData.A= fitPars.GetParValue(k,"A");
			sourceData.A_err= fitPars.GetParError(k,"A");

			//- Flux density
			sourceData.fittedFlux= fitPars.GetComponentFluxDensity(k);
			sourceData.fittedFlux_err= fitPars.GetComponentFluxDensityErr(k);
			sourceData.beamArea= source->GetBeamFluxIntegral();

			//- Fit ellipse pars
			sourceData.X0= 0;
			sourceData.Y0= 0;
			sourceData.Bmaj= 0;
			sourceData.Bmin= 0;
			sourceData.Pa= 0;
			if(fitPars.GetComponentFitEllipsePars(k,sourceData.X0,sourceData.Y0,sourceData.Bmaj,sourceData.Bmin,sourceData.Pa)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}
		
			sourceData.X0_err= 0;
			sourceData.Y0_err= 0;
			sourceData.Bmaj_err= 0;
			sourceData.Bmin_err= 0;
			sourceData.Pa_err= 0;
			if(fitPars.GetComponentFitEllipseParErrors(k,sourceData.X0_err,sourceData.Y0_err,sourceData.Bmaj_err,sourceData.Bmin_err,sourceData.Pa_err)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Fit ellipse eccentricity, area & rot angle vs beam
			double E= fitPars.GetComponentFitEllipseEccentricity(k);
			double Area= fitPars.GetComponentFitEllipseArea(k);
			

			//- WCS fit ellipse pars
			sourceData.X0_wcs= 0;
			sourceData.Y0_wcs= 0;
			sourceData.Bmaj_wcs= 0;
			sourceData.Bmin_wcs= 0;
			sourceData.Pa_wcs= 0;
			if(fitPars.GetComponentFitWCSEllipsePars(k,sourceData.X0_wcs,sourceData.Y0_wcs,sourceData.Bmaj_wcs,sourceData.Bmin_wcs,sourceData.Pa_wcs)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve WCS ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}
			
			sourceData.X0_err_wcs= 0;
			sourceData.Y0_err_wcs= 0;
			sourceData.Bmaj_err_wcs= 0;
			sourceData.Bmin_err_wcs= 0;
			sourceData.Pa_err_wcs= 0;
			if(fitPars.GetComponentFitWCSEllipseParErrors(k,sourceData.X0_err_wcs,sourceData.Y0_err_wcs,sourceData.Bmaj_err_wcs,sourceData.Bmin_err_wcs,sourceData.Pa_err_wcs)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve WCS ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Beam ellipse pars
			sourceData.Bmaj_beam= 0;
			sourceData.Bmin_beam= 0;
			sourceData.Pa_beam= 0;
			if(fitPars.GetComponentBeamEllipsePars(k,sourceData.Bmaj_beam,sourceData.Bmin_beam,sourceData.Pa_beam)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve beam ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Beam eccentricity & area
			double E_beam= fitPars.GetComponentBeamEllipseEccentricity(k);
			double Area_beam= fitPars.GetComponentBeamEllipseArea(k);

			//- Ratio between fit ellipse pars & beam ellipse pars
			sourceData.eccentricityRatio= 0;
			if(E_beam!=0) sourceData.eccentricityRatio= E/E_beam;
			
			sourceData.areaRatio= 0;
			if(Area_beam!=0) sourceData.areaRatio= Area/Area_beam;

			double dtheta= fitPars.GetComponentFitEllipseRotAngleVSBeam(k);
			sourceData.fitVSBeamRotAngle= MathUtils::GetAngleInRange(dtheta,90.);

			//- WCS beam-deconvolved ellipse pars
			sourceData.Bmaj_deconv_wcs= 0;
			sourceData.Bmin_deconv_wcs= 0;
			sourceData.Pa_deconv_wcs= 0;
			if(fitPars.GetComponentFitWCSDeconvolvedEllipsePars(k,sourceData.Bmaj_deconv_wcs,sourceData.Bmin_deconv_wcs,sourceData.Pa_deconv_wcs)<0){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Failed to retrieve WCS beam-deconvolved ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
				#endif
			}

			//- Compute IAU name
			std::stringstream ssname;
			ssname<<source->GetName()<<"_fitcomp"<<k+1;
			std::string iau_default= ssname.str();		
			sourceData.iau= iau_default;
			if(wcs){
				std::string wcspos_str= "";
				AstroUtils::PixelToWCSStrCoords(wcspos_str,wcs,sourceData.X0,sourceData.Y0);	
				int status= AstroUtils::GetIAUCoords(sourceData.iau,wcspos_str);
				if(status<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute IAU name for component no. "<<k+1<<", assume default name "<<iau_default<<" ...");
					#endif
					sourceData.iau= iau_default;		
				}
			}

			//Get component flag
			sourceData.flag= -1;
			if(fitPars.GetComponentFlag(sourceData.flag,k)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve flag for component no. "<<k+1<<"!");
				#endif
			}

			//Get component type
			sourceData.type= -1;
			if(fitPars.GetComponentType(sourceData.type,k)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to retrieve type for component no. "<<k+1<<"!");
				#endif
			}
			
			//- Fit chi2/ndf, quality
			sourceData.chi2= fitPars.GetChi2();
			sourceData.ndf= fitPars.GetNDF();
			sourceData.fitQuality= fitPars.GetFitQuality();			
	
			//- Spectral index data
			if(compSpectralIndexData.size()==nComponents)
			{
				bool hasComponentSpectralIndexData= ((source->HasComponentSpectralIndexData()) && (compSpectralIndexData[k].hasSpectralIndex));
				sourceData.hasSpectralIndexData= hasComponentSpectralIndexData;
				sourceData.isMultiSourceMatchIndex= compSpectralIndexData[k].isMultiSourceMatchIndex;
				sourceData.spectralIndex= compSpectralIndexData[k].spectralIndex;
				sourceData.spectralIndexErr= compSpectralIndexData[k].spectralIndexErr;
				sourceData.isSpectralIndexFit= compSpectralIndexData[k].isSpectralIndexFit;
				sourceData.spectralFitChi2 = compSpectralIndexData[k].spectralFitChi2;
				sourceData.spectralFitNDF = compSpectralIndexData[k].spectralFitNDF;
			}
			else{
				SpectralIndexData dummySid;
				sourceData.hasSpectralIndexData= 0;
				sourceData.isMultiSourceMatchIndex= dummySid.isMultiSourceMatchIndex;
				sourceData.spectralIndex= dummySid.spectralIndex;
				sourceData.spectralIndexErr= dummySid.spectralIndexErr;
				sourceData.isSpectralIndexFit= dummySid.isSpectralIndexFit;
				sourceData.spectralFitChi2 = dummySid.spectralFitChi2;
				sourceData.spectralFitNDF= dummySid.spectralFitNDF;
			}

			//- Astro object cross matches
			if(source->componentObjLocationIds.size()==nComponents){
			  sourceData.objLocationId= source->componentObjLocationIds[k];
			}
			else{
				sourceData.objLocationId= eUNKNOWN_OBJECT_LOCATION;
			}
			if(source->componentObjClassIds.size()==nComponents){
				sourceData.objClassId= source->componentObjClassIds[k];
			}
			else{
				sourceData.objClassId= eUNKNOWN_OBJECT;
			}
			if(source->componentObjClassSubIds.size()==nComponents){
				sourceData.objClassSubId= source->componentObjClassSubIds[k];
			}
			else{
				sourceData.objClassSubId= eUNKNOWN_OBJECT;
			}
			if(source->componentObjConfirmed.size()==nComponents){
				sourceData.objConfirmed= source->componentObjConfirmed[k];
			}
			else{
				sourceData.objConfirmed= 0;
			}

			bool hasComponentAstroObjs= source->HasComponentAstroObjects();
		
			sourceData.hasAstroObjs= hasComponentAstroObjs;

			if(hasComponentAstroObjs && !compAstroObjects.empty() && !compAstroObjects[k].empty()){
				std::stringstream sstream;
				if(compAstroObjects[k].size()==1){
					sstream<<compAstroObjects[k][0].name;
				}
				else{
					for(size_t j=0;j<compAstroObjects[k].size()-1;j++){
						sstream<<compAstroObjects[k][j].name<<";";
					}
					sstream<<compAstroObjects[k][compAstroObjects[k].size()-1].name;		
				}

				sourceData.objName= sstream.str();

			}//close if
			else{
				sourceData.objName= "XXX";
			}

			//## Fill TTree
			dataTree->Fill();

		}//end loop components
	}//close if has fit

	
	//Store nested components
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++){
			FillSourceComponentTree(dataTree,sourceData,nestedSources[k],dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo);
		}
	}//close if has nested sources 
	
	
	//Delete WCS
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	return 0;

}//close FillSourceComponentTree()

//=================================================
//==        DS9 EXPORTER
//=================================================
int SourceExporter::WriteToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS,int ds9WCSType,int ds9RegionFormat,WCS* wcs,std::string ds9RegionColor)
{
	//## Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//## Saving DS9 file region
	std::string ds9WCSTypeHeader= "image";
	if(convertDS9RegionsToWCS) ds9WCSTypeHeader= AstroUtils::GetDS9WCSTypeHeader(ds9WCSType);

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Saving DS9 region header...");
	#endif
	fprintf(fout,"global color=%s font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n",ds9RegionColor.c_str());
	fprintf(fout,"%s\n",ds9WCSTypeHeader.c_str());

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Saving "<<sources.size()<<" sources to file...");
	#endif

	for(size_t k=0;k<sources.size();k++){
		int source_type= sources[k]->Type;
		bool isAtEdge= sources[k]->IsAtEdge();

		//If WCS is not computed, compute it
		if(convertDS9RegionsToWCS && !wcs){
			wcs= sources[k]->GetWCS(ds9WCSType);
			if(!wcs) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to compute WCS from source no "<<k<<"!");
				#endif
			}
		}
	
		//Get DS9 regions
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Dumping DS9 region info for source no. "<<k<<" ...");
		#endif
		std::string regionInfo= "";
		if(ds9RegionFormat==ePolygonRegion) {
			regionInfo= SourceToDS9Region(sources[k],true,convertDS9RegionsToWCS,wcs,ds9WCSType);
		}
		else if(ds9RegionFormat==eEllipseRegion) {
			regionInfo= SourceToDS9EllipseRegion(sources[k],true);
		}
		else {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid DS9RegionType given ("<<ds9RegionFormat<<")");
			#endif
			return -1;
		}

		//Write source region to file
		fprintf(fout,"%s\n",regionInfo.c_str());
	  	
	}//end loop sources
		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Closing DS9 file region...");
	#endif
	fclose(fout);

	return 0;

}//close WriteToDS9()


int SourceExporter::WriteComponentsToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS,int ds9WCSType,WCS* wcs,std::string ds9RegionColor)
{
	//## Open file
	FILE* fout_fit= fopen(filename.c_str(),"w");

	std::string ds9WCSTypeHeader= "image";
	if(convertDS9RegionsToWCS) ds9WCSTypeHeader= AstroUtils::GetDS9WCSTypeHeader(ds9WCSType);

	//## Saving DS9 file region
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Saving DS9 region header for fitted source catalog...");
	#endif
	fprintf(fout_fit,"global color=%s font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n",ds9RegionColor.c_str());
	fprintf(fout_fit,"%s\n",ds9WCSTypeHeader.c_str());

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Saving "<<sources.size()<<" sources to file...");
	#endif
	bool useFWHM= true;

	for(size_t k=0;k<sources.size();k++){
		#ifdef LOGGING_ENABLED		
			DEBUG_LOG("Dumping DS9 region fitting info for source no. "<<k<<" ...");
		#endif

		//If WCS is not computed, compute it
		if(convertDS9RegionsToWCS && !wcs){
			wcs= sources[k]->GetWCS(ds9WCSType);
			if(!wcs) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to compute WCS from source no "<<k<<"!");
				#endif
			}
		}

		//Get DS9 regions for fitted components
		std::string regionInfo= SourceToDS9FittedEllipseRegion(sources[k],useFWHM,true,convertDS9RegionsToWCS,wcs,ds9WCSType);

		if(regionInfo!="") fprintf(fout_fit,"%s\n",regionInfo.c_str());
	}//end loop sources
		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Closing DS9 file region for fitted sources...");
	#endif
	fclose(fout_fit);

	return 0;

}//close WriteComponentsToDS9()

std::string SourceExporter::GetDS9RegionColor(Source* source)
{
	std::string colorStr= "white";
	if(source->Type==eExtended) colorStr= "green";
	else if(source->Type==eCompactPlusExtended) colorStr= "magenta";
	else if(source->Type==ePointLike) colorStr= "red";
	else if(source->Type==eCompact) colorStr= "blue";
	else colorStr= "white";
			
	return colorStr;
		
}//close GetDS9RegionColor()

const std::string SourceExporter::SourceToDS9Region(Source* source,bool dumpNestedSourceInfo,bool convertToWCS,WCS* wcs,int coordSystem)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
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
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to convert contours in WCS, region will be empty!");
			#endif
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
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to convert contours in WCS, region will be empty!");
					#endif
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
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
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


const std::string SourceExporter::SourceToDS9FittedEllipseRegion(Source* source,bool useFWHM,bool dumpNestedSourceInfo,bool convertToWCS,WCS* wcs,int coordSystem)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null input source ptr given!");
		#endif
		return std::string("");
	}

	//Check WCS & metadata
	ImgMetaData* metadata= source->GetImageMetaData();
	if(convertToWCS && !wcs && !metadata){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Requested to convert ellipse to WCS but no wcs was provided and no metadata to compute it are available!");
		#endif
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
	int nSelFitComponents= source->GetNSelFitComponents();

	if(hasFitInfo && nSelFitComponents>0){
		//Get fit pars
		SourceFitPars fitPars= source->GetFitPars();

		//Get fit ellipses
		std::vector<TEllipse*> ellipses;
		
		if(source->GetFitEllipses(ellipses,useFWHM,convertToWCS,wcs,coordSystem,pixOffset)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to get WorldCoord system from metadata!");
			#endif
			return std::string("");
		}

		//Get fit info
		//- fit quality flag
		int fitQuality= fitPars.GetFitQuality();
		std::string fitQualityFlagStr= GetSourceFitQualityStr(fitQuality);
	
		//Loop over fit ellipses and convert to DS9 regions
		for(size_t i=0;i<ellipses.size();i++){
			if(!ellipses[i]) continue;

			//Check if component is selected
			bool isSelected= fitPars.IsSelectedComponent(i);
			if(!isSelected) continue;

			//Get fit component flag
			int fitComponentFlag= -1;
			fitPars.GetComponentFlag(fitComponentFlag,i);
			std::string fitComponentFlagStr= GetSourceFlagStr(fitComponentFlag);

			//Get fit component type
			int fitComponentType= -1;
			fitPars.GetComponentType(fitComponentType,i);
			std::string fitComponentTypeStr= GetSourceTypeStr(fitComponentType);

			//Get encoded string region
			std::string regionText(Form("%s_fitcomp%d",source->GetName(),(int)(i+1)));
			std::string regionColor= "red";
			
			std::vector<std::string> regionTags {fitComponentTypeStr,"fit-component",fitComponentFlagStr,fitQualityFlagStr};
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

