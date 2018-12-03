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
//==        ASCII EXPORTER
//=================================================
int SourceExporter::WriteToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WorldCoor* wcs)
{
	//Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Write header
	std::stringstream ss;
	ss<<"#======== CAESAR RESULTS ======\n";
	ss<<"# This ascii file contains source islands (or blobs) found by CAESAR.\n";
	ss<<"# Source parameters are described in the header below.\n";
	ss<<"# NB: Fitted components, if available, are reported in a different ascii file. \n";
	ss<<"#     We include here the total fitted flux density summed over all components.\n";
	ss<<"#\n";
	ss<<"# ---- HEADER ----\n";
	ss<<"# name - Source name assigned by finder\n";
	ss<<"# IAU name - Source name in IAU notation\n";
	ss<<"# Npix - Number of pixels in source\n";
	ss<<"# nFittedComponents - Number of fitted components (=0 if fit not performed or failed)\n";
	ss<<"# nNestedSources - Number of nested sources detected\n";
	ss<<"# X0 - Source centroid in image coordinates along x axis \n";
	ss<<"# Y0 - Source centroid in image coordinates along y axis \n";
	ss<<"# X0w - Source centroid in image coordinates along x axis, weighted by pixel fluxes \n";
	ss<<"# Y0w - Source centroid in image coordinates along y axis, weighted by pixel fluxes \n";
	ss<<"# X0_wcs - Source centroid in world coordinates (deg) along x axis\n";
	ss<<"# Y0_wcs - Source centroid in world coordinates (deg) along y axis\n";
	ss<<"# X0w_wcs - Source centroid in world coordinates (deg) along x axis, weighted by pixel fluxes\n";
	ss<<"# Y0w_wcs - Source centroid in world coordinates (deg) along y axis, weighted by pixel fluxes\n";
	ss<<"# Xmin - Source minimum pixel image coordinate along x axis\n";
	ss<<"# Xmax - Source maximum pixel image coordinate along x axis\n";
	ss<<"# Ymin - Source minimum pixel image coordinate along y axis\n";
	ss<<"# Ymax - Source maximum pixel image coordinate along y axis\n";
	ss<<"# Xmin_wcs - Source minimum pixel WCS coordinate (deg) along x axis\n";
	ss<<"# Xmax_wcs - Source maximum pixel WCS coordinate (deg) along x axis\n";
	ss<<"# Ymin_wcs - Source minimum pixel WCS coordinate (deg) along y axis\n";
	ss<<"# Ymax_wcs - Source maximum pixel WCS coordinate (deg) along y axis\n";
	ss<<"# Nu - Spectral axis value present in image header. If frequency it is given in GHz units.\n";
	//ss<<"# dNu - Spectral axis width present in image header. If frequency it is given in GHz units.\n";
	ss<<"# S - Sum of source pixel fluxes in Jy/beam units.\n";
	ss<<"# Smax - Maximum source pixel flux in Jy/beam units.\n";
	ss<<"# FluxDensity - Fitted source flux density in Jy/beam units (not corrected by beam area).\n";
	ss<<"# FluxDensityErr - Fitted source flux density error in Jy/beam units (not corrected by beam area).\n";
	ss<<"# BeamArea- Number of pixels in beam. Used to convert flux parameters from Jy/beam to Jy/pixel (e.g. Jy/pixel=Jy/beam/beamarea).\n";
	ss<<"# BkgSum- Background estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# RMSSum- Noise (rms) estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# Type- Source type tag (eUnknown=0,eCompact=1,ePointLike=2,eExtended=3,eCompactPlusExtended=4)\n";
	ss<<"# Flag- Source flag (eReal=1,eCandidate=2,eFake=3)\n";
	ss<<"# IsGoodSource- Bool flag indicating if source was tagged as good (true) or bad (false) in the finding process\n";
	ss<<"# DepthLevel- Source depth level (0=mother source,1=nested source,...)\n";
	ss<<"# ---------------\n";
	ss<<"#\n";
	fprintf(fout,"%s",ss.str().c_str());

	//Loop sources
	bool deleteWCS= false;

	for(size_t k=0;k<sources.size();k++){
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				WARN_LOG("No metadata are available to retrieve WCS!");
				return -1;
			}
			wcs= metadata->GetWorldCoord(wcsType);
			if(!wcs){
				ERROR_LOG("Failed to get WCS from metadata!");
				return -1;
			}
			deleteWCS= true;
		}

		//Print source info to ascii
		/*
		std::string s= SourceToAscii(sources[k],wcsType,wcs);
		fprintf(fout,"%s\n",s.c_str());
		
		//Add nested sources
		if(dumpNestedSourceInfo){
			std::vector<Source*> nestedSources= sources[k]->GetNestedSources();
			for(size_t j=0;j<nestedSources.size();j++){
				std::string s_nested= SourceToAscii(nestedSources[j],wcsType,wcs);
				fprintf(fout,"%s\n",s_nested.c_str());
			}//end loop sources
		}//close if
		*/

		std::vector<std::string> slist= SourceToAscii(sources[k],dumpNestedSourceInfo,wcsType,wcs);
		for(size_t j=0;j<slist.size();j++){
			fprintf(fout,"%s\n",slist[j].c_str());
		}
		
	}//end loop sources

	//Delete WCS
	if(deleteWCS) CodeUtils::DeletePtr<WorldCoor>(wcs);

	//Close file
	fclose(fout);

	return 0;

}//close WriteToAscii()


const std::vector<std::string> SourceExporter::SourceToAscii(Source* source,bool dumpNestedSourceInfo,int wcsType,WorldCoor* wcs)
{
	//Init string list
	std::vector<std::string> sourceStrList;

	//Check source
	if(!source){
		WARN_LOG("Null input source ptr given!");
		return sourceStrList;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWorldCoord(wcsType);
			if(wcs) deleteWCS= true;
			else WARN_LOG("Failed to get WCS from metadata!");
		}
		else {
			WARN_LOG("No metadata are available to retrieve WCS!");
		}		
	}
	

	//Compute IAU name
	bool useWeightedPos= false;
	std::string iauName= source->GetName();
	if(wcs) source->GetIAUName(useWeightedPos,wcs,wcsType);

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
	std::stringstream ss;
	
	//- Source name
	ss<<source->GetName()<<"\t";
	ss<<iauName<<"\t";

	
	//- Number of pixels
	ss<<source->NPix<<"\t";

	//- Number of sub-components
	ss<<source->GetNFitComponents()<<"\t";
	ss<<source->GetNestedSourceNumber()<<"\t";

	//- Pixel centroids
	ss<<source->X0<<"\t"<<source->Y0<<"\t";
	ss<<source->GetSx()<<"\t"<<source->GetSy()<<"\t";
	ss<<X0_wcs<<"\t"<<Y0_wcs<<"\t"<<X0_weighted_wcs<<"\t"<<Y0_weighted_wcs<<"\t";

	//- Bounding box
	ss<<xmin<<"\t"<<xmax<<"\t"<<ymin<<"\t"<<ymax<<"\t";
	ss<<xmin_wcs<<"\t"<<xmax_wcs<<"\t"<<ymin_wcs<<"\t"<<ymax_wcs<<"\t";
	
	//- Frequency
	ss<<Nu<<"\t";

	//- Flux 
	ss<<source->GetS()<<"\t";
	ss<<source->GetSmax()<<"\t";
	ss<<fluxDensity<<"\t";
	ss<<fluxDensityErr<<"\t";
	ss<<source->GetBeamFluxIntegral()<<"\t";

	//Bkg/noise estimators
	ss<<source->GetBkgSum()<<"\t";
	ss<<source->GetBkgRMSSum()<<"\t";

	//- Source flags
	ss<<source->Type<<"\t"<<source->Flag<<"\t"<<source->IsGoodSource()<<"\t";
	ss<<source->GetDepthLevel();
	//ss<<source->HasFitInfo()<<"\t";
	//ss<<source->HasNestedSources();

	//Add string to list
	sourceStrList.push_back(ss.str());

	//Store nested sources
	bool hasNestedSources= source->HasNestedSources();
	if(hasNestedSources && dumpNestedSourceInfo){
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++)
		{
			std::vector<std::string> sourceStrList_nested= SourceToAscii(nestedSources[k],dumpNestedSourceInfo,wcsType,wcs);
			if(!sourceStrList_nested.empty()){
				sourceStrList.insert(sourceStrList.end(),sourceStrList_nested.begin(),sourceStrList_nested.end());
			}
		}
	}//close if has nested sources 
		

	//Delete WCS
	if(deleteWCS) CodeUtils::DeletePtr<WorldCoor>(wcs);

	return sourceStrList;

}//close SourceToAscii()


int SourceExporter::WriteComponentsToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo,int wcsType,WorldCoor* wcs)
{
	//Open output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Write header
	std::stringstream ss;
	ss<<"#======== CAESAR RESULTS ======\n";
	ss<<"# This ascii file contains source components fitted to islands/blobs detected by CAESAR.\n";
	ss<<"# Source parameters are described in the header below.\n";
	ss<<"#\n";
	ss<<"# ---- HEADER ----\n";
	ss<<"# name - Island source name assigned by finder\n";
	ss<<"# Npix - Number of pixels in source\n";
	ss<<"# componentId - Fitted component id\n";
	ss<<"# IAU name - Fitted component name in IAU notation\n";
	ss<<"# X0 - Fitted component centroid in image coordinates along x axis \n";
	ss<<"# Y0 - Fitted component centroid in image coordinates along y axis \n";
	ss<<"# X0_err - Fitted component centroid error in image coordinates along x axis \n";
	ss<<"# Y0_err - Fitted component centroid error in image coordinates along y axis \n";
	ss<<"# X0_wcs - Fitted component centroid in world coordinates (deg) along x axis\n";
	ss<<"# Y0_wcs - Fitted component centroid in world coordinates (deg) along y axis\n";
	ss<<"# X0_err_wcs - Fitted component centroid error in world coordinates (deg) along x axis\n";
	ss<<"# Y0_err_wcs - Fitted component centroid error in world coordinates (deg) along y axis\n";
	ss<<"# Nu - Spectral axis value present in image header. If frequency it is given in GHz units.\n";
	
	ss<<"# A - Fitted component amplitude in Jy/beam units (not corrected by beam area).\n";
	ss<<"# A_err - Fitted component amplitude error in Jy/beam units (not corrected by beam area).\n";
	ss<<"# FluxDensity - Fitted component flux density in Jy/beam units (not corrected by beam area).\n";
	ss<<"# FluxDensity_err - Fitted component flux density error in Jy/beam units (not corrected by beam area).\n";
	ss<<"# Island FluxDensity - Fitted source flux density in Jy/beam units (not corrected by beam area).\n";
	ss<<"# Island FluxDensity_err - Fitted source flux density error in Jy/beam units (not corrected by beam area).\n";
	ss<<"# BeamArea- Number of pixels in beam. Used to convert flux parameters from Jy/beam to Jy/pixel (e.g. Jy/pixel=Jy/beam/beamarea).\n";
	
	ss<<"# Bmaj - Fitted component ellipse major axis in image coordinates\n";
	ss<<"# Bmin - Fitted component ellipse major axis in image coordinates\n";
	ss<<"# Pa - Fitted component ellipse position angles in deg (measured counterclock-wise from North)\n";
	ss<<"# Bmaj_err - Fitted component ellipse major axis error in image coordinates\n";
	ss<<"# Bmin_err - Fitted component ellipse major axis error in image coordinates\n";
	ss<<"# Pa_err - Fitted component ellipse position angles error in deg\n";
	ss<<"# Bmaj_wcs - Fitted component ellipse major axis in world coordinates (arcsec)\n";
	ss<<"# Bmin_wcs - Fitted component ellipse major axis in world coordinates (arcsec)\n";
	ss<<"# Pa_wcs - Fitted component ellipse position angles in deg (measured counterclock-wise from North)\n";
	ss<<"# Bmaj_wcs_err - Fitted component ellipse major axis error in world coordinates (arcsec)\n";
	ss<<"# Bmin_wcs_err - Fitted component ellipse major axis error in world coordinates (arcsec)\n";
	ss<<"# Pa_wcs_err - Fitted component ellipse position angles error in deg\n";
	ss<<"# Bmaj_beam - Beam ellipse major axis (arcsec)\n";
	ss<<"# Bmin_beam - Beam ellipse minor axis (arcsec)\n";
	ss<<"# Pa_beam - Beam ellipse position angles in deg (measured counterclock-wise from North)\n";
	ss<<"# Bmaj_deconv_wcs - Fitted component ellipse major axis in world coordinates, deconvolved by beam (arcsec)\n";
	ss<<"# Bmin_deconv_wcs - Fitted component ellipse major axis in world coordinates, deconvolved by beam (arcsec)\n";
	ss<<"# Pa_deconv_wcs - Fitted component ellipse position angles in deg, deconvolved by beam (measured counterclock-wise from North)\n";
	ss<<"# BkgSum- Background estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# RMSSum- Noise (rms) estimator summed over all source pixels (in Jy/beam).\n";
	ss<<"# Chi2- Fit chisquare.\n";
	ss<<"# NDF - Fit number of degrees of freedom.\n";
	ss<<"# ---------------\n";
	ss<<"#\n";
	fprintf(fout,"%s",ss.str().c_str());

	//Loop sources
	bool deleteWCS= false;

	for(size_t k=0;k<sources.size();k++){
		//If wcs is not given, retrieve it from metadata
		if(!wcs){
			ImgMetaData* metadata= sources[k]->GetImageMetaData();
			if(!metadata){
				WARN_LOG("No metadata are available to retrieve WCS!");
				return -1;
			}
			wcs= metadata->GetWorldCoord(wcsType);
			if(!wcs){
				ERROR_LOG("Failed to get WCS from metadata!");
				return -1;
			}
			deleteWCS= true;
		}

		//Print source info to ascii
		std::vector<std::string> slist= SourceComponentsToAscii(sources[k],dumpNestedSourceInfo,wcsType,wcs);
		for(size_t j=0;j<slist.size();j++){
			fprintf(fout,"%s\n",slist[j].c_str());
		}
		
	}//end loop sources

	//Delete WCS
	if(deleteWCS) CodeUtils::DeletePtr<WorldCoor>(wcs);

	//Close file
	fclose(fout);

	return 0;

}//close WriteComponentsToAscii()


const std::vector<std::string> SourceExporter::SourceComponentsToAscii(Source* source,bool dumpNestedSourceInfo,int wcsType,WorldCoor* wcs)
{
	//Init vector
	std::vector<std::string> fitComponentStrList;
	
	//Check source
	if(!source){
		WARN_LOG("Null input source ptr given!");
		return fitComponentStrList;
	}

	//Get metadata
	ImgMetaData* metadata= source->GetImageMetaData();
		
	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(metadata){
			wcs= metadata->GetWorldCoord(wcsType);
			if(wcs) deleteWCS= true;
			else WARN_LOG("Failed to get WCS from metadata!");
		}
		else {
			WARN_LOG("No metadata are available to retrieve WCS!");
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
				WARN_LOG("Failed to retrieve ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
			}
		
			double x0_err= 0;
			double y0_err= 0;
			double bmaj_err= 0;
			double bmin_err= 0;
			double pa_err= 0;
			if(fitPars.GetComponentFitEllipseParErrors(k,x0_err,y0_err,bmaj_err,bmin_err,pa_err)<0){
				WARN_LOG("Failed to retrieve ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
			}

			//- WCS fit ellipse pars
			double x0_wcs= 0;
			double y0_wcs= 0;
			double bmaj_wcs= 0;
			double bmin_wcs= 0;
			double pa_wcs= 0;
			if(fitPars.GetComponentFitWCSEllipsePars(k,x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs)<0){
				WARN_LOG("Failed to retrieve WCS ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
			}
			
			double x0_wcs_err= 0;
			double y0_wcs_err= 0;
			double bmaj_wcs_err= 0;
			double bmin_wcs_err= 0;
			double pa_wcs_err= 0;
			if(fitPars.GetComponentFitWCSEllipseParErrors(k,x0_wcs_err,y0_wcs_err,bmaj_wcs_err,bmin_wcs_err,pa_wcs_err)<0){
				WARN_LOG("Failed to retrieve WCS ellipse par errors for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
			}

			//- Beam ellipse pars
			double bmaj_beam= 0;
			double bmin_beam= 0;
			double pa_beam= 0;
			if(fitPars.GetComponentBeamEllipsePars(k,bmaj_beam,bmin_beam,pa_beam)<0){
				WARN_LOG("Failed to retrieve beam ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
			}

			//- WCS beam-deconvolved ellipse pars
			double bmaj_deconv_wcs= 0;
			double bmin_deconv_wcs= 0;
			double pa_deconv_wcs= 0;
			if(fitPars.GetComponentFitWCSDeconvolvedEllipsePars(k,bmaj_deconv_wcs,bmin_deconv_wcs,pa_deconv_wcs)<0){
				WARN_LOG("Failed to retrieve WCS beam-deconvolved ellipse pars for component no. "<<k+1<<" (hint: check if they are computed correctly), setting dummy values!");
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
					WARN_LOG("Failed to compute IAU name for component no. "<<k+1<<", assume default name "<<iau_default<<" ...");
					iau= iau_default;		
				}
			}

			//## Fill source component
			std::stringstream ss;

			//- Mother source name & npix
			ss<<source->GetName()<<"\t";
			ss<<source->NPix<<"\t";

			//- Component id
			ss<<k+1<<"\t";

			//- Component IAU name
			ss<<iau<<"\t";
				
			//- Position in pixel coordinates
			ss<<x0<<"\t"<<y0<<"\t";
			ss<<x0_err<<"\t"<<y0_err<<"\t";

			//- Position in WCS coordinates
			ss<<x0_wcs<<"\t"<<y0_wcs<<"\t";
			ss<<x0_wcs_err<<"\t"<<y0_wcs_err<<"\t";

			//- Frequency
			ss<<Nu<<"\t";
		
			//- Flux amplitude & density
			ss<<A<<"\t"<<A_err<<"\t";
			ss<<fluxDensity<<"\t"<<fluxDensity_err<<"\t";
			ss<<fluxDensityTot<<"\t"<<fluxDensityTot_err<<"\t";
			ss<<source->GetBeamFluxIntegral()<<"\t";
	
			//- Ellipse pars in pixel coordinates
			ss<<bmaj<<"\t"<<bmin<<"\t"<<pa<<"\t";
			ss<<bmaj_err<<"\t"<<bmin_err<<"\t"<<pa_err<<"\t";
				
			//- Ellipse pars in WCS coordinates
			ss<<bmaj_wcs<<"\t"<<bmin_wcs<<"\t"<<pa_wcs<<"\t";
			ss<<bmaj_wcs_err<<"\t"<<bmin_wcs_err<<"\t"<<pa_wcs_err<<"\t";
		
			//Beam ellipse pars
			ss<<bmaj_beam<<"\t"<<bmin_beam<<"\t"<<pa_beam<<"\t";

			//- Beam-deconvolved ellipse pars in WCS coordinates
			ss<<bmaj_deconv_wcs<<"\t"<<bmin_deconv_wcs<<"\t"<<pa_deconv_wcs<<"\t";
			//ss<<bmaj_deconv_wcs_err<<"\t"<<bmin_deconv_wcs_err<<"\t"<<pa_deconv_wcs_err<<"\t";
		
			//Bkg/noise estimators
			ss<<source->GetBkgSum()<<"\t";
			ss<<source->GetBkgRMSSum()<<"\t";

			//Fit chi2/ndf
			ss<<fitPars.GetChi2()<<"\t"<<fitPars.GetNDF();

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
			std::vector<std::string> fitComponentStrList_nested= SourceComponentsToAscii(nestedSources[k],dumpNestedSourceInfo,wcsType,wcs);
			if(!fitComponentStrList_nested.empty()){
				fitComponentStrList.insert(fitComponentStrList.end(),fitComponentStrList_nested.begin(),fitComponentStrList_nested.end());
			}
		}
	}//close if has nested sources 
	
	
	//Delete WCS
	if(deleteWCS) CodeUtils::DeletePtr<WorldCoor>(wcs);

	return fitComponentStrList;

}//close SourceComponentsToAscii()

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
