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
class WCS;


struct SourceTreeData 
{
	std::string name;//source name
	std::string iau;//source iau name
	long long nPix;//# pixels
	int nFitComponents;//number of fit components
	int nSelFitComponents;//number of selected fit components
	int nNestedSources;//number of nested sources
	double X0;//pos X0 in image coords
	double Y0;//pos Y0 in image coords
	double X0w;//signal-weighted pos X0 in image coords
	double Y0w;//signal-weighted pos Y0 in image coords
	double X0_wcs;//pos X0 in WCS coords
	double Y0_wcs;//pos Y0 in WCS coords
	double X0w_wcs;//signal-weighted pos X0 in WCS coords
	double Y0w_wcs;//signal-weighted pos Y0 in WCS coords
	double Xmin;//bounding box Xmin
	double Xmax;//bounding box Xmax
	double Ymin;//bounding box Ymin
	double Ymax;//bounding box Ymax
	double Xmin_wcs;//bounding box Xmin in WCS coords
	double Xmax_wcs;//bounding box Xmax in WCS coords
	double Ymin_wcs;//bounding box Ymin in WCS coords
	double Ymax_wcs;//bounding box Ymax in WCS coords
	double Nu;//frequency
	double S;//flux sum over pixels
	double Smax;//flux max
	double fittedFlux;//flux density
	double fittedFluxErr;//flux density err
	double beamArea;//beamArea
	double bkgSum;//bkg sum over pixels
	double rmsSum;;//noise rms sum over pixels
	int type;//source type
	int flag;//source flag
	int good;//source isGoodFlag
	int depthLevel;//source depth level

	double chi2;//Fit chisquare
	double ndf;//Fit number of degrees of freedom
	int fitQuality;//Fit quality flag (eBadFit=0,eLQFit=1,eMQFit=2,eHQFit=3)
	double residualMean;
	double residualRMS;
	double residualMedian;
	double residualMAD;
	double residualMin;
	double residualMax;

	/*
	//Spectral index data
	bool hasSpectralIndexData;
	bool isMultiSourceMatchIndex;
	double spectralIndex;
	double spectralIndexErr;
	bool isSpectralIndexFit;
	double spectralFitChi2;
	double spectralFitNDF;

	//Astro object IDs
	int objLocationId;
	int objClassId;
	int objClassSubId;

	//Astro object cross-matched
	//...
	*/

};//close SourceTreeData


struct SourceComponentTreeData 
{
	std::string name;//Island source name assigned by finder
	long long nPix;//Number of pixels in source
	int componentId;//Fitted component id
	std::string iau;//Fitted component name in IAU notation
	int nFitComponents;//number of fit components in island
	int nSelFitComponents;//number of selected fit components in island
	double X0;//Fitted component centroid in image coordinates along x axis
	double Y0;//Fitted component centroid in image coordinates along y axis
	double X0_err;//Fitted component centroid error in image coordinates along x axis
	double Y0_err;//Fitted component centroid error in image coordinates along y axis
	double X0_wcs;//Fitted component centroid in world coordinates (deg) along x axis
	double Y0_wcs;//Fitted component centroid in world coordinates (deg) along y axis
	double X0_err_wcs;//Fitted component centroid error in world coordinates (deg) along x axis
	double Y0_err_wcs;//Fitted component centroid error in world coordinates (deg) along y axis
	double Nu;//Spectral axis value present in image header. If frequency it is given in GHz units.
	double S;//Sum of island pixel fluxes in Jy/beam units
	double Smax;//Max island pixel flux in Jy/beam units
	double A;//Fitted component amplitude in Jy/beam units (not corrected by beam area)
	double A_err;//Fitted component amplitude error in Jy/beam units (not corrected by beam area)
	double fittedFlux;//Fitted component flux density in Jy/beam units (not corrected by beam area)
	double fittedFlux_err;//Fitted component flux density error in Jy/beam units (not corrected by beam area)
	double fittedIslandFlux;//Fitted source flux density in Jy/beam units (not corrected by beam area)
	double fittedIslandFlux_err;//Fitted source flux density error in Jy/beam units (not corrected by beam area)
	double beamArea;//Number of pixels in beam. Used to convert flux parameters from Jy/beam to Jy/pixel (e.g. Jy/pixel=Jy/beam/beamarea).
	double Bmaj;//Fitted component ellipse major axis in image coordinates
	double Bmin;//Fitted component ellipse major axis in image coordinates
	double Pa;//Fitted component ellipse position angles in deg (measured counterclock-wise from North)
	double Bmaj_err;//Fitted component ellipse major axis error in image coordinates
	double Bmin_err;//Fitted component ellipse major axis error in image coordinates
	double Pa_err;//Fitted component ellipse position angles error in deg
	double Bmaj_wcs;//Fitted component ellipse major axis in world coordinates (arcsec)
	double Bmin_wcs;//Fitted component ellipse major axis in world coordinates (arcsec)
	double Pa_wcs;//Fitted component ellipse position angles in deg (measured counterclock-wise from North);
	double Bmaj_err_wcs;//Fitted component ellipse major axis error in world coordinates (arcsec);
	double Bmin_err_wcs;//Fitted component ellipse major axis error in world coordinates (arcsec);
	double Pa_err_wcs;//Fitted component ellipse position angles error in deg
	double Bmaj_beam;//Beam ellipse major axis (arcsec)
	double Bmin_beam;//Beam ellipse minor axis (arcsec)
	double Pa_beam;//Beam ellipse position angles in deg (measured counterclock-wise from North)
	double Bmaj_deconv_wcs;//Fitted component ellipse major axis in world coordinates, deconvolved by beam (arcsec)
	double Bmin_deconv_wcs;//Fitted component ellipse major axis in world coordinates, deconvolved by beam (arcsec)
	double Pa_deconv_wcs;//Fitted component ellipse position angles in deg, deconvolved by beam (measured counterclock-wise from North)
	double eccentricityRatio;//Eccentricity_fit/Eccentricity_beam - Ratio between eccentricities of fitted and beam ellipses
	double areaRatio;//Area_fit/Area_beam - Ratio between areas of fitted ellipse and beam ellipse 
	double fitVSBeamRotAngle;//Rotation angle in degrees (range 0-180) between fit ellipse and beam ellipse 
	double bkgSum;//Background estimator summed over all source pixels (in Jy/beam)
	double rmsSum;//Noise (rms) estimator summed over all source pixels (in Jy/beam)
	double chi2;//Fit chisquare
	double ndf;//Fit number of degrees of freedom
	int fitQuality;//Fit quality flag (eBadFit=0,eLQFit=1,eMQFit=2,eHQFit=3)
	int flag;//Fitted component flag (eReal=1,eCandidate=2,eFake=3)
	int type;//Fitted component type (eUnknown=0,eCompact=1,ePoint-Like=2,eExtended=3)
	
};//close SourceComponentTreeData()


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
		* \brief Write source collection to ascii file
		*/
		static int WriteToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0,bool writeAdditionalSourceInfo=false,char delimiter='\t');
		
		/**
		* \brief Get source ascii string
		*/
		static const std::vector<std::string> SourceToAscii(Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0,bool writeAdditionalSourceInfo=false,char delimiter='\t');

		/**
		* \brief Write ascii file from source component collection
		*/
		static int WriteComponentsToAscii(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0,bool writeAdditionalSourceInfo=false,char delimiter='\t');
		
		/**
		* \brief Get source component ascii string
		*/
		static const std::vector<std::string> SourceComponentsToAscii(Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0,bool writeAdditionalSourceInfo=false,char delimiter='\t');
	
		//=======================================
		//==      JSON EXPORT
		//=======================================
		
		/**
		* \brief Write json file from source collection
		*/
		static int WriteToJson(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		
		/**
		* \brief Get source json object
		*/
		static int SourceToJson(std::vector<Json::Value>& jsonValues,Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);

		/**
		* \brief Write json file from source component collection
		*/
		static int WriteComponentsToJson(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		
		/**
		* \brief Get source component json object
		*/
		static int SourceComponentsToJson(std::vector<Json::Value>& jsonValues,Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
	
		//=======================================
		//==      ROOT EXPORT
		//=======================================
		/**
		* \brief Write ROOT file with TTree from source collection
		*/
		static int WriteToROOT(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		/**
		* \brief Write ROOT file with TTree from source fit component collection
		*/
		static int WriteComponentsToROOT(std::string filename,const std::vector<Source*>& sources,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		

		//=======================================
		//==      DS9 EXPORT
		//=======================================
		/**
		* \brief Write DS9 regions from source collection
		*/
		static int WriteToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS=false,int ds9WCSType=eJ2000,int ds9RegionFormat=ePolygonRegion,WCS* wcs=0);

		/**
		* \brief Write DS9 regions for source fitted components
		*/
		static int WriteComponentsToDS9(std::string filename,const std::vector<Source*>& sources,bool convertDS9RegionsToWCS=false,int ds9WCSType=eJ2000,WCS* wcs=0);


		/**
		* \brief Get DS9 region info
		*/
		static const std::string SourceToDS9Region(Source* source,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1);

		/**
		* \brief Get DS9 ellipse info
		*/
		static const std::string SourceToDS9EllipseRegion(Source* source,bool dumpNestedSourceInfo=false);
	
		/**
		* \brief Get DS9 fitted ellipse info
		*/
		static const std::string SourceToDS9FittedEllipseRegion(Source* source,bool useFWHM=true,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1);

	
	private:

		/**
		* \brief Get DS9 region color according to source type
		*/
		static std::string GetDS9RegionColor(Source* source);

		/**
		* \brief Fill ROOT TTree with source info
		*/
		static int FillSourceTTree(TTree* dataTree,SourceTreeData& sourceTreeData,Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
		/**
		* \brief Fill ROOT TTree with source fit component info
		*/
		static int FillSourceComponentTree(TTree* dataTree,SourceComponentTreeData& sourceData,Source* source,bool dumpNestedSourceInfo=true,int wcsType=eJ2000,WCS* wcs=0);
	

	private:
	
		ClassDef(SourceExporter,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class SourceExporter+;
#endif

}//close namespace 


#endif
