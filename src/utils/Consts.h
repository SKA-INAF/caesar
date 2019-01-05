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
* @file Consts.h
* @class Consts
* @brief Definitions & constants
*
* Definitions & constants
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef _CONSTS_h
#define _CONSTS_h 1

#include <math.h>
#include <string>

namespace Caesar {

//====================================
//===          ENUM              =====
//====================================
enum WCSType {
	eJ2000= 0,
	eB1950= 1,
	eGALACTIC= 2,
	eECLIPTIC= 3,
	eALTAZ= 4,
	eLINEAR= 5
};

enum ImgFilters {
	eGausFilter= 1,
	eGuidedFilter= 2,
	eWaveletFilter= 3,
	eLoGFilter= 4
};

enum ImgType {
	//eUnknown=0,
	eUnknownMap=0,
	eFluxMap=1,
	eMeanFluxMap=2,
	eSignificanceMap=3,
	eBkgMap=4,
	eNoiseMap=5,
	eBinaryMap=6,
	eResidualMap=7,
	ePullMap=8,
	eCurvatureMap=9
};

enum ImgSmoothFilter {
	eGaus=1,
	eGuided=2,
	eWT=3,
	eGradient=4,
	eLaplacian=5,
	eLoG=6,	
	eNormLoG=7
};

enum SegmAlgo {
	eWaveletTransform= 1,
	eHClust= 2,
	eActiveContour= 3,
	eSaliencyThr= 4,
};

enum BlobMaskMethod {
	eCurvMask=1,
	eMultiScaleLoGMask=2
};

enum BkgEstimator {
	eMeanBkg=1,
	eMedianBkg=2,
	eBiWeightBkg=3,
	eMedianClippedBkg= 4
};
		
enum BkgMethod {
	eGridBkg=1,
	eSuperpixelBkg=2
};

/**
* \brief Source type enumeration
*/
enum SourceType {
	eUnknownType=0,
	eCompact=1,
	ePointLike=2,
	eExtended=3,
	eCompactPlusExtended=4
};

/**
* \brief Source flag enumeration
*/
enum SourceFlag {
	eReal=1,
	eCandidate=2,
	eFake=3
};

/**
* \brief Convert source flag enumeration to string
*/
inline std::string GetSourceFlagStr(int sourceFlag)
{
	std::string flagStr= "";
	if(sourceFlag==eReal) flagStr= "real";
	else if(sourceFlag==eCandidate) flagStr= "candidate";
	else if(sourceFlag==eFake) flagStr= "fake";
	else flagStr= "unknown";
	return flagStr;
}



/**
* \brief Simulated source type enumeration
*/
enum SimSourceType {
	eUnknownSimClass=0,
	eRingLike=1,
	eBubbleLike=2,
	eEllipseLike=3,
	eDiskLike=4,	
	eBlobLike=5
};

/**
* \brief Source fit quality enumeration
*/
enum SourceFitQuality {
	eBadFit=0,//not converged
	eLQFit=1,//converged with pars at boundary, error matrix not posdef
	eMQFit=2,//converged with good error estimates, not passing quality selection cuts
	eHQFit=3//passing quality selection
};

/**
* \brief Convert source fit quality enumeration to string
*/
inline std::string GetSourceFitQualityStr(int fitQuality)
{
	std::string flagStr= "";
	if(fitQuality==eBadFit) flagStr= "bad-fit";
	else if(fitQuality==eLQFit) flagStr= "lq-fit";
	else if(fitQuality==eMQFit) flagStr= "mq-fit";
	else if(fitQuality==eHQFit) flagStr= "hq-fit";
	else flagStr= "uq-fit";//unknown
	return flagStr;
}
		

enum ColorPaletteStyle {
	eRAINBOW= 0,
	eBLACKWHITE= 1,
	eBLACKBODY= 2,
	eHOT2COLD= 3,
	eCOLD2HOT= 4,
	eTHERMAL= 5
};

enum FileType{
	eROOT= 0,
	eFITS= 1
};

enum SLICEdgeModel {
	eKirschEdge= 1,
	eChanVeseEdge= 2
};

enum activeContourMethod {
	eChanVeseAC= 1,
	eLRAC= 2
};

enum activeContourInitLevelSetMethod {
	eCircleLevelSet= 1,
	eCheckerboardLevelSet= 2,
	eSaliencyLevelSet= 3
};

enum DS9RegionFormat {
	eEllipseRegion= 1,
	ePolygonRegion= 2
};

enum MorphStructElement {
	eMORPH_RECT= 0,
	eMORPH_ELLIPSE= 1,
	eMORPH_CROSS= 2
};

enum MorphOperation {
	eMORPH_OPENING= 0,
	eMORPH_CLOSING= 1,
	eMORPH_GRADIENT= 2,
	eMORPH_TOPHAT= 3,
	eMORPH_BLACKHAT= 4,
	eMORPH_EROSION= 5,
	eMORPH_DILATION= 6
};

enum DilationModel {
	eDilateWithBkg=1,
	eDilateWithSourceMedian=2
};

enum PSSubtractionMethod {
	ePS_DILATION= 1,
	ePS_MODELSUBTRACTION= 2
};


enum ContourOverlapFlag {
	eCONT_OVERLAP_FAILED=-1,
	eCONT_NOT_OVERLAPPING=0,
	eCONT_OVERLAPPING=1,
	eCONT1_INSIDE_CONT2=2,
	eCONT2_INSIDE_CONT1=3
};

enum FitMinimizer {
	eMIGRAD= 1,
	eHESS= 2,
	eMINOS= 3,
	eSIMPLEX= 4,	
	eMINIMIZE= 5
};

#ifdef __MAKECINT__
#pragma link C++ enum WCSType+;
#pragma link C++ enum ImgFilters+;
#pragma link C++ enum ImgType+;
#pragma link C++ enum ImgSmoothFilter+;
#pragma link C++ enum SegmAlgo+;
#pragma link C++ enum BlobMaskMethod+;
#pragma link C++ enum BkgEstimator+;
#pragma link C++ enum BkgMethod+;
#pragma link C++ enum ColorPaletteStyle+;
#pragma link C++ enum FileType+;
#pragma link C++ enum SLICEdgeModel+;
#pragma link C++ enum activeContourMethod+;
#pragma link C++ enum DS9RegionFormat+;
#pragma link C++ enum MorphStructElement+;
#pragma link C++ enum MorphOperation+;
#pragma link C++ enum ContourOverlapFlag+;
#pragma link C++ enum FitMinimizer+;
#endif

//====================================
//===          CONSTANTS         =====
//====================================
static const double GausSigma2FWHM= 2.*sqrt(2*log(2.));



}//close namespace


#endif 
