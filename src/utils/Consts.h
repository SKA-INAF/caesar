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

#include <TMath.h>
#include <math.h>
#include <string>

namespace Caesar {

//====================================
//===          ENUM              =====
//====================================
enum WCSType {
	eUNKNOWN_CS= -1,
	eJ2000= 0,
	eB1950= 1,
	eGALACTIC= 2,
	eECLIPTIC= 3,
	eALTAZ= 4,
	eLINEAR= 5,
	eIMG_CS= 6,
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
  eUnknownSourceFlag=-1,
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
	else if(sourceFlag==eUnknownSourceFlag) flagStr= "unknown-flag";
	else flagStr= "unknown-flag";
	return flagStr;
}

/**
* \brief Convert source flag enumeration from string
*/
inline int GetSourceFlag(std::string flagStr)
{
	int sourceFlag= eUnknownSourceFlag;
	if(flagStr=="real") sourceFlag= eReal;
	else if(flagStr=="candidate") sourceFlag= eCandidate;
	else if(flagStr=="fake") sourceFlag= eFake;
	else if(flagStr=="unknown-flag") sourceFlag= eUnknownSourceFlag;
	else sourceFlag= eUnknownSourceFlag;
	return sourceFlag;
}

/**
* \brief Convert source type enumeration to string
*/
inline std::string GetSourceTypeStr(int sourceType)
{
	std::string typeStr= "";
	if(sourceType==eUnknownType) typeStr= "unknown-type";
	else if(sourceType==eCompact) typeStr= "compact";
	else if(sourceType==ePointLike) typeStr= "point-like";
	else if(sourceType==eExtended) typeStr= "extended";
	else if(sourceType==eCompactPlusExtended) typeStr= "compact-extended";
	else typeStr= "unknown-type";
	return typeStr;
}

/**
* \brief Convert source type string to enumeration
*/
inline int GetSourceType(std::string typeStr)
{
	int type= eUnknownType;
	if(typeStr=="" || typeStr=="unknown-type") type= eUnknownType;
	else if(typeStr=="compact") type= eCompact;
	else if(typeStr=="point-like") type= ePointLike;
	else if(typeStr=="extended") type= eExtended;
	else if(typeStr=="compact-extended") type= eCompactPlusExtended;
	else type= eUnknownType;
	return type;
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
enum SourceFitQuality 
{
	eUnknownFitQuality=-1,//unknown
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

/**
* \brief Convert source fit quality enumeration to string
*/
inline int GetSourceFitQuality(std::string flagStr)
{
	int fitQuality= eUnknownFitQuality;
	if(flagStr=="bad-fit") fitQuality= eBadFit;
	else if(flagStr=="lq-fit") fitQuality= eLQFit;
	else if(flagStr=="mq-fit") fitQuality= eMQFit;
	else if(flagStr=="hq-fit") fitQuality= eHQFit;
	else if(flagStr=="uq-fit") fitQuality= eUnknownFitQuality;
	else fitQuality= eUnknownFitQuality;
	return fitQuality;
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

/**
* \brief Fit status enum flag
*/
enum SourceFitStatus {
	eFitUnknownStatus= 0,
	eFitAborted= 1,
	eFitNotConverged= 2,
	eFitConverged= 3,
	eFitConvergedWithWarns= 4
};

/**
* \brief Logger target enumerations
*/
enum LoggerTarget {
	eCONSOLE_TARGET= 1,
	eFILE_TARGET= 2,
	eSYSLOG_TARGET= 3
};

enum CutType {
	eUNKNOWN_CUT=0,
	eEQUALITY_CUT=1,
	eBOUND_CUT=2,
	eSINGLE_BOUND_CUT=3
};



enum AstroObjectType 
{
	//Major categories
	eMULTI_CLASS_OBJECT=-1,
	eUNKNOWN_OBJECT=0,
	eSTAR=1,
	eGALAXY=2,
	eGALAXY_C1=2001,
	eGALAXY_C2=2002,
	eGALAXY_C3=2003,
	eGALAXY_C4=2004,
	eGALAXY_C5=2005,
	eGALAXY_C6=2006,
	eSTARFORMING_GALAXY=3000,
	eSTARBURST_GALAXY=4000,
	eAGN=5000,
	eQSO=6000,	
	ePN=3,
	eSNR=4,
	eBUBBLE=5,
	eHII=6,
	eMOLECULAR_CLOUD=7,
	eGALAXY_CLUSTER=8,
	eGALAXY_GROUP=9,
	eGALAXY_SUPERCLUSTER=10,
	eGLOBULAR_CLUSTER=11,
	ePLANET=12,
	eSTAR_CLUSTER=13,
	eX_BINARY=14,
	eNEBULA=15,
	eSTAR_FORMING_REGION=16,
	eCLOUD=17,
	eGAMMA_RAY_BURST=18,
	eSTAR_BINARY=19,
	eHI_SHELL=20,
	eNOVA=21,
	eTRANSIENT_EVENT=22,	
	ePULSAR=23,
	eYSO=24,
	
	//Generic unclassified objects by wavelength
	eRADIO_OBJ=100,
	eIR_OBJ=200,
	eX_OBJ=300,
	eUV_OBJ=400,
	eGAMMA_OBJ=500,

	//Stars sub-types
	eSTAR_IN_CLUSTER=1001,
	eSTAR_IN_NEBULA=1002,
	eSTAR_IN_ASSOCIATION=1003,
	eSTAR_IN_DOUBLE_SYSTEM=1004,
	eSTAR_VARIABLE=1005,
	eSTAR_PECULIAR=1006,
	eSTAR_HB=1007,
	eSTAR_YSO=1008,
	eSTAR_HERBIG=1009,
	eSTAR_EMISS_LINE=1010,
	eSTAR_BE=1011,
	eSTAR_BS=1012,
	eSTAR_RG=1013,
	eSTAR_AB=1014,
	eSTAR_C=1015,
	eSTAR_S=1016,
	eSTAR_ESG=1017,
	eSTAR_RSG=1018,
	eSTAR_YSG=1019,
	eSTAR_BSG=1020,
	eSTAR_HSD=1021,
	eSTAR_PAGB=1022,
	eSTAR_WD=1023,
	eSTAR_PWD=1024,
	eSTAR_LM=1025,
	eSTAR_BD=1026,
	eSTAR_NEUTRON=1027,
	eSTAR_OH=1028,
	eSTAR_CH=1029,
	eSTAR_PMS=1030,
	eSTAR_TTAU=1031,
	eSTAR_WR=1032,
	eSTAR_PM=1033,
	eSTAR_HV=1034,
	eSTAR_FLARE=1035,
	eSTAR_PULSAR=1036,
	eSTAR_SUPERNOVA=1037,
	eSTAR_SUBSTELLAR=1038,
	
};//close AstroObjectType enum codes

enum AstroObjectLocation 
{
	eUNKNOWN_OBJECT_LOCATION=0,
	eGALACTIC_OBJECT=1,
	eEXTRAGALACTIC_OBJECT=2
};


/**
* \brief Convert astro object type enumeration to string
*/
inline int GetAstroObjectType(std::string typeStr)
{
	int type= eUNKNOWN_OBJECT;
	if(typeStr=="STAR") type= eSTAR;	
	else if(typeStr=="YSO") type= eYSO;
	else if(typeStr=="PULSAR") type= ePULSAR;
	else if(typeStr=="HII") type= eHII;
	else if(typeStr=="PN") type= ePN;	
	else if(typeStr=="SNR") type= eSNR;
	/*
	else if(typeStr=="GALAXY_C1") type= eGALAXY_C1;
	else if(typeStr=="GALAXY_C2") type= eGALAXY_C2;
	else if(typeStr=="GALAXY_C3") type= eGALAXY_C3;
	else if(typeStr=="GALAXY_C4") type= eGALAXY_C4;
	else if(typeStr=="GALAXY_C5") type= eGALAXY_C5;
	else if(typeStr=="GALAXY_C6") type= eGALAXY_C6;
	*/
	else if(typeStr=="GALAXY_C1") type= eGALAXY;
	else if(typeStr=="GALAXY_C2") type= eGALAXY;
	else if(typeStr=="GALAXY_C3") type= eGALAXY;
	else if(typeStr=="GALAXY_C4") type= eGALAXY;
	else if(typeStr=="GALAXY_C5") type= eGALAXY;
	else if(typeStr=="GALAXY_C6") type= eGALAXY;
	else type= eUNKNOWN_OBJECT;

	return type;
}

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
#pragma link C++ enum AstroObjectType+;
#pragma link C++ enum AstroObjectLocation+;
#endif

//====================================
//===          CONSTANTS         =====
//====================================
static const double GausSigma2FWHM= 2.*sqrt(2*log(2.));
static const double SrToDeg2= pow(180/TMath::Pi(),2);
static const double SrToArcmin2= SrToDeg2*60*60;
static const double SrToArcsec2= SrToDeg2*3600*3600;
static const double Arcsec2ToDeg2= 7.71615e-8;
static const double Arcsec2ToArcmin2= 2.77778e-4;
static const double Arcsec2ToSr= 1./SrToArcsec2;

}//close namespace


#endif 
