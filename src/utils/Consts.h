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
	eUnknown=0,
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

enum DS9RegionFormat {
	eEllipseRegion= 1,
	ePolygonRegion= 2
};

#ifdef __MAKECINT__
#pragma link C++ enum WCSType+;
#pragma link C++ enum ImgFilters+;
#pragma link C++ enum ImgType+;
#pragma link C++ enum ImgSmoothFilter+;
#pragma link C++ enum SegmAlgo+;
#pragma link C++ enum BkgEstimator+;
#pragma link C++ enum BkgMethod+;
#pragma link C++ enum ColorPaletteStyle+;
#pragma link C++ enum FileType+;
#pragma link C++ enum SLICEdgeModel+;
#endif

//====================================
//===          CONSTANTS         =====
//====================================
static const double GausSigma2FWHM= 2.*sqrt(2*log(2.));






}//close namespace


#endif 
