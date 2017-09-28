
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
* @file ImgStats.h
* @class ImgStats
* @brief ImgStats
*
* Image stats class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _IMG_STATS_h
#define _IMG_STATS_h 1

#include <TObject.h>
#include <string>

namespace Caesar {


class ImgStats : public TObject {

	public:

		/**
		* \brief Constructor
		*/
		ImgStats();
		/**
		* \brief Destructor
		*/
		virtual ~ImgStats();

	public:
		/**
		* \brief Reset stats params
		*/
		void Reset();

		/**
		* \brief Log stats
		*/
		void Log(std::string level="INFO");

		/**
		* \brief Print stats to stdout
		*/
		void Print();

		/**
		* \brief Get a printable string
		*/
		std::string GetPrintable();

	public:		
	
		int n;
		double min;
		double max;
		double mean;
		double meanErr;
		double rms;
		double rmsErr;
  	double skewness;	
		double skewnessErr;			
		double kurtosis;	
		double kurtosisErr;	
		double median;
		double medianRMS;
		double bwLocation;
		double bwScale;
		double bwLocationIter;
		double bwScaleIter;
		double clippedMedian;
		double clippedRMS;
		
	ClassDef(ImgStats,2)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class ImgStats+;
#endif

}//close namespace

#endif
