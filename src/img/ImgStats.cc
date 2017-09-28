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
* @file ImgStats.cc
* @class ImgStats
* @brief ImgStats
*
* Image stats class
* @author S. Riggi
* @date 20/01/2015
*/

#include <ImgStats.h>
#include <Logger.h>


ClassImp(Caesar::ImgStats)

namespace Caesar {

ImgStats::ImgStats()
{
	Reset();
}

ImgStats::~ImgStats()
{

}//close destructor

void ImgStats::Reset(){
			
	n= 0;
	min= 0;
	max= 0;
	mean= 0;
	meanErr= 0;	
	rms= 0;
	rmsErr= 0;
	skewness= 0;
	skewnessErr= 0;
	kurtosis= 0;
	kurtosisErr= 0;
	median= 0;
	medianRMS= 0;
	bwLocation= 0;
	bwScale= 0;
	bwLocationIter= 0;
	bwScaleIter= 0;
	clippedMedian= 0;
	clippedRMS= 0;

}//close Reset()


void ImgStats::Log(std::string level)
{
	//Log stats
	LOG(level,GetPrintable());

}

void ImgStats::Print()
{
			
	cout<<"*** IMG STATS ***"<<endl;
	cout<<"N="<<n<<" min/max="<<min<<"/"<<max<<endl;
	cout<<"Mean: "<<mean<<" +- "<<meanErr<<endl;
	cout<<"RMS: "<<rms<<" +- "<<rmsErr<<endl;
	cout<<"Skewness: "<<skewness<<" +- "<<skewnessErr<<endl;
	cout<<"Kurtosis: "<<kurtosis<<" +- "<<kurtosisErr<<endl;
  cout<<"Median: "<<median<<", MAD: "<<medianRMS<<endl;
	cout<<"BiWeight Location: "<<bwLocation<<", Scale: "<<bwScale<<endl;		
	cout<<"Clipped Median: "<<clippedMedian<<" MAD: "<<clippedRMS<<endl;
	cout<<"*****************"<<endl;

}//close Print()


std::string ImgStats::GetPrintable()
{
			
	std::stringstream ss;
	ss<<"IMG STATS: ";
	ss<<"N="<<n<<" min/max="<<min<<"/"<<max<<", ";
	ss<<"Mean: "<<mean<<" +- "<<meanErr<<", ";
	ss<<"RMS: "<<rms<<" +- "<<rmsErr<<", ";
	ss<<"Skewness: "<<skewness<<" +- "<<skewnessErr<<", ";
	ss<<"Median: "<<median<<", MAD: "<<medianRMS<<", ";
	ss<<"BiWeight Location: "<<bwLocation<<", Scale: "<<bwScale<<", ";
	ss<<"Clipped Median: "<<clippedMedian<<" MAD: "<<clippedRMS;
	return ss.str();

}//close GetPrintable()

}//close namespace
