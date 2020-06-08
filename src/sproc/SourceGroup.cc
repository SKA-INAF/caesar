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
* @file SourceGroup.cc
* @class SourceGroup
* @brief SourceGroup class
*
* Class representing source group
* @author S. Riggi
* @date 26/03/2018
*/

#include <SourceGroup.h>
#include <Source.h>
#include <CodeUtils.h>
#include <GraphicsUtils.h>

#include <Image.h>
#include <Contour.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TPad.h>
#include <TExec.h>
#include <TStyle.h>
#include <TFitResult.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;

ClassImp(Caesar::SourceGroup)

namespace Caesar {


SourceGroup::SourceGroup() 
	: TObject()
{
	//Initialize
	Init();

}//close costructor


SourceGroup::~SourceGroup()
{
	//Delete data
	Clear();

}//close destructor

/*
SourceGroup::SourceGroup(const SourceGroup& sgroup) 
	: TObject(sgroup) 
{
  // Copy constructor
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy constuctor called...");
	#endif
  Init();
  ((SourceGroup&)sgroup).Copy(*this);
}

SourceGroup& SourceGroup::operator=(const SourceGroup& sgroup) { 
	// Operator =
  if (this != &sgroup) ((SourceGroup&)sgroup).Copy(*this);
  return *this;
}

void SourceGroup::Copy(TObject& obj) const 
{
	//Copy object
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy TObject ...");
	#endif
	TObject::Copy((Pixel&)obj);

	// Copy variables
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy source variables...");
	#endif
  ((SourceGroup&)obj).m_owned = m_owned;

	//Delete first a previously existing vector
	#ifdef LOGGING_ENABLED
		INFO_LOG("Delete existing source list ...");
	#endif
	for(size_t i=0;i<(((SourceGroup&)obj).m_sources).size();i++){
		if( (((SourceGroup&)obj).m_sources)[i] ){
			delete (((SourceGroup&)obj).m_sources)[i];
			(((SourceGroup&)obj).m_sources)[i]= 0;
		}
	}
	(((SourceGroup&)obj).m_sources).clear();

	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy source list ...");
	#endif
	Source* source= 0;
	for(unsigned int i=0;i<m_sources.size();i++){
		source= new Source;
		*source= *(m_sources[i]);
		(((SourceGroup&)obj).m_sources).push_back(source);
	}
	#ifdef LOGGING_ENABLED
		INFO_LOG("End copy source group...");
	#endif

}//close Copy()
	
*/


void SourceGroup::Init()
{
	//Init list
	m_sources.clear();	
	m_owned= false;	

}//close Init()

void SourceGroup::Clear()
{
	//Delete collection of source
	if(m_owned){
		CodeUtils::DeletePtrCollection<Source>(m_sources);
	}

}//close Clear()

int SourceGroup::AddSource(Source* aSource,int clone)
{
	if(!aSource) return -1;
	if(clone){
		Source* aNewSource= new Source;
		*aNewSource= *aSource;
		m_sources.push_back(aNewSource);
		m_owned= true;
	}
	else{
		m_sources.push_back(aSource);
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("Added source "<<aSource->GetName()<< " (ref="<<aSource<<") to group (N="<<m_sources.size()<<")");	
	#endif

	return 0;

}//close AddSource()

int SourceGroup::GetSourceNames(std::vector<std::string>& snames)
{
	snames.clear();

	if(m_sources.empty()) return -1;

	for(size_t i=0;i<m_sources.size();i++){
		if(!m_sources[i]){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Matched source no. "<<i+1<<" in group is nullptr, this should not occur!");	
			#endif
			return -1;
		}
		std::string sname= m_sources[i]->GetName();
		snames.push_back(sname);
	}
		
	return 0;
		
}//close GetSourceNames()


int SourceGroup::GetFlux(double& flux,double& fluxErr,bool& summed)
{
	flux= 0;
	fluxErr= 0;
	summed= false;
	if(m_sources.empty()) return -1;
	if(m_sources.size()>1) summed= true;

	double fluxErrSum2= 0;

	for(size_t i=0;i<m_sources.size();i++){
		double fluxDensity= 0;
		double fluxDensityErr= 0;
		int status= m_sources[i]->GetFluxDensity(fluxDensity);
				
		if(status<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot get flux density of source no. "<<i+1<<" in group (no fit info?), using Smax...");	
			#endif
			fluxDensity= m_sources[i]->GetS();
			fluxDensityErr= 0;
		}
		else{
			m_sources[i]->GetFluxDensityErr(fluxDensityErr);
		}
		double beamArea= m_sources[i]->GetBeamFluxIntegral();
		if(beamArea<=0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Cannot get beam area for source no. "<<i+1<<" in group (not computed?), cannot compute flux!");	
			#endif	
			return -1;
		}
		double thisFlux= fluxDensity/beamArea;
		double thisFluxErr= fluxDensityErr/beamArea;	
		flux+= thisFlux;
		fluxErrSum2+= thisFluxErr*thisFluxErr;

	}//end loop sources

	fluxErr= sqrt(fluxErrSum2);

	return 0;
		
}//close GetFlux()


int SourceGroup::GetFrequency(double& freq,double& dfreq)
{
	freq= 0;
	dfreq= 0;
	if(m_sources.empty()) return -1;

	//Get metadata
	ImgMetaData* metadata= m_sources[0]->GetImageMetaData();
	if(!metadata){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get source metadata!");
		#endif
		return -1;
	}

	//Get frequency	
	double Nu= metadata->Freq;
	double dNu= metadata->dFreq;
	std::string FreqUnits= metadata->FreqUnit;
	if(FreqUnits=="Hz"){//convert to GHz
		Nu/= 1.e+9;
		dNu/= 1.e+9;
	}

	freq= Nu;
	dfreq= dNu;
			
	return 0;

}//close GetFrequency()

}//close namespace

