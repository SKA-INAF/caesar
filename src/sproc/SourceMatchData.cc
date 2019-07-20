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
* @file SourceMatchData.cc
* @class SourceMatchData
* @brief SourceMatchData class
*
* Class representing source match data
* @author S. Riggi
* @date 26/03/2018
*/

#include <SourceMatchData.h>
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

ClassImp(Caesar::SourceMatchData)
ClassImp(Caesar::SourceGroup)
ClassImp(Caesar::ComponentMatchIndex)

namespace Caesar {


SourceMatchData::SourceMatchData() 
	: TObject()
{
	//Initialize
	m_source= nullptr;
	m_nCatalogs= 0;
	Init();

}//close costructor

SourceMatchData::SourceMatchData(Source* aSource,int nCatalogs) 
	: TObject(), m_source(aSource), m_nCatalogs(nCatalogs)
{
	//Initialize 
	Init();

}//close costructor


SourceMatchData::~SourceMatchData()
{
	//Delete source 
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source ...");
	#endif
	CodeUtils::DeletePtr<Source>(m_source);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif

	//Delete source collection
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source match collection added ...");
	#endif
	CodeUtils::DeletePtrCollection<SourceGroup>(m_matchedSources);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif

	//Delete source SED graph
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source SED graph ...");
	#endif
	CodeUtils::DeletePtr<TGraphAsymmErrors>(m_sourceSED);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif

	//Delete source component SED collection
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source component SED graph added ...");
	#endif
	CodeUtils::DeletePtrCollection<TGraphAsymmErrors>(m_sourceComponentSED);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif

	
}//close destructor


SourceMatchData::SourceMatchData(const SourceMatchData& sourceMatchData) 
{
  // Copy constructor
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy constuctor called...");
	#endif
  Init();
  ((SourceMatchData&)sourceMatchData).Copy(*this);
}

void SourceMatchData::Copy(TObject &obj) const 
{
	//Copy vars
	((SourceMatchData&)obj).m_nCatalogs = m_nCatalogs;
	((SourceMatchData&)obj).m_nMatches = m_nMatches;
	((SourceMatchData&)obj).m_hasSpectralIndex = m_hasSpectralIndex;
	((SourceMatchData&)obj).m_isMultiSourceMatchIndex = m_isMultiSourceMatchIndex;
	((SourceMatchData&)obj).m_spectralIndex = m_spectralIndex;
	((SourceMatchData&)obj).m_spectralIndexErr = m_spectralIndexErr;
	((SourceMatchData&)obj).m_isSpectralIndexFit = m_isSpectralIndexFit;
	((SourceMatchData&)obj).m_spectralFitChi2 = m_spectralFitChi2;
	((SourceMatchData&)obj).m_spectralFitNDF = m_spectralFitNDF;
	
	//Copy source SED graph
	if(((SourceMatchData&)obj).m_sourceSED){
		delete ((SourceMatchData&)obj).m_sourceSED;
		((SourceMatchData&)obj).m_sourceSED= 0;
	}
	if(m_sourceSED){
		((SourceMatchData&)obj).m_sourceSED= new TGraphAsymmErrors;
		*((SourceMatchData&)obj).m_sourceSED = *m_sourceSED;
	}

	//Copy source SED component collection
	//Delete first any existing collection
	for(size_t i=0;i<(((SourceMatchData&)obj).m_sourceComponentSED).size();i++){
		if( (((SourceMatchData&)obj).m_sourceComponentSED)[i] ){
			delete (((SourceMatchData&)obj).m_sourceComponentSED)[i];
			(((SourceMatchData&)obj).m_sourceComponentSED)[i]= 0;
		}
	}
	(((SourceMatchData&)obj).m_sourceComponentSED).clear();

	TGraphAsymmErrors* aGraph= 0;
	for(size_t i=0;i<m_matchedSources.size();i++){
		aGraph= new TGraphAsymmErrors;
		*aGraph= *(m_sourceComponentSED[i]);
		(((SourceMatchData&)obj).m_sourceComponentSED).push_back(aGraph);
	}

	//Copy source
	if(((SourceMatchData&)obj).m_source){
		delete ((SourceMatchData&)obj).m_source;
		((SourceMatchData&)obj).m_source= 0;
	}
	if(m_source){
		((SourceMatchData&)obj).m_source= new Source;
		*((SourceMatchData&)obj).m_source = *m_source;
	}

	//Copy source collection
	//Delete first any existing collection
	for(size_t i=0;i<(((SourceMatchData&)obj).m_matchedSources).size();i++){
		if( (((SourceMatchData&)obj).m_matchedSources)[i] ){
			delete (((SourceMatchData&)obj).m_matchedSources)[i];
			(((SourceMatchData&)obj).m_matchedSources)[i]= 0;
		}
	}
	(((SourceMatchData&)obj).m_matchedSources).clear();

	SourceGroup* aSourceGroup= 0;
	for(size_t i=0;i<m_matchedSources.size();i++){
		aSourceGroup= new SourceGroup;
		*aSourceGroup= *(m_matchedSources[i]);
		(((SourceMatchData&)obj).m_matchedSources).push_back(aSourceGroup);
	}

}//close Copy()

SourceMatchData& SourceMatchData::operator=(const SourceMatchData& sourceMatchData) 
{ 
	// Operator =
  if (this != &sourceMatchData) ((SourceMatchData&)sourceMatchData).Copy(*this);
  return *this;
}


void SourceMatchData::Init()
{
	//Init spectral index
	m_hasSpectralIndex= false;
	m_spectralIndex= -999;
	m_spectralIndexErr= 0;
	m_isMultiSourceMatchIndex= false;
	m_spectralFitChi2= -999;
	m_spectralFitNDF= 0;
	m_isSpectralIndexFit= false;

	//Init source SEDs
	m_sourceSED= nullptr;
	m_sourceComponentSED.clear();

	//Init source group data
	m_nMatches= 0;
	m_componentMatchIndexes.clear();

	SourceGroup* aSourceGroup= 0;
	for(int i=0;i<m_nCatalogs;i++){
		aSourceGroup= new SourceGroup;
		m_matchedSources.push_back(aSourceGroup);
	}

	//Init component match index
	if(m_source){
		int nComponents= m_source->GetNFitComponents();
		//if(nComponents>0) m_componentMatchIndexes.push_back( std::vector<std::pair<size_t,size_t>>(nComponents) );
		if(nComponents>0) m_componentMatchIndexes.push_back( std::vector<ComponentMatchIndex>(nComponents) );
	}

}//close Init()



int SourceMatchData::AddMatchedComponentIndex(int componentIndex,int catalogIndex,int matchedSourceIndex,int matchedComponentIndex)
{
	//Check if cube component id was allocated
	int nComponents= static_cast<int>(m_componentMatchIndexes.size());
	if(nComponents<=0 || nComponents<componentIndex+1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Component with index "<<componentIndex<<" was not allocated!");
		#endif
		return -1;
	}
	
	//Fill component indexes
	//m_componentMatchIndexes[componentIndex].push_back(std::make_pair(matchedSourceIndex,matchedComponentIndex));
	m_componentMatchIndexes[componentIndex].push_back( ComponentMatchIndex(catalogIndex,matchedSourceIndex,matchedComponentIndex) );
				
	return 0;
		
}//close AddMatchedComponentIndex()



int SourceMatchData::ComputeSourceSEDs()
{
	//Do nothing if no source set
	if(!m_source) return -1;
	if(m_matchedSources.empty()) return 0;//nothing to be done without matched sources
	
	//Delete existing source SEDs
	#ifdef LOGGING_ENABLED
		INFO_LOG("Delete existing SED graphs...");
	#endif
	CodeUtils::DeletePtr<TGraphAsymmErrors>(m_sourceSED);
	CodeUtils::DeletePtrCollection<TGraphAsymmErrors>(m_sourceComponentSED);
	m_isMultiSourceMatchIndex= false;
	m_hasSpectralIndex= false;

	//## Get source spectral info
	std::vector<double> lgNu_list;
	std::vector<double> lgFlux_list;
	std::vector<double> lgFluxErr_list;

	//Get source metadata
	ImgMetaData* metadata= m_source->GetImageMetaData();
	if(!metadata){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get source metadata!");
		#endif
		return -1;
	}
	
	//Get source frequency
	double Nu= metadata->Freq;
	double dNu= metadata->dFreq;
	std::string FreqUnits= metadata->FreqUnit;
	if(FreqUnits=="Hz"){//convert to GHz
		Nu/= 1.e+9;
		dNu/= 1.e+9;
	}
	double lgNu= log10(Nu);
	double dlgNu= log10(TMath::E())*dNu/Nu;	

	
	
	//Get source flux
	double flux= 0;
	double fluxErr= 0;
	if(m_source->GetFluxDensity(flux)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get source flux density!");
		#endif
		return -1;
	}
	if(m_source->GetFluxDensityErr(fluxErr)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get source flux density error!");
		#endif
		return -1;
	}
	double lgFlux= log10(flux);
	double lgFluxErr= log10(TMath::E())*fluxErr/flux;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<": Nu="<<Nu<<", dNu="<<dNu<<", FreqUnits="<<FreqUnits<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
	#endif
	
	lgNu_list.push_back(lgNu);
	lgFlux_list.push_back(lgFlux);
	lgFluxErr_list.push_back(lgFluxErr);


	//## Get matched source info
	for(size_t i=0;i<m_matchedSources.size();i++)
	{
		/*
		//Get source names
		std::vector<std::string> snames= m_matchedSources[i]->GetSourceNames();
		std::stringstream ss;
		ss<<"{";
		for(size_t j=0;j<snames.size()-1;j++) ss<<snames[j]<<",";
		ss<<snames[snames.size()-1]<<"}"; 
		*/

		//Get frequency
		if(m_matchedSources[i]->GetFrequency(Nu,dNu)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to get frequency for matched source no. "<<i+1<<"!");
			#endif
			return -1;
		}
		lgNu= log10(Nu);
		dlgNu= log10(TMath::E())*dNu/Nu;	
		
		//Get flux
		bool fluxSummed= false;
		if(m_matchedSources[i]->GetFlux(flux,fluxErr,fluxSummed)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to get flux for matched source no. "<<i+1<<"!");
			#endif
			return -1;
		}
		lgFlux= log10(flux);
		lgFluxErr= log10(TMath::E())*fluxErr/flux;
		if(fluxSummed) m_isMultiSourceMatchIndex= true;

		#ifdef LOGGING_ENABLED
			//INFO_LOG("Matched source group no. "<<i+1<<": snames="<<ss.str()<<", Nu="<<Nu<<", dNu="<<dNu<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
			INFO_LOG("Matched source group no. "<<i+1<<": Nu="<<Nu<<", dNu="<<dNu<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
		#endif

		lgNu_list.push_back(lgNu);
		lgFlux_list.push_back(lgFlux);
		lgFluxErr_list.push_back(lgFluxErr);

	}//end loop matched sources

	//# Sort data by ascending frequencies
	std::vector<double> lgNu_list_sorted;
	std::vector<double> lgFlux_list_sorted;
	std::vector<double> lgFluxErr_list_sorted;
	std::vector<size_t> sorting_indexes;
	CodeUtils::sort(lgNu_list,lgNu_list_sorted,sorting_indexes);
	for(size_t i=0;i<sorting_indexes.size();i++)
	{
		int index= sorting_indexes[i];
		lgFlux_list_sorted.push_back(lgFlux_list[index]);
		lgFluxErr_list_sorted.push_back(lgFluxErr_list[index]);
	}

	//# Fill source SED
	int N= static_cast<int>(lgNu_list_sorted.size());
	m_sourceSED= new TGraphAsymmErrors(N);
	
	for(int i=0;i<N;i++){
		m_sourceSED->SetPoint(i,lgNu_list_sorted[i],lgFlux_list_sorted[i]);
		m_sourceSED->SetPointError(i,0,0,lgFluxErr_list_sorted[i],lgFluxErr_list_sorted[i]);
		#ifdef LOGGING_ENABLED
			INFO_LOG("Source SED: P"<<i+1<<"("<<lgNu_list_sorted[i]<<","<<lgFlux_list_sorted[i]<<")");
		#endif
	}

	//# Find spectral index
	if(N==2){
		m_isSpectralIndexFit= false;
		double x1= 0; 
		double x2= 0;	
		double y1= 0;
		double y2= 0;
		m_sourceSED->GetPoint(0,x1,y1);
		m_sourceSED->GetPoint(1,x2,y2);
		m_spectralIndex= (y2-y1)/(x2-x1);
		m_spectralIndexErr= 0;
		m_spectralFitChi2= 0;
		m_spectralFitNDF= 0;
		m_hasSpectralIndex= true;
	}
	else{
		m_isSpectralIndexFit= true;
		TFitResultPtr fitRes= m_sourceSED->Fit("pol1","S");
		m_spectralFitChi2= fitRes->Chi2();
		m_spectralFitNDF= fitRes->Ndf();
		m_spectralIndex= fitRes->Value(1);
		m_spectralIndexErr= fitRes->ParError(1);
		int fitStatus= fitRes;
		if(fitStatus==0) m_hasSpectralIndex= true;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<": hasSpectralIndex?"<<m_hasSpectralIndex<<", gamma="<<m_spectralIndex<<" +- "<<m_spectralIndexErr<<", chi2/ndf="<<m_spectralFitChi2<<"/"<<m_spectralFitNDF);
	#endif
	

	return 0;

}//close ComputeSourceSEDs()

}//close namespace 


