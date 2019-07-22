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
ClassImp(Caesar::SpectralIndexData)
ClassImp(Caesar::ComponentMatchIndex)
ClassImp(Caesar::ComponentMatchIndexGroup)

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
	//CodeUtils::DeletePtrCollection<SourceGroup>(m_matchedSources);
	for(size_t i=0;i<m_matchedSources.size();i++){
		for(size_t j=0;j<m_matchedSources[i].size();j++){
			if(m_matchedSources[i][j]){
				delete m_matchedSources[i][j];
				m_matchedSources[i][j]= 0;
			}
		}	
	}
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

	((SourceMatchData&)obj).m_spectralIndexData = m_spectralIndexData;
 	(((SourceMatchData&)obj).m_componentSpectralIndexData).clear();
	for(size_t i=0;i<m_componentSpectralIndexData.size();i++){
		(((SourceMatchData&)obj).m_componentSpectralIndexData).push_back(SpectralIndexData());
	}

	/*
	((SourceMatchData&)obj).m_hasSpectralIndex = m_hasSpectralIndex;
	((SourceMatchData&)obj).m_isMultiSourceMatchIndex = m_isMultiSourceMatchIndex;
	((SourceMatchData&)obj).m_spectralIndex = m_spectralIndex;
	((SourceMatchData&)obj).m_spectralIndexErr = m_spectralIndexErr;
	((SourceMatchData&)obj).m_isSpectralIndexFit = m_isSpectralIndexFit;
	((SourceMatchData&)obj).m_spectralFitChi2 = m_spectralFitChi2;
	((SourceMatchData&)obj).m_spectralFitNDF = m_spectralFitNDF;
	*/

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
	for(size_t i=0;i<m_sourceComponentSED.size();i++){
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
		for(size_t j=0;j<(((SourceMatchData&)obj).m_matchedSources[i]).size();j++){
			if( (((SourceMatchData&)obj).m_matchedSources)[i][j] ){
				delete (((SourceMatchData&)obj).m_matchedSources)[i][j];
				(((SourceMatchData&)obj).m_matchedSources)[i][j]= 0;
			}
		}
	}
	(((SourceMatchData&)obj).m_matchedSources).clear();

	Source* aSource= 0;
	for(size_t i=0;i<m_matchedSources.size();i++){
		for(size_t j=0;j<m_matchedSources[i].size();j++){
			aSource= new Source;
			*aSource= *(m_matchedSources[i][j]);
			(((SourceMatchData&)obj).m_matchedSources)[i].push_back(aSource);
		}
	}

	/*
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
	*/

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
	m_componentSpectralIndexData.clear();

	//Init source SEDs
	m_sourceSED= nullptr;
	m_sourceComponentSED.clear();

	//Init source group data
	m_nMatches= 0;
	m_componentMatchIndexes.clear();
	for(int i=0;i<m_nCatalogs;i++){
		m_matchedSources.push_back( std::vector<Source*>());
	}
	
	//Init component match index
	if(m_source){
		int nComponents= m_source->GetNFitComponents();
		if(nComponents>0) {
			for(int i=0;i<nComponents;i++){
				//m_componentMatchIndexes.push_back( std::vector<ComponentMatchIndex>() );
				m_componentMatchIndexes.push_back( std::vector<ComponentMatchIndexGroup>() );
				for(int j=0;j<m_nCatalogs;j++){
					m_componentMatchIndexes[i].push_back(ComponentMatchIndexGroup());
				}
			}
			m_componentSpectralIndexData.push_back( SpectralIndexData() );
		}
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
	//m_componentMatchIndexes[componentIndex].push_back( ComponentMatchIndex(catalogIndex,matchedSourceIndex,matchedComponentIndex) );
	m_componentMatchIndexes[componentIndex][catalogIndex].AddIndex(catalogIndex,matchedSourceIndex,matchedComponentIndex);				

	return 0;
		
}//close AddMatchedComponentIndex()



int SourceMatchData::ComputeSourceSEDs()
{
	//Do nothing if no source set
	if(!m_source) return -1;
	if(m_matchedSources.empty()) return 0;

	//Delete existing source SEDs
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Delete existing SED graphs...");
	#endif
	CodeUtils::DeletePtr<TGraphAsymmErrors>(m_sourceSED);
	//m_isMultiSourceMatchIndex= false;
	//m_hasSpectralIndex= false;
	m_spectralIndexData.isMultiSourceMatchIndex= false;
	m_spectralIndexData.hasSpectralIndex= false;
	

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
	double fluxDensity= 0;
	double fluxDensityErr= 0;
	if(m_source->GetFluxDensity(fluxDensity)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get source flux density!");
		#endif
		return -1;
	}
	if(m_source->GetFluxDensityErr(fluxDensityErr)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get source flux density error!");
		#endif
		return -1;
	}
	double beamArea= m_source->GetBeamFluxIntegral();
	if(beamArea<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get beam area (cannot compute flux)!");
		#endif
		return -1;
	}	

	double flux= fluxDensity/beamArea;
	double fluxErr= fluxDensityErr/beamArea;
	double lgFlux= log10(flux);
	double lgFluxErr= log10(TMath::E())*fluxErr/flux;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<": Nu="<<Nu<<", dNu="<<dNu<<", FreqUnits="<<FreqUnits<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
	#endif
	
	lgNu_list.push_back(lgNu);
	lgFlux_list.push_back(lgFlux);
	lgFluxErr_list.push_back(lgFluxErr);

	//Check if we have matched sources
	if(!HasSourceMatch()) {
		m_sourceSED= new TGraphAsymmErrors(1);
		m_sourceSED->SetPoint(0,lgNu,lgFlux);
		m_sourceSED->SetPointError(0,0,0,lgFluxErr,lgFluxErr);
		return 0;//nothing to be done without matched sources
	}
	

	//## Get matched source info
	for(size_t i=0;i<m_matchedSources.size();i++)
	{
		//Skip if no matches
		if(m_matchedSources[i].empty()) continue;

		//Get source names
		std::vector<std::string> snames= GetMatchedSourceNames(i);
		std::stringstream ss;
		ss<<"{";
		for(size_t j=0;j<snames.size()-1;j++) ss<<snames[j]<<",";
		ss<<snames[snames.size()-1]<<"}"; 
		

		//Get frequency
		if(GetMatchedSourceFrequency(i,Nu,dNu)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to get frequency for matched source no. "<<i+1<<"!");
			#endif
			return -1;
		}
		lgNu= log10(Nu);
		dlgNu= log10(TMath::E())*dNu/Nu;

		//Get flux
		bool fluxSummed= false;
		if(GetMatchedSourceFlux(i,flux,fluxErr,fluxSummed)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to get flux for matched source no. "<<i+1<<"!");
			#endif
			return -1;
		}
		lgFlux= log10(flux);
		lgFluxErr= log10(TMath::E())*fluxErr/flux;
		//if(fluxSummed) m_isMultiSourceMatchIndex= true;
		if(fluxSummed) m_spectralIndexData.isMultiSourceMatchIndex= true;

		#ifdef LOGGING_ENABLED
			INFO_LOG("Matched source group no. "<<i+1<<": snames="<<ss.str()<<", Nu="<<Nu<<", dNu="<<dNu<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
		#endif

		lgNu_list.push_back(lgNu);
		lgFlux_list.push_back(lgFlux);
		lgFluxErr_list.push_back(lgFluxErr);

	}//end loop matched sources
	
	/*
	for(size_t i=0;i<m_matchedSources.size();i++)
	{
		
		//Get source names
		//std::vector<std::string> snames= m_matchedSources[i]->GetSourceNames();
		//std::stringstream ss;
		//ss<<"{";
		//for(size_t j=0;j<snames.size()-1;j++) ss<<snames[j]<<",";
		//ss<<snames[snames.size()-1]<<"}"; 
		

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
	*/

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
		//m_isSpectralIndexFit= false;
		m_spectralIndexData.isSpectralIndexFit= false;
		double x1= 0; 
		double x2= 0;	
		double y1= 0;
		double y2= 0;
		m_sourceSED->GetPoint(0,x1,y1);
		m_sourceSED->GetPoint(1,x2,y2);
		//m_spectralIndex= (y2-y1)/(x2-x1);
		//m_spectralIndexErr= 0;
		//m_spectralFitChi2= 0;
		//m_spectralFitNDF= 0;
		//m_hasSpectralIndex= true;
		m_spectralIndexData.spectralIndex= (y2-y1)/(x2-x1);
		m_spectralIndexData.spectralIndexErr= 0;
		m_spectralIndexData.spectralFitChi2= 0;
		m_spectralIndexData.spectralFitNDF= 0;
		m_spectralIndexData.hasSpectralIndex= true;
	}
	else{
		//m_isSpectralIndexFit= true;
		m_spectralIndexData.isSpectralIndexFit= true;
		TFitResultPtr fitRes= m_sourceSED->Fit("pol1","S");
		//m_spectralFitChi2= fitRes->Chi2();
		//m_spectralFitNDF= fitRes->Ndf();
		//m_spectralIndex= fitRes->Value(1);
		//m_spectralIndexErr= fitRes->ParError(1);

		m_spectralIndexData.spectralFitChi2= fitRes->Chi2();
		m_spectralIndexData.spectralFitNDF= fitRes->Ndf();
		m_spectralIndexData.spectralIndex= fitRes->Value(1);
		m_spectralIndexData.spectralIndexErr= fitRes->ParError(1);
		int fitStatus= fitRes;
		//if(fitStatus==0) m_hasSpectralIndex= true;
		if(fitStatus==0) m_spectralIndexData.hasSpectralIndex= true;
	}

	#ifdef LOGGING_ENABLED
		//INFO_LOG("Source "<<m_source->GetName()<<": hasSpectralIndex?"<<m_hasSpectralIndex<<", gamma="<<m_spectralIndex<<" +- "<<m_spectralIndexErr<<", chi2/ndf="<<m_spectralFitChi2<<"/"<<m_spectralFitNDF);
		INFO_LOG("Source "<<m_source->GetName()<<": hasSpectralIndex?"<<m_spectralIndexData.hasSpectralIndex<<", gamma="<<m_spectralIndexData.spectralIndex<<" +- "<<m_spectralIndexData.spectralIndexErr<<", chi2/ndf="<<m_spectralIndexData.spectralFitChi2<<"/"<<m_spectralIndexData.spectralFitNDF);
	#endif
	

	return 0;

}//close ComputeSourceSEDs()


int SourceMatchData::ComputeSourceComponentSEDs()
{
	//Do nothing if no source set
	if(!m_source) return -1;
	if(m_matchedSources.empty()) return 0;
	if(m_componentMatchIndexes.empty()) return 0;

	//Get source pars
	if(!m_source->HasFitInfo()) return 0;//no components, nothing to be done

	std::string sourceName= m_source->GetName();
	double beamArea= m_source->GetBeamFluxIntegral();
	if(beamArea<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cannot get beam area for source (not computed?), cannot compute flux!");	
		#endif	
		return -1;
	}

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

	
	//Get source fit component pars
	SourceFitPars fitPars= m_source->GetFitPars();

	
	//Delete existing source SEDs
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Delete existing SED graphs...");
	#endif
	CodeUtils::DeletePtrCollection<TGraphAsymmErrors>(m_sourceComponentSED);
	
	
	//Loop over source components
	TGraphAsymmErrors* sed= 0;

	for(size_t i=0;i<m_componentMatchIndexes.size();i++)
	{
		std::vector<double> lgNu_list;
		std::vector<double> lgFlux_list;
		std::vector<double> lgFluxErr_list;

		//Skip source component if not selected
		bool isSelected= fitPars.IsSelectedComponent(i);
		if(!isSelected) continue;

		//Get source component pars
		std::string sname= sourceName + std::string(Form("_fitcomp%d",i+1));
		double fluxDensity= fitPars.GetComponentFluxDensity(i);
		double fluxDensityErr= fitPars.GetComponentFluxDensityErr(i);
		double flux= fluxDensity/beamArea;
		double fluxErr= fluxDensityErr/beamArea;
		double lgFlux= log10(flux);
		double lgFluxErr= log10(TMath::E())*fluxErr/flux;
		
		lgNu_list.push_back(lgNu);
		lgFlux_list.push_back(lgFlux);
		lgFluxErr_list.push_back(lgFluxErr);
		bool isMatchComponentFluxSummed= false;	

		//Skip if no component match	
		if(!HasSourceComponentMatch(i)){
			sed= new TGraphAsymmErrors(1);
			sed->SetPoint(0,lgNu,lgFlux);
			sed->SetPointError(0,0,0,lgFluxErr,lgFluxErr);
			m_sourceComponentSED.push_back(sed);
			continue;
		}

		//Loop over catalogs
		for(size_t j=0;j<m_componentMatchIndexes[i].size();j++)
		{
			//Loop over matched components per catalog
			std::vector<ComponentMatchIndex> index_list= m_componentMatchIndexes[i][j].GetIndexes();
			if(index_list.empty()){
				continue;
			}

			double flux_matched= 0;
			double fluxErr_matched= 0;	
			double lgNu_matched= 0;
			double dlgNu_matched= 0;
			double fluxErrSum2_matched= 0;
			int nMatches= 0;
		
			for(size_t k=0;k<index_list.size();k++)
			{
				int catalogIndex= index_list[k].catalogIndex;
				int matchedSourceIndex= index_list[k].sourceGroupIndex;
				int matchedComponentIndex= index_list[k].componentIndex;
				Source* matchedSource= m_matchedSources[catalogIndex][matchedSourceIndex];

				//Skip if no fit pars
				if(!matchedSource->HasFitInfo()) continue;	
				SourceFitPars fitPars_matched= matchedSource->GetFitPars();

				//Get beam area	
				double beamArea_matched= matchedSource->GetBeamFluxIntegral();
				if(beamArea_matched<=0){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Cannot get beam area for matched source (not computed?), cannot compute flux!");	
					#endif	
					return -1;
				}		
			
				//Get metadata
				ImgMetaData* metadata= matchedSource->GetImageMetaData();
				if(!metadata){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to get matched source metadata!");
					#endif
					return -1;
				}

				//Get component frequency	
				Nu= metadata->Freq;
				dNu= metadata->dFreq;
				std::string FreqUnits= metadata->FreqUnit;
				if(FreqUnits=="Hz"){//convert to GHz
					Nu/= 1.e+9;
					dNu/= 1.e+9;
				}
				lgNu_matched= log10(Nu);
				dlgNu_matched= log10(TMath::E())*dNu/Nu;

				//Get component flux density
				double fluxDensity_matched= fitPars_matched.GetComponentFluxDensity(matchedComponentIndex);
				double fluxDensityErr_matched= fitPars_matched.GetComponentFluxDensityErr(matchedComponentIndex);
				double thisFlux= fluxDensity_matched/beamArea_matched;
				double thisFluxErr= fluxDensityErr_matched/beamArea_matched;
				flux_matched+= thisFlux;
				fluxErrSum2_matched+= thisFluxErr*thisFluxErr;
				nMatches++;

			}//end loop matched components per catalog


			fluxErr_matched= sqrt(fluxErrSum2_matched);
			bool fluxSummed= (nMatches>1);
			double lgFlux_matched= log10(flux_matched);
			double lgFluxErr_matched= log10(TMath::E())*fluxErr_matched/flux_matched;
		
			lgNu_list.push_back(lgNu_matched);
			lgFlux_list.push_back(lgFlux_matched);
			lgFluxErr_list.push_back(lgFluxErr_matched);
			if(fluxSummed) isMatchComponentFluxSummed= true;

		}//end loop catalogs


		//# Sort data by ascending frequencies
		std::vector<double> lgNu_list_sorted;
		std::vector<double> lgFlux_list_sorted;
		std::vector<double> lgFluxErr_list_sorted;
		std::vector<size_t> sorting_indexes;
		CodeUtils::sort(lgNu_list,lgNu_list_sorted,sorting_indexes);
		for(size_t k=0;k<sorting_indexes.size();k++)
		{
			int index= sorting_indexes[k];
			lgFlux_list_sorted.push_back(lgFlux_list[index]);
			lgFluxErr_list_sorted.push_back(lgFluxErr_list[index]);
		}

		//# Fill source SED
		int N= static_cast<int>(lgNu_list_sorted.size());
		sed= new TGraphAsymmErrors(N);
		
		for(int k=0;k<N;k++){
			sed->SetPoint(k,lgNu_list_sorted[k],lgFlux_list_sorted[k]);
			sed->SetPointError(k,0,0,lgFluxErr_list_sorted[k],lgFluxErr_list_sorted[k]);
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source SED: P"<<k+1<<"("<<lgNu_list_sorted[k]<<","<<lgFlux_list_sorted[k]<<")");
			#endif
		}
		m_sourceComponentSED.push_back(sed);
	

		//# Find spectral index
		SpectralIndexData spectralIndexData;
		spectralIndexData.isMultiSourceMatchIndex= isMatchComponentFluxSummed;
		if(N==2){
			spectralIndexData.isSpectralIndexFit= false;
			double x1= 0; 
			double x2= 0;	
			double y1= 0;
			double y2= 0;
			sed->GetPoint(0,x1,y1);
			sed->GetPoint(1,x2,y2);
			spectralIndexData.spectralIndex= (y2-y1)/(x2-x1);
			spectralIndexData.spectralIndexErr= 0;
			spectralIndexData.spectralFitChi2= 0;
			spectralIndexData.spectralFitNDF= 0;
			spectralIndexData.hasSpectralIndex= true;
		}
		else{
			spectralIndexData.isSpectralIndexFit= true;
			TFitResultPtr fitRes= sed->Fit("pol1","S");
			spectralIndexData.spectralFitChi2= fitRes->Chi2();
			spectralIndexData.spectralFitNDF= fitRes->Ndf();
			spectralIndexData.spectralIndex= fitRes->Value(1);
			spectralIndexData.spectralIndexErr= fitRes->ParError(1);
			int fitStatus= fitRes;
			if(fitStatus==0) spectralIndexData.hasSpectralIndex= true;
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("Source "<<m_source->GetName()<<", component no. "<<i+1<<": hasSpectralIndex?"<<spectralIndexData.hasSpectralIndex<<", gamma="<<spectralIndexData.spectralIndex<<" +- "<<spectralIndexData.spectralIndexErr<<", chi2/ndf="<<spectralIndexData.spectralFitChi2<<"/"<<spectralIndexData.spectralFitNDF);
		#endif

		m_componentSpectralIndexData[i]= spectralIndexData;
		

	}//end loop components
	

	return 0;

}//close ComputeSourceComponentSEDs()


int SourceMatchData::DrawSED()
{
	//Do nothing if no source is present
	if(!m_source || !m_sourceSED){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No source and/or sed present, nothing to be drawn.");
		#endif	
		return 0;
	}	
	if(m_sourceSED->GetN()<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No points stored in source SED, nothing to be drawn.");
		#endif	
		return 0;
	}

	//Find SED min/max
	double lgNu_min= 1.e+99;
	double lgNu_max= -1.e+99;
	double lgFlux_min= 1.e+99;
	double lgFlux_max= -1.e+99;
	for(int i=0;i<m_sourceSED->GetN();i++){
		double x;
		double y;
		m_sourceSED->GetPoint(i,x,y);
		if(x<lgNu_min) lgNu_min= x;
		if(x>lgNu_max) lgNu_max= x;
		if(y<lgFlux_min) lgFlux_min= y;
		if(y>lgFlux_max) lgFlux_max= y;
	}

	//Get metadata info
	std::string sname= m_source->GetName();
	std::string units= "";
	TString freq= "";
	ImgMetaData* metadata= m_source->GetImageMetaData();
	if(metadata){
		double Nu= metadata->Freq;
		std::string FreqUnits= metadata->FreqUnit;
		if(FreqUnits=="Hz"){//convert to GHz
			Nu/= 1.e+9;
			FreqUnits= "GHz";
		}
		freq= Form("%1.2f %s",Nu,(FreqUnits).c_str());
	}

	//Get WCS position
	int wcsType= eJ2000;
	WCS* wcs= m_source->GetWCS(wcsType);
	double X0= m_source->X0;
	double Y0= m_source->Y0;
	if(wcs){
		m_source->GetWCSPos(X0,Y0,wcs,wcsType);
		WCSUtils::DeleteWCS(&wcs);
	}
	

	//Get flux
	double fluxDensity= 0;
	int status= m_source->GetFluxDensity(fluxDensity);
	double beamArea= m_source->GetBeamFluxIntegral();
	double flux= fluxDensity;
	TString fluxText= "";
	std::string fluxUnits= "mJy/beam";
	if(beamArea>0) {
		flux/= fluxDensity;
		fluxUnits= "mJy";
	}
	fluxText= Form("%1.2f %s",flux,(fluxUnits).c_str());

	//Create canvas
	TString canvasName= Form("SEDPlot_%s",sname.c_str());
	TCanvas* Plot= new TCanvas(canvasName,canvasName,800,800);
	Plot->cd();

	double step= 0.1;
	double x_min= lgNu_min - fabs(step*(lgNu_min-lgNu_max));
	double x_max= lgNu_min + fabs(step*(lgNu_min-lgNu_max));
	double y_min= lgFlux_min - fabs(step*(lgFlux_min-lgFlux_max));
	double y_max= lgFlux_min + fabs(step*(lgFlux_min-lgFlux_max));

	TH2D* PlotBkg= new TH2D("PlotBkg","",100,x_min,x_max,100,y_min,y_max);
	PlotBkg->GetXaxis()->SetTitle("log_{10}(#nu/Hz)");
	PlotBkg->GetYaxis()->SetTitle("log_{10}(S/Jy)");	
	PlotBkg->SetStats(0);
	PlotBkg->Draw();

	//Draw source SED
	m_sourceSED->SetMarkerStyle(8);
	m_sourceSED->SetMarkerColor(kBlack);
	m_sourceSED->SetLineColor(kBlack);
	m_sourceSED->Draw("EPZL same");
	
	//Draw main source info
	TPaveText* sourceInfoText = new TPaveText(0.4,0.15,0.8,0.3,"NDC");
	sourceInfoText->AddText(Form("Name: %s, Freq: %s",sname.c_str(),freq.Data()));
	sourceInfoText->AddText(Form("C(%1.2f,%1.2f), S(mJy): %s",X0,Y0,fluxText.Data()));
	sourceInfoText->SetTextAlign(12);
	sourceInfoText->SetTextSize(0.02);
	sourceInfoText->SetTextFont(52);
	sourceInfoText->SetFillColor(0);
	sourceInfoText->SetBorderSize(1);	
	sourceInfoText->Draw("same");

	//Draw source component SEDs
	std::vector<int> componentColors {kRed,kGreen+1,kBlue,kMagenta,kYellow+1,kOrange,kGray};
		
	if(!m_sourceComponentSED.empty()){
		for(size_t i=0;i<m_sourceComponentSED.size();i++){	
			int color= kBlack;
			if(m_sourceComponentSED.size()<componentColors.size()){
				color= componentColors[i];
			}
			m_sourceComponentSED[i]->SetMarkerStyle(21);
			m_sourceComponentSED[i]->SetMarkerColor(color);
			m_sourceComponentSED[i]->SetLineColor(color);
			m_sourceComponentSED[i]->Draw("EPZL same");
		}
	}

	return 0;

}//close DrawSED()

}//close namespace 


