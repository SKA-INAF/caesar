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
#include <SourceGroup.h>
#include <SpectralIndexData.h>

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
#include <TF1.h>
#include <TSystem.h>
#include <TPad.h>
#include <TExec.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TLegend.h>

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
ClassImp(Caesar::ComponentMatchIndex)
ClassImp(Caesar::ComponentMatchIndexGroup)
ClassImp(Caesar::SourceMatchPars)
ClassImp(Caesar::SourceComponentMatchPars)
ClassImp(Caesar::SourceMatchParsGroup)
ClassImp(Caesar::SourceComponentMatchParsGroup)


namespace Caesar {


SourceMatchData::SourceMatchData() 
	: TObject()
{
	//Initialize
	m_source= nullptr;
	m_ownedSource= false;
	m_ownedMatchedSources= false;
	m_nCatalogs= 0;
	Init();

}//close costructor

SourceMatchData::SourceMatchData(Source* aSource,int nCatalogs,bool cloneSource) 
	: TObject(), m_nCatalogs(nCatalogs)
{
	//Init source
	if(aSource){	
		if(cloneSource){
			m_source= new Source;
			*(m_source)= *(aSource);
			m_ownedSource= true;
		}	
		else{
			m_source= aSource;
		}
	}

	//Initialize 
	Init();

}//close costructor


SourceMatchData::~SourceMatchData()
{
	//Delete source (if pointer is owned)
	if(m_ownedSource){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting source ...");
		#endif
	
		CodeUtils::DeletePtr<Source>(m_source);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("done!");
		#endif
	}

	//Delete source collection
	if(m_ownedMatchedSources){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting source match collection added ...");
		#endif
		CodeUtils::DeletePtrCollection<SourceGroup>(m_matchedSources);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("done!");
		#endif
	}
	
	
	//Delete source SED graph
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source SED graph ...");
	#endif
	CodeUtils::DeletePtr<TGraphErrors>(m_sourceSED);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif

	//Delete source component SED collection
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source component SED graph added ...");
	#endif
	CodeUtils::DeletePtrCollection<TGraphErrors>(m_sourceComponentSED);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif

	//Delete source SED fit fcn
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Deleting source SED fit fcn ...");
	#endif
	CodeUtils::DeletePtr<TF1>(m_sedPLFitFcn);
	CodeUtils::DeletePtrCollection<TF1>(m_sedComponentPLFitFcns);
	CodeUtils::DeletePtr<TF1>(m_sedSSAFitFcn);
	CodeUtils::DeletePtrCollection<TF1>(m_sedComponentSSAFitFcns);
	CodeUtils::DeletePtr<TF1>(m_sedSIFFAFitFcn);
	CodeUtils::DeletePtrCollection<TF1>(m_sedComponentSIFFAFitFcns);
	CodeUtils::DeletePtr<TF1>(m_sedSEFFAFitFcn);
	CodeUtils::DeletePtrCollection<TF1>(m_sedComponentSEFFAFitFcns);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("done!");
	#endif
	
}//close destructor

/*
SourceMatchData::SourceMatchData(const SourceMatchData& sourceMatchData) 
{
  // Copy constructor
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy constructor called...");
	#endif
  Init();
  ((SourceMatchData&)sourceMatchData).Copy(*this);
	#ifdef LOGGING_ENABLED
		INFO_LOG("done!");
	#endif
}
*/



/*
void SourceMatchData::Copy(TObject &obj) const 
{
	//Copy object
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy object...");
	#endif
	TObject::Copy((SourceMatchData&)obj);

	//Copy vars
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy vars...");
	#endif
	((SourceMatchData&)obj).m_nCatalogs = m_nCatalogs;
	((SourceMatchData&)obj).m_nMatches = m_nMatches;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy spectral index data...");
	#endif
	((SourceMatchData&)obj).m_spectralIndexData = m_spectralIndexData;
 	(((SourceMatchData&)obj).m_componentSpectralIndexData).clear();
	for(size_t i=0;i<m_componentSpectralIndexData.size();i++){
		(((SourceMatchData&)obj).m_componentSpectralIndexData).push_back(m_componentSpectralIndexData[i]);
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy component match indexes...");
	#endif
	(((SourceMatchData&)obj).m_componentMatchIndexes).clear();
	for(size_t i=0;i<m_componentMatchIndexes.size();i++){
		(((SourceMatchData&)obj).m_componentMatchIndexes).push_back( std::vector<ComponentMatchIndexGroup>() );
		for(size_t j=0;j<m_componentMatchIndexes[i].size();j++) (((SourceMatchData&)obj).m_componentMatchIndexes)[i].push_back(m_componentMatchIndexes[i][j]);
	}

	//Copy source SED graph
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy source SED...");
	#endif
	if(((SourceMatchData&)obj).m_sourceSED){
		delete ((SourceMatchData&)obj).m_sourceSED;
		((SourceMatchData&)obj).m_sourceSED= 0;
	}
	if(m_sourceSED){
		((SourceMatchData&)obj).m_sourceSED= new TGraphErrors;
		*((SourceMatchData&)obj).m_sourceSED = *m_sourceSED;
	}

	//Copy source SED component collection
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy source component SED...");
	#endif
	//Delete first any existing collection
	for(size_t i=0;i<(((SourceMatchData&)obj).m_sourceComponentSED).size();i++){
		if( (((SourceMatchData&)obj).m_sourceComponentSED)[i] ){
			delete (((SourceMatchData&)obj).m_sourceComponentSED)[i];
			(((SourceMatchData&)obj).m_sourceComponentSED)[i]= 0;
		}
	}
	(((SourceMatchData&)obj).m_sourceComponentSED).clear();

	
	TGraphErrors* aGraph= 0;
	for(size_t i=0;i<m_sourceComponentSED.size();i++){
		aGraph= new TGraphErrors;
		*aGraph= *(m_sourceComponentSED[i]);
		(((SourceMatchData&)obj).m_sourceComponentSED).push_back(aGraph);
	}

	//Copy source
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy source ...");
	#endif
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
	#ifdef LOGGING_ENABLED
		INFO_LOG("Deleting existing matched source collection ...");
	#endif
	for(size_t i=0;i<(((SourceMatchData&)obj).m_matchedSources).size();i++){
		if( (((SourceMatchData&)obj).m_matchedSources)[i] ){
			delete (((SourceMatchData&)obj).m_matchedSources)[i];
			(((SourceMatchData&)obj).m_matchedSources)[i]= 0;
		}
	}
	(((SourceMatchData&)obj).m_matchedSources).clear();
	

	#ifdef LOGGING_ENABLED
		INFO_LOG("Copying matched source collection ...");
	#endif

	SourceGroup* aSourceGroup= 0;
	for(size_t i=0;i<m_matchedSources.size();i++){
		aSourceGroup= new SourceGroup;
		*aSourceGroup= *(m_matchedSources[i]);
		(((SourceMatchData&)obj).m_matchedSources).push_back(aSourceGroup);
	}
	
}//close Copy()
*/

/*
SourceMatchData& SourceMatchData::operator=(const SourceMatchData& sourceMatchData) 
{ 
	#ifdef LOGGING_ENABLED
		INFO_LOG("Copy operator called ...");
	#endif
	// Operator =
  if (this != &sourceMatchData) ((SourceMatchData&)sourceMatchData).Copy(*this);
  return *this;
}
*/

void SourceMatchData::Init()
{
	
	//Init source SEDs
	m_sourceSED= new TGraphErrors;
	
	
	//Init source group data
	m_nMatches= 0;
	m_nMatchesPerCatalog.clear();
	m_componentMatchIndexes.clear();
	m_sourceMatchPars.clear();
	m_sourceMatchPars.resize(m_nCatalogs);
	SourceGroup* aSourceGroup= 0;

	for(int i=0;i<m_nCatalogs;i++){
		m_nMatchesPerCatalog.push_back(0);
		aSourceGroup= new SourceGroup;
		m_matchedSources.push_back(aSourceGroup);
	}
	
	
	//Init component match index
	m_componentSpectralIndexData.clear();
	m_componentPolSpectralIndexData.clear();
	m_componentSSASpectralIndexData.clear();
	m_componentSIFFASpectralIndexData.clear();
	m_componentSEFFASpectralIndexData.clear();
	m_sourceComponentSED.clear();
	TGraphErrors* sed= 0;

	if(m_source){
		int nComponents= m_source->GetNFitComponents();
		if(nComponents>0) {
			for(int i=0;i<nComponents;i++){
				sed= new TGraphErrors;
				m_sourceComponentSED.push_back(sed);
				m_componentMatchIndexes.push_back( std::vector<ComponentMatchIndexGroup>() );
				m_componentMatchIndexes[i].resize(m_nCatalogs);

				m_sourceComponentMatchPars.push_back( std::vector<SourceComponentMatchParsGroup>() );	
				m_sourceComponentMatchPars[i].resize(m_nCatalogs);
			}
			m_componentSpectralIndexData.resize(nComponents);
			m_componentPolSpectralIndexData.resize(nComponents);
			m_componentSSASpectralIndexData.resize(nComponents);
			m_componentSIFFASpectralIndexData.resize(nComponents);
			m_componentSEFFASpectralIndexData.resize(nComponents);
		}
	}

	m_shiftFlux= false;
	m_fluxShift= 0.;

	//Init fit functions
	m_sedPLFitFcn= 0;
	m_sedComponentPLFitFcns.clear();	
	m_sedSSAFitFcn= 0;
	m_sedComponentSSAFitFcns.clear();
	m_sedSIFFAFitFcn= 0;
	m_sedComponentSIFFAFitFcns.clear();
	m_sedSEFFAFitFcn= 0;
	m_sedComponentSEFFAFitFcns.clear();

}//close Init()




bool SourceMatchData::HasSourceMatch()
{
	if(m_matchedSources.empty()) return false;
	bool hasMatch= false;
	for(size_t i=0;i<m_matchedSources.size();i++){
		if(m_matchedSources[i]->GetNSources()>0){
		//if(m_matchedSources[i].GetNSources()>0){
			hasMatch= true;
			break;
		}
	}
	return hasMatch;

}//close HasSourceMatch()


bool SourceMatchData::HasSourceComponentMatch(int componentId)
{
	if(m_componentMatchIndexes.empty()) return false;
	if(componentId<0 || componentId>=(int)(m_componentMatchIndexes.size())) return false;
	if(m_componentMatchIndexes[componentId].empty()) return false;
			
	bool hasMatch= false;
	for(size_t j=0;j<m_componentMatchIndexes[componentId].size();j++){
		int nIndexes= m_componentMatchIndexes[componentId][j].GetNIndexes();
		//std::vector<ComponentMatchIndex> index_list= m_componentMatchIndexes[componentId][j].GetIndexes();
		//if(!index_list.empty()){
		if(nIndexes>0){
			hasMatch= true;
			break;
		}
	}
	return hasMatch;

}//close HasSourceComponentMatch()



int SourceMatchData::GetMatchedSourcesPerCatalog(std::vector<Source*>& sources,int catalogIndex)
{
	if(m_matchedSources.empty()) return -1;
	if(catalogIndex<0 || catalogIndex>=(int)(m_matchedSources.size())) return -1;
	if(!m_matchedSources[catalogIndex]){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Matched source group (CAT="<<catalogIndex<<") is nullptr, this should not occur!");
		#endif
		return -1;
	}
	sources= m_matchedSources[catalogIndex]->GetSources();
				
	return 0;

}//close GetMatchedSourcesPerCatalog()





int SourceMatchData::AddMatchedSourceToGroup(int catalogIndex,Source* aSource,bool clone)
{
	if(!aSource) return -1;
	if(catalogIndex<0 || catalogIndex>=(int)(m_matchedSources.size())) return -1;
	if(!m_matchedSources[catalogIndex]){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Matched source group (CAT="<<catalogIndex<<") is nullptr, this should not occur!");
		#endif
		return -1;
	}

	//Add source
	int status= m_matchedSources[catalogIndex]->AddSource(aSource,clone);

	//Update number of matched sources
	m_nMatches= 0;
	for(size_t i=0;i<m_matchedSources.size();i++){
		int nMatchesInCatalog= m_matchedSources[i]->GetNSources();
		m_nMatchesPerCatalog[i]= nMatchesInCatalog;
		if(nMatchesInCatalog>0) m_nMatches++;
	}
	return status;	

}//close AddMatchedSourceToGroup()


int SourceMatchData::AddMatchedSourcePars(int catalogIndex,SourceMatchPars* smatchpars)
{
	if(!smatchpars) return -1;
	if(catalogIndex<0 || catalogIndex>=(int)(m_sourceMatchPars.size())) return -1;
	//m_sourceMatchPars[catalogIndex].push_back(*smatchpars);
	m_sourceMatchPars[catalogIndex].AddPars(smatchpars);

	return 0;

}//close AddMatchedSourcePars()


int SourceMatchData::AddMatchedSourceComponentPars(int componentIndex,int catalogIndex,SourceComponentMatchPars* cmatchpars)
{
	if(!cmatchpars) return -1;
	if(componentIndex<0 || componentIndex>=(int)(m_sourceComponentMatchPars.size())){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Component with index "<<componentIndex<<" was not allocated!");
		#endif
		return -1;
	}
	if(catalogIndex<0 || catalogIndex>=(int)(m_sourceComponentMatchPars[componentIndex].size())) return -1;

	m_sourceComponentMatchPars[componentIndex][catalogIndex].AddPars(cmatchpars);	

	return 0;

}//close	AddMatchedSourceComponentPars()	


int SourceMatchData::AddMatchedComponentIndex(int componentIndex,int catalogIndex,int matchedSourceIndex,int matchedComponentIndex)
{
	//Check if component id was allocated
	int nComponents= static_cast<int>(m_componentMatchIndexes.size());
	//if(nComponents<=0 || nComponents<componentIndex+1){
	if(nComponents<=0 || componentIndex<0 || componentIndex>=nComponents){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Component with index "<<componentIndex<<" was not allocated!");
		#endif
		return -1;
	}
	
	//Fill component indexes
	m_componentMatchIndexes[componentIndex][catalogIndex].AddIndex(catalogIndex,matchedSourceIndex,matchedComponentIndex);				

	return 0;
		
}//close AddMatchedComponentIndex()


int SourceMatchData::GetComponentCatalogMatchMultiplicity(int componentIndex)
{
	//Check if component id was allocated
	int nComponents= static_cast<int>(m_componentMatchIndexes.size());
	if(nComponents<=0 || componentIndex<0 || componentIndex>=nComponents){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Component with index "<<componentIndex<<" was not allocated!");
		#endif
		return -1;
	}

	//Loop over match index and count
	int nMatches= 0;
	for(size_t j=0;j<m_componentMatchIndexes[componentIndex].size();j++)
	{
		int nCatalogMatches= m_componentMatchIndexes[componentIndex][j].GetNIndexes();
		if(nCatalogMatches>0) nMatches++;

	}//end loop catalog

	return nMatches;

}//close GetComponentCatalogMatchMultiplicity()



int SourceMatchData::GetMatchedSourceNames(std::vector<std::string>& snames,int catalogIndex)
{
	snames.clear();
	if(catalogIndex<0 || catalogIndex>=(int)(m_matchedSources.size())) return -1;
	if(!m_matchedSources[catalogIndex]){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Matched source group (CAT="<<catalogIndex<<") is nullptr, this should not occur!");
		#endif
		return -1;
	}

	if(m_matchedSources[catalogIndex]->GetSourceNames(snames)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get matched source names!");
		#endif
		return -1;
	}
	
	return 0;

}//close GetMatchedSourceNames()



int SourceMatchData::GetMatchedSourceFrequency(int catalogIndex,double& freq,double& dfreq)
{
	freq= 0;
	dfreq= 0;	
	if(m_matchedSources.empty()) return -1;				
	if(catalogIndex<0 || catalogIndex>=(int)(m_matchedSources.size())) return -1;
	if(!m_matchedSources[catalogIndex]){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Matched source group (CAT="<<catalogIndex<<") is nullptr, this should not occur!");
		#endif
		return -1;
	}
			
	int status= m_matchedSources[catalogIndex]->GetFrequency(freq,dfreq);
			
	return status;
}



int SourceMatchData::GetMatchedSourceFlux(int catalogIndex,double& flux,double& fluxErr,bool& summed)
{
	flux= 0;
	fluxErr= 0;
	summed= false;
			
	if(m_matchedSources.empty()) return -1;
	if(catalogIndex<0 || catalogIndex>=(int)(m_matchedSources.size())) return -1;
	if(!m_matchedSources[catalogIndex]){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Matched source group (CAT="<<catalogIndex<<") is nullptr, this should not occur!");
		#endif
		return -1;
	}

	int status= m_matchedSources[catalogIndex]->GetFlux(flux,fluxErr,summed);
			
	return status;
}



int SourceMatchData::ComputeSourceSEDs()
{
	//Do nothing if no source set
	if(!m_source) return -1;
	if(m_matchedSources.empty()) return 0;
	if(!m_sourceSED) return -1;
	
	//Delete existing source SEDs
	#ifdef LOGGING_ENABLED
		INFO_LOG("Clear existing SED graph...");
	#endif
	TString graphName= Form("sourceSED_%s",m_source->GetName());
	m_sourceSED->SetNameTitle(graphName,graphName);
	m_sourceSED->Set(0);//remove all existing points	
	m_spectralIndexData.Reset();//reset spectral index data
	//m_spectralIndexData.isMultiSourceMatchIndex= false;
	//m_spectralIndexData.hasSpectralIndex= false;
	
	std::string sname=  m_source->GetName();

	//## Get source spectral info
	std::vector<double> lgNu_list;
	std::vector<double> lgFlux_list;
	std::vector<double> lgFluxErr_list;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Get source metadata ...");
	#endif

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
	if(m_shiftFlux){
		flux*= (1. + m_fluxShift);
		fluxErr*= (1. + m_fluxShift);
	}

	double lgFlux= log10(flux);
	double lgFluxErr= log10(TMath::E())*fluxErr/flux;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<sname<<": Nu="<<Nu<<", dNu="<<dNu<<", FreqUnits="<<FreqUnits<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
	#endif
	
	lgNu_list.push_back(lgNu);
	lgFlux_list.push_back(lgFlux);
	lgFluxErr_list.push_back(lgFluxErr);

	//Check if we have matched sources
	if(!HasSourceMatch()) {
		m_sourceSED->SetPoint(0,lgNu,lgFlux);
		m_sourceSED->SetPointError(0,0,lgFluxErr);
		return 0;//nothing to be done without matched sources
	}
	

	//## Get matched source info
	bool isMultiSourceMatchIndex= false;

	for(size_t i=0;i<m_matchedSources.size();i++)
	{
		//Skip if no matches
		if(m_matchedSources[i]->GetNSources()<=0) continue;
		
		//Get source names
		std::vector<std::string> snames;
		if(GetMatchedSourceNames(snames,i)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to get names of matched source group no. "<<i+1<<"!");
			#endif
		}
		std::stringstream ss;
		if(!snames.empty()){
			ss<<"{";
			for(size_t j=0;j<snames.size()-1;j++) ss<<snames[j]<<",";
			ss<<snames[snames.size()-1]<<"}";
		}

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
		//if(fluxSummed) m_spectralIndexData.isMultiSourceMatchIndex= true;
		if(fluxSummed) isMultiSourceMatchIndex= true;

		#ifdef LOGGING_ENABLED
			INFO_LOG("Matched source group no. "<<i+1<<": snames="<<ss.str()<<", Nu="<<Nu<<", dNu="<<dNu<<", lgNu="<<lgNu<<", flux="<<flux<<", lgFlux="<<lgFlux);
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
	
	for(int i=0;i<N;i++){
		m_sourceSED->SetPoint(i,lgNu_list_sorted[i],lgFlux_list_sorted[i]);
		m_sourceSED->SetPointError(i,0,lgFluxErr_list_sorted[i]);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Source SED: P"<<i+1<<"("<<lgNu_list_sorted[i]<<","<<lgFlux_list_sorted[i]<<")");
		#endif
	}


	//# Compute spectral index from power law fit
	m_spectralIndexData.isMultiSourceMatchIndex= isMultiSourceMatchIndex;
	FitPLToSourceSED(m_spectralIndexData,m_sourceSED);

	//# Fit pol3 to SED
	m_polSpectralIndexData.isMultiMatch= isMultiSourceMatchIndex;
	FitPolToSourceSED(m_polSpectralIndexData,m_sourceSED);

	//# Fit SSA to SED
	m_ssaSpectralIndexData.isMultiMatch= isMultiSourceMatchIndex;
	FitSSAToSourceSED(m_ssaSpectralIndexData,m_sourceSED);

	//# Fit SIFFA to SED
	m_siffaSpectralIndexData.isMultiMatch= isMultiSourceMatchIndex;
	FitSIFFAToSourceSED(m_siffaSpectralIndexData,m_sourceSED);

	//# Fit SEFFA to SED
	m_seffaSpectralIndexData.isMultiMatch= isMultiSourceMatchIndex;
	FitSEFFAToSourceSED(m_seffaSpectralIndexData,m_sourceSED);

	return 0;

}//close ComputeSourceSEDs()


int SourceMatchData::FitSourceSED(bool useRobustFitter,double rob)
{
	//## Check SED
	if(!m_sourceSED){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Source SED was not computed yet!");
		#endif
		return -1;
	}
	if(m_sourceSED->GetN()<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source SED has no points, nothing to be done...");
		#endif
		return 0;
	}

	//## Fit SED
	//- Fit power law
	FitPLToSourceSED(m_spectralIndexData,m_sourceSED,useRobustFitter,rob);

	//- Fit pol3 
	FitPolToSourceSED(m_polSpectralIndexData,m_sourceSED,useRobustFitter,rob);

	//- Fit SSA 
	FitSSAToSourceSED(m_ssaSpectralIndexData,m_sourceSED,useRobustFitter,rob);

	//- Fit SIFFA 
	FitSIFFAToSourceSED(m_siffaSpectralIndexData,m_sourceSED,useRobustFitter,rob);

	//- Fit SEFFA
	FitSEFFAToSourceSED(m_seffaSpectralIndexData,m_sourceSED,useRobustFitter,rob);

	return 0;

}//close FitSourceSED()


int SourceMatchData::FitPLToSourceSED(SpectralIndexData& spectralIndexData,TGraphErrors* sed,bool useRobustFitter,double rob)
{
	//- Find spectral index
	int N= sed->GetN();
	
	if(N==2){
		spectralIndexData.isSpectralIndexFit= false;
		double x1= 0; 
		double x2= 0;	
		double y1= 0;
		double y2= 0;	
		double x1_err= 0;
		double x2_err= 0;
		double y1_err= 0;
		double y2_err= 0;

		sed->GetPoint(0,x1,y1);
		x1_err= sed->GetErrorX(0);
		y1_err= sed->GetErrorY(0);	
		sed->GetPoint(1,x2,y2);
		x2_err= sed->GetErrorX(1);
		y2_err= sed->GetErrorY(1);

		double y2_up= y2 + y2_err;
		double y2_down= y2 - y2_err;
		double y1_up= y1 + y1_err;
		double y1_down= y1 - y1_err;
		double alpha= (y2-y1)/(x2-x1);
		double alpha_up= (y2_up-y1_down)/(x2-x1);
		double alpha_down= (y2_down-y1_up)/(x2-x1);
		double dalpha_up= alpha_up-alpha;
		double dalpha_down= alpha-alpha_down;		
		double alphaErr= fabs(dalpha_up);		
		double norm= (x2*y1-x1*y2)/(x2-x1);
		double normErr= sqrt( pow(x2*y1_err/(x2-x1),2) + pow(x1*y2_err/(x2-x1),2) );

		spectralIndexData.spectralIndex= alpha;
		spectralIndexData.spectralIndexErr= alphaErr;
		spectralIndexData.norm= norm;	
		spectralIndexData.normErr= normErr;		
		spectralIndexData.spectralFitChi2= 0;
		spectralIndexData.spectralFitNDF= 0;
		spectralIndexData.hasSpectralIndex= true;
	}
	else{
		spectralIndexData.isSpectralIndexFit= true;

		TString fitOptions= "RMS";
		if(useRobustFitter) fitOptions= Form("SN rob=%1.2f",rob);
		TFitResultPtr fitRes= sed->Fit("pol1",fitOptions);
		spectralIndexData.spectralFitChi2= fitRes->Chi2();
		spectralIndexData.spectralFitNDF= fitRes->Ndf();
		spectralIndexData.norm= fitRes->Value(0);
		spectralIndexData.normErr= fitRes->ParError(0);
		spectralIndexData.spectralIndex= fitRes->Value(1);
		spectralIndexData.spectralIndexErr= fitRes->ParError(1);
		spectralIndexData.spectralIndexFitStatus= fitRes->Status();
		spectralIndexData.hasSpectralIndex= true;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<": hasData?"<<spectralIndexData.hasSpectralIndex<<", alpha="<<spectralIndexData.spectralIndex<<" +- "<<spectralIndexData.spectralIndexErr<<", chi2/ndf="<<spectralIndexData.spectralFitChi2<<"/"<<spectralIndexData.spectralFitNDF);
	#endif
	

	return 0;

}//close FitPLToSourceSED()


int SourceMatchData::FitPolToSourceSED(PolSpectralIndexData& spectralIndexData,TGraphErrors* sed,bool useRobustFitter,double rob)
{
	//# Nothing to be done if less than N=4 points present in SED
	int N= sed->GetN();
	spectralIndexData.hasData= false;
	if(N<4) return 0; 
	
	//# Fit sed
	TString fitOptions= "RMS";
	if(useRobustFitter) fitOptions= Form("SN rob=%1.2f",rob);
	TFitResultPtr fitRes= sed->Fit("pol3",fitOptions);
	spectralIndexData.fitChi2= fitRes->Chi2();
	spectralIndexData.fitNDF= fitRes->Ndf();
	std::stringstream ss;
	ss<<"polPars{";
	for(int i=0;i<fitRes->NPar();i++){
		double p= fitRes->Value(i);
		double pErr= fitRes->ParError(i);
		(spectralIndexData.polPars).push_back(p);
		(spectralIndexData.polParErrors).push_back(pErr);
		ss<<"p"<<i+1<<"="<<p<<"+-"<<pErr<<",";
	}
	ss<<"}";

	spectralIndexData.status= fitRes->Status();
	spectralIndexData.hasData= true;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<": hasData?"<<spectralIndexData.hasData<<", "<<ss.str()<<", chi2/ndf="<<spectralIndexData.fitChi2<<"/"<<spectralIndexData.fitNDF);
	#endif

	
	
	return 0;

}//close FitPolToSourceSED()


int SourceMatchData::FitSSAToSourceSED(SASpectralIndexData& spectralIndexData,TGraphErrors* sed,bool useRobustFitter,double rob)
{
	//# Nothing to be done if less than N=4 points present in SED
	int N= sed->GetN();
	spectralIndexData.hasData= false;
	if(N<6) return 0; 
	
	//# Create fit function
	double xmin= TMath::MinElement(sed->GetN(),sed->GetX());
	double xmax= TMath::MaxElement(sed->GetN(),sed->GetX());
	double dx= fabs(xmax-xmin);
	double xmin_fcn= xmin - 0.2*dx;
	double xmax_fcn= xmax + 0.2*dx;
	double ymax= TMath::MaxElement(sed->GetN(),sed->GetY());
	double norm_start= ymax;
	double alpha_start= -0.5;
	double nu_t_start= xmin + 0.5*(xmax-xmin);
	TF1* fitFcn= new TF1("",MathUtils::SynchrotronSelfAbsSED,xmin_fcn,xmax_fcn,3);
	fitFcn->SetParameters(ymax,alpha_start,nu_t_start);

	//# Fit sed
	TString fitOptions= "RMS";
	if(useRobustFitter) fitOptions= Form("SN rob=%1.2f",rob);
	TFitResultPtr fitRes= sed->Fit(fitFcn,fitOptions);
	spectralIndexData.fitChi2= fitRes->Chi2();
	spectralIndexData.fitNDF= fitRes->Ndf();
	spectralIndexData.norm= fitRes->Value(0);
	spectralIndexData.normErr= fitRes->ParError(0);
	spectralIndexData.alpha= fitRes->Value(1);
	spectralIndexData.alphaErr= fitRes->ParError(1);
	spectralIndexData.nu_t= fitRes->Value(2);
	spectralIndexData.nuErr_t= fitRes->ParError(2);
	spectralIndexData.status= fitRes->Status();
	spectralIndexData.hasData= true;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<" SSA Fit: hasData?"<<spectralIndexData.hasData<<", alpha="<<spectralIndexData.alpha<<"+-"<<spectralIndexData.alphaErr<<", nu_t="<<spectralIndexData.nu_t<<"+-"<<spectralIndexData.nuErr_t<<", chi2/ndf="<<spectralIndexData.fitChi2<<"/"<<spectralIndexData.fitNDF);
	#endif

	//# Clear data
	delete fitFcn;
	fitFcn= 0;
	
	return 0;

}//close FitSSAToSourceSED()


int SourceMatchData::FitSIFFAToSourceSED(SASpectralIndexData& spectralIndexData,TGraphErrors* sed,bool useRobustFitter,double rob)
{
	//# Nothing to be done if less than N=4 points present in SED
	int N= sed->GetN();
	spectralIndexData.hasData= false;
	if(N<6) return 0; 
	
	//# Create fit function
	double xmin= TMath::MinElement(sed->GetN(),sed->GetX());
	double xmax= TMath::MaxElement(sed->GetN(),sed->GetX());
	double dx= fabs(xmax-xmin);
	double xmin_fcn= xmin - 0.2*dx;
	double xmax_fcn= xmax + 0.2*dx;
	double ymax= TMath::MaxElement(sed->GetN(),sed->GetY());
	double norm_start= ymax;
	double alpha_start= -0.5;
	double nu_t_start= xmin + 0.5*(xmax-xmin);
	TF1* fitFcn= new TF1("",MathUtils::SynchrotronIntFreeFreeAbsSED,xmin_fcn,xmax_fcn,3);
	fitFcn->SetParameters(ymax,alpha_start,nu_t_start);

	//# Fit sed
	TString fitOptions= "RMS";
	if(useRobustFitter) fitOptions= Form("SN rob=%1.2f",rob);
	TFitResultPtr fitRes= sed->Fit(fitFcn,fitOptions);
	spectralIndexData.fitChi2= fitRes->Chi2();
	spectralIndexData.fitNDF= fitRes->Ndf();
	spectralIndexData.norm= fitRes->Value(0);
	spectralIndexData.normErr= fitRes->ParError(0);
	spectralIndexData.alpha= fitRes->Value(1);
	spectralIndexData.alphaErr= fitRes->ParError(1);
	spectralIndexData.nu_t= fitRes->Value(2);
	spectralIndexData.nuErr_t= fitRes->ParError(2);
	spectralIndexData.status= fitRes->Status();
	spectralIndexData.hasData= true;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<" SIFFA Fit: hasData?"<<spectralIndexData.hasData<<", alpha="<<spectralIndexData.alpha<<"+-"<<spectralIndexData.alphaErr<<", nu_t="<<spectralIndexData.nu_t<<"+-"<<spectralIndexData.nuErr_t<<", chi2/ndf="<<spectralIndexData.fitChi2<<"/"<<spectralIndexData.fitNDF);
	#endif

	//# Clear data
	delete fitFcn;
	fitFcn= 0;
	
	return 0;

}//close FitSIFFAToSourceSED()

int SourceMatchData::FitSEFFAToSourceSED(SASpectralIndexData& spectralIndexData,TGraphErrors* sed,bool useRobustFitter,double rob)
{
	//# Nothing to be done if less than N=4 points present in SED
	int N= sed->GetN();
	spectralIndexData.hasData= false;
	if(N<6) return 0; 
	
	//# Create fit function
	double xmin= TMath::MinElement(sed->GetN(),sed->GetX());
	double xmax= TMath::MaxElement(sed->GetN(),sed->GetX());
	double dx= fabs(xmax-xmin);
	double xmin_fcn= xmin - 0.2*dx;
	double xmax_fcn= xmax + 0.2*dx;
	double ymax= TMath::MaxElement(sed->GetN(),sed->GetY());
	double norm_start= ymax;
	double alpha_start= -0.5;
	double nu_t_start= xmin + 0.5*(xmax-xmin);
	TF1* fitFcn= new TF1("",MathUtils::SynchrotronExtFreeFreeAbsSED,xmin_fcn,xmax_fcn,3);
	fitFcn->SetParameters(ymax,alpha_start,nu_t_start);

	//# Fit sed
	TString fitOptions= "RMS";
	if(useRobustFitter) fitOptions= Form("SN rob=%1.2f",rob);
	TFitResultPtr fitRes= sed->Fit(fitFcn,fitOptions);
	spectralIndexData.fitChi2= fitRes->Chi2();
	spectralIndexData.fitNDF= fitRes->Ndf();
	spectralIndexData.norm= fitRes->Value(0);
	spectralIndexData.normErr= fitRes->ParError(0);
	spectralIndexData.alpha= fitRes->Value(1);
	spectralIndexData.alphaErr= fitRes->ParError(1);
	spectralIndexData.nu_t= fitRes->Value(2);
	spectralIndexData.nuErr_t= fitRes->ParError(2);
	spectralIndexData.status= fitRes->Status();
	spectralIndexData.hasData= true;

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source "<<m_source->GetName()<<" SEFFA Fit: hasData?"<<spectralIndexData.hasData<<", alpha="<<spectralIndexData.alpha<<"+-"<<spectralIndexData.alphaErr<<", nu_t="<<spectralIndexData.nu_t<<"+-"<<spectralIndexData.nuErr_t<<", chi2/ndf="<<spectralIndexData.fitChi2<<"/"<<spectralIndexData.fitNDF);
	#endif

	//# Clear data
	delete fitFcn;
	fitFcn= 0;
	
	return 0;

}//close FitSEFFAToSourceSED()

int SourceMatchData::ComputeSourceComponentSEDs()
{
	//Do nothing if no source set
	if(!m_source) return -1;
	if(m_matchedSources.empty()) return 0;
	if(m_componentMatchIndexes.empty()) return 0;

	//Get source pars
	if(!m_source->HasFitInfo()) return 0;//no components, nothing to be done

	//Get source fit component pars
	SourceFitPars fitPars= m_source->GetFitPars();
	int nComponents= m_source->GetNFitComponents();
	int nComponents_index= static_cast<int>(m_componentMatchIndexes.size());

	if(nComponents!=nComponents_index){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Mismatch between number of components found in source and set in class (nComponents="<<nComponents<<"!="<<nComponents_index<<")!");	
		#endif	
		return -1;
	}
	
	std::string sourceName= m_source->GetName();

	#ifdef LOGGING_ENABLED
		INFO_LOG("sourceName="<<sourceName<<", nComponents="<<nComponents_index);	
	#endif	

	
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


	//Reset existing source SEDs
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reset existing SED graphs...");
	#endif
	for(size_t i=0;i<m_sourceComponentSED.size();i++){
		TString graphName= Form("sourceSED_%s_fitcomp%d",m_source->GetName(),(int)(i+1));
		m_sourceComponentSED[i]->SetNameTitle(graphName,graphName);
		m_sourceComponentSED[i]->Set(0);
	}

	
	//Loop over source components
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
		if(m_shiftFlux){
			flux*= (1. + m_fluxShift);
			fluxErr*= (1. + m_fluxShift);
		}
		double lgFlux= log10(flux);
		double lgFluxErr= log10(TMath::E())*fluxErr/flux;
		
		lgNu_list.push_back(lgNu);
		lgFlux_list.push_back(lgFlux);
		lgFluxErr_list.push_back(lgFluxErr);
		bool isMatchComponentFluxSummed= false;	

		//Skip if no component match	
		if(!HasSourceComponentMatch(i)){
			m_sourceComponentSED[i]->SetPoint(0,lgNu,lgFlux);
			m_sourceComponentSED[i]->SetPointError(0,0,lgFluxErr);
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
				Source* matchedSource= m_matchedSources[catalogIndex]->GetSource(matchedSourceIndex);
				
				if(!matchedSource){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Matched source (CAT="<<catalogIndex<<", matchedSourceIndex="<<matchedSourceIndex<<", matchedComponentIndex="<<matchedComponentIndex<<") is null ptr, this should not occur!");
					#endif
					continue;
				}

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
		
		for(int k=0;k<N;k++){
			m_sourceComponentSED[i]->SetPoint(k,lgNu_list_sorted[k],lgFlux_list_sorted[k]);
			m_sourceComponentSED[i]->SetPointError(k,0,lgFluxErr_list_sorted[k]);
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source SED: P"<<k+1<<"("<<lgNu_list_sorted[k]<<","<<lgFlux_list_sorted[k]<<")");
			#endif
		}
		
		//# Find spectral index
		//# Compute spectral index from power law fit
		m_componentSpectralIndexData[i].isMultiSourceMatchIndex= isMatchComponentFluxSummed;	
		FitPLToSourceSED(m_componentSpectralIndexData[i],m_sourceComponentSED[i]);

		//# Fit pol3 to SED
		m_componentPolSpectralIndexData[i].isMultiMatch= isMatchComponentFluxSummed;
		FitPolToSourceSED(m_componentPolSpectralIndexData[i],m_sourceComponentSED[i]);

		//# Fit SSA to SED
		m_componentSSASpectralIndexData[i].isMultiMatch= isMatchComponentFluxSummed;
		FitSSAToSourceSED(m_componentSSASpectralIndexData[i],m_sourceComponentSED[i]);

		//# Fit SIFFA to SED
		m_componentSIFFASpectralIndexData[i].isMultiMatch= isMatchComponentFluxSummed;
		FitSIFFAToSourceSED(m_componentSIFFASpectralIndexData[i],m_sourceComponentSED[i]);

		//# Fit SEFFA to SED
		m_componentSEFFASpectralIndexData[i].isMultiMatch= isMatchComponentFluxSummed;
		FitSEFFAToSourceSED(m_componentSEFFASpectralIndexData[i],m_sourceComponentSED[i]);
	
	}//end loop components
	

	return 0;

}//close ComputeSourceComponentSEDs()

int SourceMatchData::FitSourceComponentSED(bool useRobustFitter,double rob)
{
	//Loop over component SEDs
	for(size_t i=0;i<m_sourceComponentSED.size();i++)
	{
		//Check SED
		if(!m_sourceComponentSED[i]) continue;
		if(m_sourceComponentSED[i]->GetN()<=0) continue;

		//Fit power law fit
		FitPLToSourceSED(m_componentSpectralIndexData[i],m_sourceComponentSED[i],useRobustFitter,rob);

		//Fit pol3 to SED
		FitPolToSourceSED(m_componentPolSpectralIndexData[i],m_sourceComponentSED[i],useRobustFitter,rob);

		//Fit SSA to SED
		FitSSAToSourceSED(m_componentSSASpectralIndexData[i],m_sourceComponentSED[i],useRobustFitter,rob);

		// Fit SIFFA to SED
		FitSIFFAToSourceSED(m_componentSIFFASpectralIndexData[i],m_sourceComponentSED[i],useRobustFitter,rob);

		//# Fit SEFFA to SED
		FitSEFFAToSourceSED(m_componentSEFFASpectralIndexData[i],m_sourceComponentSED[i],useRobustFitter,rob);

	}//end loop components

	return 0;

}//close FitSourceComponentSED()


TCanvas* SourceMatchData::DrawSED()
{
	//Do nothing if no source is present
	if(!m_source || !m_sourceSED){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No source and/or sed present, nothing to be drawn.");
		#endif	
		return nullptr;
	}	
	if(m_sourceSED->GetN()<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No points stored in source SED, nothing to be drawn.");
		#endif	
		return nullptr;
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
	flux*= 1.e+3;//convert to mJy
	std::string fluxUnits= "mJy/beam";
	if(beamArea>0) {
		flux/= beamArea;
		fluxUnits= "mJy";
	}

	if(m_shiftFlux){
		flux*= (1.+m_fluxShift);
	}
	
	fluxText= Form("%1.2f %s",flux,(fluxUnits).c_str());

	//Get spectral index
	TString sindexText= "";
	if(m_spectralIndexData.hasSpectralIndex){
	//if(m_spectralIndexData.hasData){
		double alpha= m_spectralIndexData.spectralIndex;
		double alphaErr= m_spectralIndexData.spectralIndexErr;
		sindexText= Form("alpha=%1.2f#pm%1.2f",alpha,alphaErr);
	}

	//Get component spectral index
	std::vector<TString> sindexTexts;
	for(size_t i=0;i<m_componentSpectralIndexData.size();i++){
		if(!m_componentSpectralIndexData[i].hasSpectralIndex) continue;
		double alpha= m_componentSpectralIndexData[i].spectralIndex;
		double alphaErr= m_componentSpectralIndexData[i].spectralIndexErr;
		TString s= Form("alpha_comp%d=%1.2f#pm%1.2f",(int)(i+1),alpha,alphaErr);
		sindexTexts.push_back(s);
	}

	//Create canvas
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadLeftMargin(0.15);

	TString canvasName= Form("SEDPlot_%s",sname.c_str());
	TCanvas* Plot= new TCanvas(canvasName,canvasName,700,700);
	Plot->cd();

	double step= 0.1;
	double x_min= lgNu_min - fabs(step*(lgNu_min-lgNu_max));
	double x_max= lgNu_max + fabs(step*(lgNu_min-lgNu_max));
	double y_min= lgFlux_min - fabs(step*(lgFlux_min-lgFlux_max));
	double y_max= lgFlux_max + fabs(step*(lgFlux_min-lgFlux_max));

	TString histoName= Form("SEDPlotBkg_%s",sname.c_str());
	TH2D* PlotBkg= new TH2D(histoName,"",100,x_min,x_max,100,y_min,y_max);
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
	if(sindexText!="") sourceInfoText->AddText(Form("%s",sindexText.Data()));
	for(size_t i=0;i<sindexTexts.size();i++){
		if(sindexTexts[i]!="") sourceInfoText->AddText(Form("%s",sindexTexts[i].Data()));
	}
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

	//## Draw SED fit
	// - Power law fit
	if(m_sedPLFitFcn){
		delete m_sedPLFitFcn;
		m_sedPLFitFcn= 0;
	}
	TString sedLPFitLegendText= "";
	if(m_spectralIndexData.hasSpectralIndex && m_spectralIndexData.isSpectralIndexFit)
	{	
		TString fcnName= Form("sedPLFitFcn_%s",m_source->GetName());
		m_sedPLFitFcn= new TF1(fcnName,"pol1",x_min,x_max);	
		m_sedPLFitFcn->SetParameter(0,m_spectralIndexData.norm);
		m_sedPLFitFcn->SetParameter(1,m_spectralIndexData.spectralIndex);
		m_sedPLFitFcn->SetLineColor(kBlack);
		m_sedPLFitFcn->SetLineStyle(kSolid);
		m_sedPLFitFcn->Draw("l same");

		sedLPFitLegendText= Form("sed PL fit (alpha=%1.2f#pm%1.2f)",m_spectralIndexData.spectralIndex,m_spectralIndexData.spectralIndexErr);

	}//close if

	// - SSA fit
	if(m_sedSSAFitFcn){
		delete m_sedSSAFitFcn;
		m_sedSSAFitFcn= 0;
	}
	TString sedSSAFitLegendText= "";
	if(m_ssaSpectralIndexData.hasData)
	{
		TString fcnName= Form("sedSSAFitFcn_%s",m_source->GetName());
		m_sedSSAFitFcn= new TF1(fcnName,MathUtils::SynchrotronSelfAbsSED,x_min,x_max);	
		m_sedSSAFitFcn->SetParameter(0,m_ssaSpectralIndexData.norm);
		m_sedSSAFitFcn->SetParameter(1,m_ssaSpectralIndexData.alpha);		
		m_sedSSAFitFcn->SetParameter(2,m_ssaSpectralIndexData.nu_t);
		m_sedSSAFitFcn->SetLineColor(kBlack);
		m_sedSSAFitFcn->SetLineStyle(kDashed);
		m_sedSSAFitFcn->Draw("l same");

		sedSSAFitLegendText= Form("sed SSA fit (alpha=%1.2f#pm%1.2f,nu_t=%1.2f#pm%1.2f)",m_ssaSpectralIndexData.alpha,m_ssaSpectralIndexData.alphaErr,m_ssaSpectralIndexData.nu_t,m_ssaSpectralIndexData.nuErr_t);
	}//close if

	// - SIFFA fit
	if(m_sedSIFFAFitFcn){
		delete m_sedSIFFAFitFcn;
		m_sedSIFFAFitFcn= 0;
	}
	TString sedSIFFAFitLegendText= "";
	if(m_siffaSpectralIndexData.hasData)
	{
		TString fcnName= Form("sedSIFFAFitFcn_%s",m_source->GetName());
		m_sedSIFFAFitFcn= new TF1("",MathUtils::SynchrotronIntFreeFreeAbsSED,x_min,x_max);	
		m_sedSIFFAFitFcn->SetParameter(0,m_siffaSpectralIndexData.norm);
		m_sedSIFFAFitFcn->SetParameter(1,m_siffaSpectralIndexData.alpha);		
		m_sedSIFFAFitFcn->SetParameter(2,m_siffaSpectralIndexData.nu_t);
		m_sedSIFFAFitFcn->SetLineColor(kBlack);
		m_sedSIFFAFitFcn->SetLineStyle(kDotted);
		m_sedSIFFAFitFcn->Draw("l same");

		sedSIFFAFitLegendText= Form("sed SIFFA fit (alpha=%1.2f#pm%1.2f,nu_t=%1.2f#pm%1.2f)",m_siffaSpectralIndexData.alpha,m_siffaSpectralIndexData.alphaErr,m_siffaSpectralIndexData.nu_t,m_siffaSpectralIndexData.nuErr_t);
	}//close if

	// - SEFFA fit
	if(m_sedSEFFAFitFcn){
		delete m_sedSEFFAFitFcn;
		m_sedSEFFAFitFcn= 0;
	}
	TString sedSEFFAFitLegendText= "";
	if(m_seffaSpectralIndexData.hasData)
	{	
		TString fcnName= Form("sedSEFFAFitFcn_%s",m_source->GetName());
		m_sedSEFFAFitFcn= new TF1("",MathUtils::SynchrotronExtFreeFreeAbsSED,x_min,x_max);	
		m_sedSEFFAFitFcn->SetParameter(0,m_seffaSpectralIndexData.norm);
		m_sedSEFFAFitFcn->SetParameter(1,m_seffaSpectralIndexData.alpha);		
		m_sedSEFFAFitFcn->SetParameter(2,m_seffaSpectralIndexData.nu_t);
		m_sedSEFFAFitFcn->SetLineColor(kBlack);
		m_sedSEFFAFitFcn->SetLineStyle(9);
		m_sedSEFFAFitFcn->Draw("l same");

		sedSEFFAFitLegendText= Form("sed SEFFA fit (alpha=%1.2f#pm%1.2f,nu_t=%1.2f#pm%1.2f)",m_seffaSpectralIndexData.alpha,m_seffaSpectralIndexData.alphaErr,m_seffaSpectralIndexData.nu_t,m_seffaSpectralIndexData.nuErr_t);
	}//close if

	
	//## Draw SED component fit
	// - Power law fit
	TF1* sedComponentPLFitFcn= 0;
	std::vector<TString> sedComponentPLFitLegendTexts;

	for(size_t i=0;i<m_componentSpectralIndexData.size();i++){	
		if(!m_componentSpectralIndexData[i].hasSpectralIndex || !m_componentSpectralIndexData[i].isSpectralIndexFit) continue;
		int color= kBlack;
		if(m_componentSpectralIndexData.size()<componentColors.size()){
			color= componentColors[i];
		}
		
		TString fcnName= Form("sedComponentPLFitFcn_%s_comp%d",m_source->GetName(),(int)(i+1));
		sedComponentPLFitFcn= new TF1(fcnName,"pol1",x_min,x_max);	
		sedComponentPLFitFcn->SetParameter(0,m_componentSpectralIndexData[i].norm);
		sedComponentPLFitFcn->SetParameter(1,m_componentSpectralIndexData[i].spectralIndex);
		sedComponentPLFitFcn->SetLineColor(color);
		sedComponentPLFitFcn->SetLineStyle(kSolid);
		sedComponentPLFitFcn->Draw("l same");
		m_sedComponentPLFitFcns.push_back(sedComponentPLFitFcn);

		TString legendText= Form("sedComp%d PL fit (alpha=%1.2f#pm%1.2f)",(int)(i+1),m_componentSpectralIndexData[i].spectralIndex,m_componentSpectralIndexData[i].spectralIndexErr);
		sedComponentPLFitLegendTexts.push_back(legendText);
	}//end loop components	

	// - SSA fit
	TF1* sedComponentSSAFitFcn= 0;
	std::vector<TString> sedComponentSSAFitLegendTexts;

	for(size_t i=0;i<m_componentSSASpectralIndexData.size();i++){	
		if(!m_componentSSASpectralIndexData[i].hasData) continue;
		int color= kBlack;
		if(m_componentSSASpectralIndexData.size()<componentColors.size()){
			color= componentColors[i];
		}
		
		TString fcnName= Form("sedComponentSSAFitFcn_%s_comp%d",m_source->GetName(),(int)(i+1));
		sedComponentSSAFitFcn= new TF1(fcnName,MathUtils::SynchrotronSelfAbsSED,x_min,x_max);	
		sedComponentSSAFitFcn->SetParameter(0,m_componentSSASpectralIndexData[i].norm);
		sedComponentSSAFitFcn->SetParameter(1,m_componentSSASpectralIndexData[i].alpha);
		sedComponentSSAFitFcn->SetParameter(2,m_componentSSASpectralIndexData[i].nu_t);
		sedComponentSSAFitFcn->SetLineColor(color);
		sedComponentSSAFitFcn->SetLineStyle(kDashed);
		sedComponentSSAFitFcn->Draw("l same");
		m_sedComponentSSAFitFcns.push_back(sedComponentSSAFitFcn);

		TString legendText= Form("sedComp%d SSA fit (alpha=%1.2f#pm%1.2f, nu_t=%1.2f#pm%1.2f)",(int)(i+1),m_componentSSASpectralIndexData[i].alpha,m_componentSSASpectralIndexData[i].alphaErr,m_componentSSASpectralIndexData[i].nu_t,m_componentSSASpectralIndexData[i].nuErr_t);
		sedComponentSSAFitLegendTexts.push_back(legendText);
	}//end loop components

	// - SIFFA fit
	TF1* sedComponentSIFFAFitFcn= 0;
	std::vector<TString> sedComponentSIFFAFitLegendTexts;

	for(size_t i=0;i<m_componentSIFFASpectralIndexData.size();i++){	
		if(!m_componentSIFFASpectralIndexData[i].hasData) continue;
		int color= kBlack;
		if(m_componentSIFFASpectralIndexData.size()<componentColors.size()){
			color= componentColors[i];
		}
		
		TString fcnName= Form("sedComponentSIFFAFitFcn_%s_comp%d",m_source->GetName(),(int)(i+1));
		sedComponentSIFFAFitFcn= new TF1(fcnName,MathUtils::SynchrotronIntFreeFreeAbsSED,x_min,x_max);	
		sedComponentSIFFAFitFcn->SetParameter(0,m_componentSIFFASpectralIndexData[i].norm);
		sedComponentSIFFAFitFcn->SetParameter(1,m_componentSIFFASpectralIndexData[i].alpha);
		sedComponentSIFFAFitFcn->SetParameter(2,m_componentSIFFASpectralIndexData[i].nu_t);
		sedComponentSIFFAFitFcn->SetLineColor(color);
		sedComponentSIFFAFitFcn->SetLineStyle(kDotted);
		sedComponentSIFFAFitFcn->Draw("l same");
		m_sedComponentSIFFAFitFcns.push_back(sedComponentSIFFAFitFcn);

		TString legendText= Form("sedComp%d SIFFA fit (alpha=%1.2f#pm%1.2f, nu_t=%1.2f#pm%1.2f)",(int)(i+1),m_componentSIFFASpectralIndexData[i].alpha,m_componentSIFFASpectralIndexData[i].alphaErr,m_componentSIFFASpectralIndexData[i].nu_t,m_componentSIFFASpectralIndexData[i].nuErr_t);
		sedComponentSIFFAFitLegendTexts.push_back(legendText);
	}//end loop components	
	
	// - SEFFA fit
	TF1* sedComponentSEFFAFitFcn= 0;
	std::vector<TString> sedComponentSEFFAFitLegendTexts;

	for(size_t i=0;i<m_componentSEFFASpectralIndexData.size();i++){	
		if(!m_componentSEFFASpectralIndexData[i].hasData) continue;
		int color= kBlack;
		if(m_componentSEFFASpectralIndexData.size()<componentColors.size()){
			color= componentColors[i];
		}
		
		TString fcnName= Form("sedComponentSEFFAFitFcn_%s_comp%d",m_source->GetName(),(int)(i+1));
		sedComponentSEFFAFitFcn= new TF1(fcnName,MathUtils::SynchrotronExtFreeFreeAbsSED,x_min,x_max);	
		sedComponentSEFFAFitFcn->SetParameter(0,m_componentSEFFASpectralIndexData[i].norm);
		sedComponentSEFFAFitFcn->SetParameter(1,m_componentSEFFASpectralIndexData[i].alpha);
		sedComponentSEFFAFitFcn->SetParameter(2,m_componentSEFFASpectralIndexData[i].nu_t);
		sedComponentSEFFAFitFcn->SetLineColor(color);
		sedComponentSEFFAFitFcn->SetLineStyle(kDotted);
		sedComponentSEFFAFitFcn->Draw("l same");
		m_sedComponentSEFFAFitFcns.push_back(sedComponentSEFFAFitFcn);

		TString legendText= Form("sedComp%d SEFFA fit (alpha=%1.2f#pm%1.2f, nu_t=%1.2f#pm%1.2f)",(int)(i+1),m_componentSEFFASpectralIndexData[i].alpha,m_componentSEFFASpectralIndexData[i].alphaErr,m_componentSEFFASpectralIndexData[i].nu_t,m_componentSEFFASpectralIndexData[i].nuErr_t);
		sedComponentSEFFAFitLegendTexts.push_back(legendText);
	}//end loop components	

	//## Draw legend
	TLegend* PlotLegend= new TLegend(0.6,0.6,0.8,0.8);
	PlotLegend->SetBorderSize(0);
	PlotLegend->SetFillColor(0);
	PlotLegend->SetTextFont(52);
	PlotLegend->SetTextSize(0.04);
	if(m_sourceSED) PlotLegend->AddEntry(m_sourceSED,"island","PL");
	for(size_t i=0;i<m_sourceComponentSED.size();i++){	
		if(m_sourceComponentSED[i]) PlotLegend->AddEntry(m_sourceComponentSED[i],Form("component %d",i+1),"PL");
	}
	if(m_sedPLFitFcn) PlotLegend->AddEntry(m_sedPLFitFcn,sedLPFitLegendText,"L");	
	if(m_sedSSAFitFcn) PlotLegend->AddEntry(m_sedSSAFitFcn,sedSSAFitLegendText,"L");
	if(m_sedSIFFAFitFcn) PlotLegend->AddEntry(m_sedSIFFAFitFcn,sedSIFFAFitLegendText,"L");
	if(m_sedSEFFAFitFcn) PlotLegend->AddEntry(m_sedSEFFAFitFcn,sedSEFFAFitLegendText,"L");
	for(size_t i=0;i<m_sedComponentPLFitFcns.size();i++){
		if(m_sedComponentPLFitFcns[i]) PlotLegend->AddEntry(m_sedComponentPLFitFcns[i],sedComponentPLFitLegendTexts[i],"L");	
	}
	for(size_t i=0;i<m_sedComponentSSAFitFcns.size();i++){
		if(m_sedComponentSSAFitFcns[i]) PlotLegend->AddEntry(m_sedComponentSSAFitFcns[i],sedComponentSSAFitLegendTexts[i],"L");	
	}
	for(size_t i=0;i<m_sedComponentSIFFAFitFcns.size();i++){
		if(m_sedComponentSIFFAFitFcns[i]) PlotLegend->AddEntry(m_sedComponentSIFFAFitFcns[i],sedComponentSIFFAFitLegendTexts[i],"L");	
	}
	for(size_t i=0;i<m_sedComponentSEFFAFitFcns.size();i++){
		if(m_sedComponentSEFFAFitFcns[i]) PlotLegend->AddEntry(m_sedComponentSEFFAFitFcns[i],sedComponentSEFFAFitLegendTexts[i],"L");	
	}
	PlotLegend->Draw("same");

	return Plot;

}//close DrawSED()

}//close namespace 


