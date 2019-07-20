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
* @file SourceMatchData.h
* @class SourceMatchData
* @brief SourceMatchData class
*
* Class representing a source match data
* @author S. Riggi
* @date 26/03/2018
*/

#ifndef _SOURCE_MATCH_DATA_h
#define _SOURCE_MATCH_DATA_h 1

#include <Source.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

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

//class Source;

//======================================
//==      CLASS: SOURCE GROUP
//======================================
class SourceGroup : public TObject 
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceGroup()
		{
			Init();		
		}
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceGroup()
		{
			Clear();
		}

	public:
		/**
		* \brief Add source to group
		*/
		int AddSource(Source* aSource){
			if(!aSource) return -1;
			m_sources.push_back(aSource);
			return 0;
		}
		/**
		* \brief Get size of group
		*/
		int GetNSources(){return static_cast<int>(m_sources.size());}
		/**
		* \brief Get source names
		*/
		std::vector<std::string> GetSourceNames()
		{
			std::vector<std::string> snames;
			for(size_t i=0;i<m_sources.size();i++){
				std::string sname= m_sources[i]->GetName();
				snames.push_back(sname);
			}
			return snames;
		}

		/**
		* \brief Get flux and its error. If Nsources>2 return sum
		*/
		int GetFlux(double& flux,double& fluxErr,bool& summed)
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
		}

		/**
		* \brief Get frequency and relative width
		*/
		int GetFrequency(double& freq,double& dfreq)
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

	private:
		/**
		* \brief Init class data
		*/
		void Init()
		{
			m_sources.clear();	
		}
		/**
		* \brief Delete class data
		*/
		void Clear()
		{
			CodeUtils::DeletePtrCollection<Source>(m_sources);
		}
		

	private:
		// - Source collection
		std::vector<Source*> m_sources;

	ClassDef(SourceGroup,1)


};//close class SourceGroup

#ifdef __MAKECINT__
#pragma link C++ class SourceGroup+;
#pragma link C++ class vector<SourceGroup>+;
#pragma link C++ class vector<SourceGroup*>+;
#endif

class ComponentMatchIndex : public TObject
{
	public:	
		ComponentMatchIndex()
			: TObject()
		{
			catalogIndex= -1;
			sourceGroupIndex= -1;
			componentIndex= -1;
		}

		ComponentMatchIndex(int catindex,int sindex,int cindex) :
			TObject(), catalogIndex(catindex), sourceGroupIndex(sindex), componentIndex(cindex) 
		{}	

		virtual ~ComponentMatchIndex(){}
	
	public:
		int catalogIndex;
		int sourceGroupIndex;
		int componentIndex;

	ClassDef(ComponentMatchIndex,1)

};

#ifdef __MAKECINT__
#pragma link C++ class ComponentMatchIndex+;
#pragma link C++ class vector<ComponentMatchIndex>+;
#pragma link C++ class vector<ComponentMatchIndex*>+;
#endif

//======================================
//==      CLASS: SOURCE MATCH GROUP
//======================================
class SourceMatchData : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceMatchData();
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceMatchData(Source* aSource,int nCatalogs);
			
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceMatchData();

		/**
		* \brief Copy constructor
		*/
		SourceMatchData(const SourceMatchData& smatch);
		
		/**
		* \brief Assignment Operator
		*/
		SourceMatchData& operator=(const SourceMatchData& smatch);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& smatch) const;
	
		
	public:

		
		/**
		* \brief Add matched source to group
		*/
		int AddMatchedSourceToGroup(int catalogIndex,Source* aSource){
			if(!aSource) return -1;
			if(catalogIndex<0 || catalogIndex>=(int)(m_matchedSources.size())) return -1;
			int status= m_matchedSources[catalogIndex]->AddSource(aSource);
			if(status<0) return -1;
			m_nMatches= 0;
			for(size_t i=0;i<m_matchedSources.size();i++){
				if(m_matchedSources[i]->GetNSources()>0) m_nMatches++;
			}
			return 0;
		}

		/**
		* \brief Add matched component index info
		*/
		int AddMatchedComponentIndex(int componentIndex,int catalogIndex,int matchedSourceIndex,int matchedComponentIndex);
		/**
		* \brief Compute source SED and spectral index
		*/
		int ComputeSourceSEDs();
		
	private:
		
		/**
		* \brief Initialize class members
		*/
		void Init();
		

	protected:
	
		// - Source  to be crossmatched
		Source* m_source;

		//- List of matched sources per catalog
		int m_nCatalogs;
		int m_nMatches;//number of catalog matches
		std::vector<SourceGroup*> m_matchedSources;

		//- List of component matches
		std::vector<std::vector<ComponentMatchIndex>> m_componentMatchIndexes; 	

		// - Spectral index (lgS vs lgNu fit)
		bool m_hasSpectralIndex;
		bool m_isMultiSourceMatchIndex;
		double m_spectralIndex;
		double m_spectralIndexErr;
		bool m_isSpectralIndexFit;
		double m_spectralFitChi2;
		double m_spectralFitNDF;
			
		//- Graph with source SED
		TGraphAsymmErrors* m_sourceSED;
		std::vector<TGraphAsymmErrors*> m_sourceComponentSED;

	ClassDef(SourceMatchData,1)


};//close class SourceMatchData

#ifdef __MAKECINT__
#pragma link C++ class SourceMatchData+;
#pragma link C++ class vector<SourceMatchData>+;
#pragma link C++ class vector<SourceMatchData*>+;
#endif


}//close namespace

#endif


