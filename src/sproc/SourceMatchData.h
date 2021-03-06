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
#include <SourceGroup.h>
#include <SpectralIndexData.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

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

//=======================================
//==  CLASS: SourceMatchPars
//=======================================
class SourceMatchPars : public TObject
{
	public:	
		SourceMatchPars()
			: TObject()
		{
			Init();
		}
		virtual ~SourceMatchPars(){}
		
	private:
		
		void Init()
		{
			posDist_Euclidean= -999;
			posDist_Haversine= -999;
			contourOverlapArea= -999;
			contourArea1= -999;
			contourArea2= -999;
			contourOverlapFlag= -999;
			contourOverlapAreaRatio1= -999;
			contourOverlapAreaRatio2= -999;
			fluxRelDiff= -999;
		}


	public:

		double posDist_Euclidean;
		double posDist_Haversine;
		double contourOverlapArea;
		double contourArea1;
		double contourArea2;
		int contourOverlapFlag;
		double contourOverlapAreaRatio1;
		double contourOverlapAreaRatio2;
		double fluxRelDiff;

	ClassDef(SourceMatchPars,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SourceMatchPars+;
#pragma link C++ class vector<SourceMatchPars>+;
#pragma link C++ class vector<SourceMatchPars*>+;
#endif


//=======================================
//==  CLASS: SourceMatchParsGroup
//=======================================
class SourceMatchParsGroup : public TObject
{
	public:	
		SourceMatchParsGroup()
			: TObject()
		{
			Init();
		}
		virtual ~SourceMatchParsGroup(){}
		
	public:
		int AddPars(SourceMatchPars* pars)
		{
			if(!pars) return -1;
			
			//Add pars to group
			sourceMatchPars.push_back(*pars);
			int index= sourceMatchPars.size()-1;
			
			//Update pars
			if(pars->posDist_Euclidean>posDist_Euclidean_max)	{
				posDist_Euclidean_max= pars->posDist_Euclidean;
				index_posDist_Euclidean_max= index;
			}
		  if(pars->posDist_Euclidean<posDist_Euclidean_min)	{
				posDist_Euclidean_min= pars->posDist_Euclidean;
				index_posDist_Euclidean_min= index;	
			}
			if(pars->posDist_Haversine>posDist_Haversine_max)	{
				posDist_Haversine_max= pars->posDist_Haversine;
				index_posDist_Haversine_max= index;
			}
		  if(pars->posDist_Haversine<posDist_Haversine_min)	{
				posDist_Haversine_min= pars->posDist_Haversine;
				index_posDist_Haversine_min= index;	
			}

			if(pars->contourOverlapAreaRatio1>contourOverlapAreaRatio1_max){
				contourOverlapAreaRatio1_max= pars->contourOverlapAreaRatio1;
				index_contourOverlapAreaRatio1_max= index;
			} 
			if(pars->contourOverlapAreaRatio1<contourOverlapAreaRatio1_min){
				contourOverlapAreaRatio1_min= pars->contourOverlapAreaRatio1;
				index_contourOverlapAreaRatio1_min= index;
			}
			if(pars->contourOverlapAreaRatio2>contourOverlapAreaRatio2_max){
				contourOverlapAreaRatio2_max= pars->contourOverlapAreaRatio2;
				index_contourOverlapAreaRatio2_max= index;
			}
			if(pars->contourOverlapAreaRatio2<contourOverlapAreaRatio2_min){
				contourOverlapAreaRatio2_min= pars->contourOverlapAreaRatio2;
				index_contourOverlapAreaRatio2_min= index;
			}

			return 0;
		}

	private:
		
		void Init()
		{
			sourceMatchPars.clear();
			index_posDist_Euclidean_min= -1;
			index_posDist_Euclidean_max= -1;
			posDist_Euclidean_min= 1.e+99;
			posDist_Euclidean_max= -1.e+99;	
			index_posDist_Haversine_min= -1;
			index_posDist_Haversine_max= -1;
			posDist_Haversine_min= 1.e+99;
			posDist_Haversine_max= -1.e+99;

			index_contourOverlapAreaRatio1_min= -1;
			index_contourOverlapAreaRatio1_max= -1;
			contourOverlapAreaRatio1_min= 1.e+99;
			contourOverlapAreaRatio1_max= -1.e+99;
			index_contourOverlapAreaRatio2_min= -1;
			index_contourOverlapAreaRatio2_max= -1;
			contourOverlapAreaRatio2_min= 1.e+99;
			contourOverlapAreaRatio2_max= -1.e+99;			
		}


	public:

		std::vector<SourceMatchPars> sourceMatchPars;

		int index_posDist_Euclidean_min;
		int index_posDist_Euclidean_max;
		double posDist_Euclidean_min;
		double posDist_Euclidean_max;	
		int index_posDist_Haversine_min;
		int index_posDist_Haversine_max;		
		double posDist_Haversine_min;
		double posDist_Haversine_max;
		
		int index_contourOverlapAreaRatio1_min;
		int index_contourOverlapAreaRatio1_max;
		double contourOverlapAreaRatio1_min;
		double contourOverlapAreaRatio1_max;
		int index_contourOverlapAreaRatio2_min;
		int index_contourOverlapAreaRatio2_max;
		double contourOverlapAreaRatio2_min;
		double contourOverlapAreaRatio2_max;

	ClassDef(SourceMatchParsGroup,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SourceMatchParsGroup+;
#pragma link C++ class vector<SourceMatchParsGroup>+;
#pragma link C++ class vector<SourceMatchParsGroup*>+;
#endif



//=======================================
//==  CLASS: SourceComponentMatchPars
//=======================================
class SourceComponentMatchPars : public TObject
{
	public:	
		SourceComponentMatchPars()
			: TObject()
		{
			Init();
		}
		virtual ~SourceComponentMatchPars(){}
		
	private:
		
		void Init()
		{
			posDist_Euclidean= -999;
			posDist_Haversine= -999;
			ellipseOverlapArea= -999;
			ellipseArea1= -999;
			ellipseArea2= -999;
			ellipseOverlapFlag= -999;
			ellipseOverlapAreaRatio1= -999;
			ellipseOverlapAreaRatio2= -999;
			fluxRelDiff= -999;
		}


	public:

		double posDist_Euclidean;
		double posDist_Haversine;
		double ellipseOverlapArea;
		double ellipseArea1;
		double ellipseArea2;
		int ellipseOverlapFlag;
		double ellipseOverlapAreaRatio1;
		double ellipseOverlapAreaRatio2;
		double fluxRelDiff;

	ClassDef(SourceComponentMatchPars,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SourceComponentMatchPars+;
#pragma link C++ class vector<SourceComponentMatchPars>+;
#pragma link C++ class vector<SourceComponentMatchPars*>+;
#endif


//=============================================
//==  CLASS: SourceComponentMatchParsGroup
//=============================================
class SourceComponentMatchParsGroup : public TObject
{
	public:	
		SourceComponentMatchParsGroup()
			: TObject()
		{
			Init();
		}
		virtual ~SourceComponentMatchParsGroup(){}
		
	public:
		int AddPars(SourceComponentMatchPars* pars)
		{
			if(!pars) return -1;
			
			//Add pars to group
			sourceComponentMatchPars.push_back(*pars);
			int index= sourceComponentMatchPars.size()-1;
			
			//Update pars
			if(pars->posDist_Euclidean>posDist_Euclidean_max)	{
				posDist_Euclidean_max= pars->posDist_Euclidean;
				index_posDist_Euclidean_max= index;
			}
		  if(pars->posDist_Euclidean<posDist_Euclidean_min)	{
				posDist_Euclidean_min= pars->posDist_Euclidean;
				index_posDist_Euclidean_min= index;	
			}
			if(pars->posDist_Haversine>posDist_Haversine_max)	{
				posDist_Haversine_max= pars->posDist_Haversine;
				index_posDist_Haversine_max= index;
			}
		  if(pars->posDist_Haversine<posDist_Haversine_min)	{
				posDist_Haversine_min= pars->posDist_Haversine;
				index_posDist_Haversine_min= index;	
			}
			
			if(pars->ellipseOverlapAreaRatio1>ellipseOverlapAreaRatio1_max){
				ellipseOverlapAreaRatio1_max= pars->ellipseOverlapAreaRatio1;
				index_ellipseOverlapAreaRatio1_max= index;
			} 
			if(pars->ellipseOverlapAreaRatio1<ellipseOverlapAreaRatio1_min){
				ellipseOverlapAreaRatio1_min= pars->ellipseOverlapAreaRatio1;
				index_ellipseOverlapAreaRatio1_min= index;
			}
			if(pars->ellipseOverlapAreaRatio2>ellipseOverlapAreaRatio2_max){
				ellipseOverlapAreaRatio2_max= pars->ellipseOverlapAreaRatio2;
				index_ellipseOverlapAreaRatio2_max= index;
			}
			if(pars->ellipseOverlapAreaRatio2<ellipseOverlapAreaRatio2_min){
				ellipseOverlapAreaRatio2_min= pars->ellipseOverlapAreaRatio2;
				index_ellipseOverlapAreaRatio2_min= index;
			}
			return 0;
		}

		
	private:

		void Init()
		{
			sourceComponentMatchPars.clear();
			index_posDist_Euclidean_min= -1;
			index_posDist_Euclidean_max= -1;
			posDist_Euclidean_min= 1.e+99;
			posDist_Euclidean_max= -1.e+99;	
			index_posDist_Haversine_min= -1;
			index_posDist_Haversine_max= -1;
			posDist_Haversine_min= 1.e+99;
			posDist_Haversine_max= -1.e+99;

			index_ellipseOverlapAreaRatio1_min= -1;
			index_ellipseOverlapAreaRatio1_max= -1;
			ellipseOverlapAreaRatio1_min= 1.e+99;
			ellipseOverlapAreaRatio1_max= -1.e+99;
			index_ellipseOverlapAreaRatio2_min= -1;
			index_ellipseOverlapAreaRatio2_max= -1;
			ellipseOverlapAreaRatio2_min= 1.e+99;
			ellipseOverlapAreaRatio2_max= -1.e+99;	
		}

	public:

		std::vector<SourceComponentMatchPars> sourceComponentMatchPars;

		int index_posDist_Euclidean_min;
		int index_posDist_Euclidean_max;
		double posDist_Euclidean_min;
		double posDist_Euclidean_max;	
		int index_posDist_Haversine_min;
		int index_posDist_Haversine_max;		
		double posDist_Haversine_min;
		double posDist_Haversine_max;

		int index_ellipseOverlapAreaRatio1_min;
		int index_ellipseOverlapAreaRatio1_max;
		double ellipseOverlapAreaRatio1_min;
		double ellipseOverlapAreaRatio1_max;
		int index_ellipseOverlapAreaRatio2_min;
		int index_ellipseOverlapAreaRatio2_max;
		double ellipseOverlapAreaRatio2_min;
		double ellipseOverlapAreaRatio2_max;

	ClassDef(SourceComponentMatchParsGroup,1)
};


#ifdef __MAKECINT__
#pragma link C++ class SourceComponentMatchParsGroup+;
#pragma link C++ class vector<SourceComponentMatchParsGroup>+;
#pragma link C++ class vector<SourceComponentMatchParsGroup*>+;
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


class ComponentMatchIndexGroup : public TObject
{
	public:	
		/**
		* \brief Constructor
		*/
		ComponentMatchIndexGroup()
			: TObject()
		{
			index_list.clear();
		}

		/**
		* \brief Destructor
		*/
		virtual ~ComponentMatchIndexGroup(){}

		/**
		* \brief Copy constructor
		*/
		/*
		ComponentMatchIndexGroup(const ComponentMatchIndexGroup& cmig)		
		{
			// Copy constructor
  		index_list.clear();
  		((ComponentMatchIndexGroup&)cmig).Copy(*this);
		}
		*/
		
		/**
		* \brief Assignment Operator
		*/
		/*
		ComponentMatchIndexGroup& operator=(const ComponentMatchIndexGroup& cmig)
		{
			// Operator =
  		if (this != &cmig) ((ComponentMatchIndexGroup&)cmig).Copy(*this);
  		return *this;
		}
		*/
		/**
		* \brief Copy method
		*/
		/*
		void Copy(TObject& obj) const	
		{
			(((ComponentMatchIndexGroup&)obj).index_list).clear();
			(((ComponentMatchIndexGroup&)obj).index_list)= index_list;
		}
		*/
	
	public:
		/**
		* \brief Add index
		*/
		void AddIndex(int catindex,int sindex,int cindex)
		{
			index_list.push_back(ComponentMatchIndex(catindex,sindex,cindex));
		}

		/**
		* \brief Add index
		*/
		void AddIndex(ComponentMatchIndex cmi)
		{
			index_list.push_back(cmi);
		}

		/**
		* \brief Get indexes
		*/
		std::vector<ComponentMatchIndex>& GetIndexes(){return index_list;}
		
		/**
		* \brief Get number of indexes
		*/
		int GetNIndexes(){return static_cast<int>(index_list.size());}
		

	public:
		std::vector<ComponentMatchIndex> index_list;
		
	ClassDef(ComponentMatchIndexGroup,1)

};
#ifdef __MAKECINT__
#pragma link C++ class ComponentMatchIndexGroup+;
#pragma link C++ class vector<ComponentMatchIndexGroup>+;
#pragma link C++ class vector<ComponentMatchIndexGroup*>+;
#endif


//======================================
//==      CLASS: SOURCE MATCH DATA
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
		SourceMatchData(Source* aSource,int nCatalogs,bool cloneSource=true);
		
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceMatchData();

		/**
		* \brief Copy constructor
		*/
		//SourceMatchData(const SourceMatchData& sourceMatchData);
		
		/**
		* \brief Assignment Operator
		*/
		//SourceMatchData& operator=(const SourceMatchData& sourceMatchData);
		/**
		* \brief Copy method
		*/
		//void Copy(TObject& obj) const;
	
		
	public:

		/**
		* \brief Has source match
		*/			
		bool HasSourceMatch();

		/**
		* \brief Has source component match
		*/
		bool HasSourceComponentMatch(int componentId);

		/**
		* \brief Add matched source to group
		*/
		int AddMatchedSourceToGroup(int catalogIndex,Source* aSource,bool clone=true);
		/**
		* \brief Add matched source pars
		*/
		int AddMatchedSourcePars(int catalogIndex,SourceMatchPars* smatchpars);
		
		/**
		* \brief Add matched component index info
		*/
		int AddMatchedComponentIndex(int componentIndex,int catalogIndex,int matchedSourceIndex,int matchedComponentIndex);

		/**
		* \brief Add matched source component pars
		*/
		int AddMatchedSourceComponentPars(int componentIndex,int catalogIndex,SourceComponentMatchPars* cmatchpars);
		
		
		/**
		* \brief Compute source SED and spectral index
		*/
		int ComputeSourceSEDs();
		/**
		* \brief Compute source component SED and spectral index
		*/
		int ComputeSourceComponentSEDs();
		
		/**
		* \brief Get source spectral index data
		*/
		SpectralIndexData& GetSpectralIndexData(){return m_spectralIndexData;}
		//PLSpectralIndexData& GetSpectralIndexData(){return m_spectralIndexData;}
		/**
		* \brief Get source spectral index
		*/
		double GetSpectralIndex(){return m_spectralIndexData.spectralIndex;}
		//double GetSpectralIndex(){return m_spectralIndexData.alpha;}
		/**
		* \brief Get source spectral index error
		*/
		double GetSpectralIndexError(){return m_spectralIndexData.spectralIndexErr;}
		//double GetSpectralIndexError(){return m_spectralIndexData.alphaErr;}
		/**
		* \brief Has spectral index
		*/
		double HasSpectralIndex(){return m_spectralIndexData.hasSpectralIndex;}
		//double HasSpectralIndex(){return m_spectralIndexData.hasData;}
		/**
		* \brief Is multimatch spectral index
		*/
		double IsMultiMatchSpectralIndex(){return m_spectralIndexData.isMultiSourceMatchIndex;}
		//double IsMultiMatchSpectralIndex(){return m_spectralIndexData.isMultiMatch;}

		/**
		* \brief Get component spectral index
		*/
		std::vector<SpectralIndexData>& GetSourceComponentSpectralIndexData(){return m_componentSpectralIndexData;}

		/**
		* \brief Get 3rd deg pol spectral index data
		*/
		PolSpectralIndexData& GetPolSpectralIndexData(){return m_polSpectralIndexData;}
		/**
		* \brief Get component pol spectral index data
		*/
		std::vector<PolSpectralIndexData>& GetSourceComponentPolSpectralIndexData(){return m_componentPolSpectralIndexData;}
		/**
		* \brief Has SSA spectral index
		*/
		double HasSSASpectralIndex(){return m_ssaSpectralIndexData.hasData;}
		/**
		* \brief Get synch self abs spectral index data
		*/
		SASpectralIndexData& GetSSASpectralIndexData(){return m_ssaSpectralIndexData;}
		/**
		* \brief Get component synch self abs spectral index data
		*/
		std::vector<SASpectralIndexData>& GetSourceComponentSSASpectralIndexData(){return m_componentSSASpectralIndexData;}

		/**
		* \brief Has synch internal free-free abs  spectral index
		*/
		double HasSIFFASpectralIndex(){return m_siffaSpectralIndexData.hasData;}
		/**
		* \brief Get synch internal free-free abs spectral index data
		*/
		SASpectralIndexData& GetSIFFASpectralIndexData(){return m_siffaSpectralIndexData;}
		/**
		* \brief Get component synch self abs spectral index data
		*/
		std::vector<SASpectralIndexData>& GetSourceComponentSIFFASpectralIndexData(){return m_componentSIFFASpectralIndexData;}
		/**
		* \brief Has synch external free-free abs  spectral index
		*/
		double HasSEFFASpectralIndex(){return m_seffaSpectralIndexData.hasData;}
		/**
		* \brief Get synch external free-free abs spectral index data
		*/
		SASpectralIndexData& GetSEFFASpectralIndexData(){return m_seffaSpectralIndexData;}
		/**
		* \brief Get component synch self abs spectral index data
		*/
		std::vector<SASpectralIndexData>& GetSourceComponentSEFFASpectralIndexData(){return m_componentSEFFASpectralIndexData;}

		/**
		* \brief Get source SED graph
		*/
		TGraphErrors* GetSEDGraph(){return m_sourceSED;}
		/**
		* \brief Get source component SED graph
		*/
		std::vector<TGraphErrors*>& GetSEDComponentGraphs(){return m_sourceComponentSED;}

		/**
		* \brief Fit source SED graph
		*/
		int FitSourceSED(bool useRobustFitter=false,double rob=0.75);
		/**
		* \brief Fit source SED component graphs
		*/
		int FitSourceComponentSED(bool useRobustFitter=false,double rob=0.75);

		/**
		* \brief Get source
		*/
		Source* GetSource(){return m_source;}
		/**
		* \brief Get source name
		*/
		std::string GetSourceName(){
			if(!m_source) return std::string("");			
			return m_source->GetName();	
		}
		/**
		* \brief Get source flux
		*/
		double GetSourceFlux(){
			if(!m_source) return -999;
			if(!m_source->HasFitInfo()) return -999;
			SourceFitPars fitPars= m_source->GetFitPars();
			double beamArea= m_source->GetBeamFluxIntegral();
			double fluxDensity= 0;
			if(m_source->GetFluxDensity(fluxDensity)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to get source flux density!");
				#endif
				return -999;
			}
			double flux= fluxDensity;	
			if(beamArea>0) flux/= beamArea;
			return flux;
		}

		
		/**
		* \brief Get source match names
		*/
		int GetMatchedSourceNames(std::vector<std::string>& snames,int catalogIndex);
		
		/**
		* \brief Get source match frequencies
		*/
		int GetMatchedSourceFrequency(int catalogIndex,double& freq,double& dfreq);

		/**
		* \brief Get flux and its error. If Nsources>2 return sum
		*/		
		int GetMatchedSourceFlux(int catalogIndex,double& flux,double& fluxErr,bool& summed);

		/**
		* \brief Get matched sources
		*/
		std::vector<SourceGroup*>& GetMatchedSources(){return m_matchedSources;}
		
		/**
		* \brief Get matched sources per catalog
		*/
		int GetMatchedSourcesPerCatalog(std::vector<Source*>& sources,int catalogIndex);

		/**
		* \brief Get matched source pars
		*/
		std::vector<SourceMatchParsGroup>& GetSourceMatchParsGroup(){return m_sourceMatchPars;}

		/**
		* \brief Get matched source component pars
		*/
		std::vector<std::vector<SourceComponentMatchParsGroup>>& GetSourceComponentMatchParsGroup(){return m_sourceComponentMatchPars;}

		/**
		* \brief Get matched source component indexes
		*/
		std::vector<std::vector<ComponentMatchIndexGroup>>& GetSourceComponentMatchIndexGroup(){return m_componentMatchIndexes;}

		/**
		* \brief Get number of matched catalogs per component
		*/
		int GetComponentCatalogMatchMultiplicity(int componentId);

		/**
		* \brief Draw SED
		*/
		TCanvas* DrawSED();

		/**
		* \brief Shift flux
		*/
		void ShiftFlux(bool choice,double shift=0)
		{
			m_shiftFlux= choice;
			m_fluxShift= shift;
		}

		/**
		* \brief Is flux shifted
		*/
		bool IsFluxShifted(){return m_shiftFlux;}
		
		/**
		* \brief Shift flux
		*/
		double GetFluxShift(){return m_fluxShift;}
		

	protected:
		/**
		* \brief Fit SED with power law
		*/
		int FitPLToSourceSED(SpectralIndexData& sid,TGraphErrors* sed,bool useRobustFitter=false,double rob=0.75);
		/**
		* \brief Fit SED with 3rd deg pol
		*/
		int FitPolToSourceSED(PolSpectralIndexData& sid,TGraphErrors* sed,bool useRobustFitter=false,double rob=0.75);
		/**
		* \brief Fit SED with synchrotron self-absorption model
		*/
		int FitSSAToSourceSED(SASpectralIndexData& sid,TGraphErrors* sed,bool useRobustFitter=false,double rob=0.75);
		/**
		* \brief Fit SED with synchrotron internal free-free absorption model
		*/
		int FitSIFFAToSourceSED(SASpectralIndexData& sid,TGraphErrors* sed,bool useRobustFitter=false,double rob=0.75);
		/**
		* \brief Fit SED with synchrotron external free-free absorption model
		*/
		int FitSEFFAToSourceSED(SASpectralIndexData& sid,TGraphErrors* sed,bool useRobustFitter=false,double rob=0.75);

	private:
		
		/**
		* \brief Initialize class members
		*/
		void Init();
		

	protected:
	
		// - Source  to be crossmatched
		Source* m_source;
		bool m_ownedSource;

		//- List of matched sources per catalog
		int m_nCatalogs;
		int m_nMatches;//number of catalog matches
		std::vector<int> m_nMatchesPerCatalog;

		std::vector<SourceGroup*> m_matchedSources;
		bool m_ownedMatchedSources;

		//- List of component matches
		std::vector<std::vector<ComponentMatchIndexGroup>> m_componentMatchIndexes;
		
		//- Matched source pars
		std::vector<SourceMatchParsGroup> m_sourceMatchPars;
		std::vector<std::vector<SourceComponentMatchParsGroup>> m_sourceComponentMatchPars;

		// - Spectral index (lgS vs lgNu fit)
		SpectralIndexData m_spectralIndexData;
		std::vector<SpectralIndexData> m_componentSpectralIndexData;
		TF1* m_sedPLFitFcn;
		std::vector<TF1*> m_sedComponentPLFitFcns;

		// - 3rd degree polynomial fit (lgS vs lgNu fit)
		PolSpectralIndexData m_polSpectralIndexData;
		std::vector<PolSpectralIndexData> m_componentPolSpectralIndexData;

		// - Synchtrotron self-absorption fit (lgS vs lgNu fit)
		SASpectralIndexData m_ssaSpectralIndexData;
		std::vector<SASpectralIndexData> m_componentSSASpectralIndexData;
		TF1* m_sedSSAFitFcn;
		std::vector<TF1*> m_sedComponentSSAFitFcns;
		
		// - Synchtrotron internal free-free absorption fit (lgS vs lgNu fit)
		SASpectralIndexData m_siffaSpectralIndexData;
		std::vector<SASpectralIndexData> m_componentSIFFASpectralIndexData;
		TF1* m_sedSIFFAFitFcn;
		std::vector<TF1*> m_sedComponentSIFFAFitFcns;

		// - Synchtrotron external free-free absorption fit (lgS vs lgNu fit)
		SASpectralIndexData m_seffaSpectralIndexData;
		std::vector<SASpectralIndexData> m_componentSEFFASpectralIndexData;
		TF1* m_sedSEFFAFitFcn;
		std::vector<TF1*> m_sedComponentSEFFAFitFcns;
		
		//- Graph with source SED
		TGraphErrors* m_sourceSED;
		std::vector<TGraphErrors*> m_sourceComponentSED;

		// - Shift flux of main source?
		bool m_shiftFlux;
		double m_fluxShift;

	ClassDef(SourceMatchData,2)


};//close class SourceMatchData

#ifdef __MAKECINT__
#pragma link C++ class SourceMatchData+;
#pragma link C++ class vector<SourceMatchData>+;
#pragma link C++ class vector<SourceMatchData*>+;
#endif


}//close namespace

#endif


