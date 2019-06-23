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
* @file SourceSelector.cc
* @class SourceSelector
* @brief SourceSelector class
*
* Class to select sources using quality cuts provided in a cut file
* @author S. Riggi
* @date 21/06/2019
*/

#include <SourceSelector.h>
#include <Source.h>
#include <Contour.h>
#include <Cut.h>
#include <CutParser.h>
#include <CodeUtils.h>
#include <AstroUtils.h>
#include <WCSUtils.h>
#include <Consts.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TEllipse.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::SourceSelector)

namespace Caesar {

//Initialize cut fcn registry
SourceSelector::CutFcnRegistry SourceSelector::m_cutFcnRegistry = 
{ 
	{"nPixels",NPixelsCut},
	{"minBoundingBox",MinBoundingBoxCut},
	{"goodSourceFlag",GoodSourceFlagCut},
	{"circRatio",CircularityRatioCut},
	{"elongation",ElongationCut},
	{"ellipseAreaRatio",EllipseAreaRatioCut},
	{"beamAreaRatio",BeamAreaRatioCut},
	{"sourceType",SourceTypeCut},
	{"sourceFlag",SourceFlagCut},
	{"sourceSimType",SourceSimTypeCut},
	{"hasFit",HasFitCut},	
	{"fitStatus",SourceFitStatusCut},
	{"hasGoodFitChi2",HasGoodFitChi2Cut},
	{"fitQualityFlag",SourceFitQualityFlagCut},
	{"fitFlux",SourceFluxCut},
	{"fitFluxToIslandRatio",SourceFluxToIslandRatioCut},
	{"fitComponentFlux",SourceComponentFluxCut},
	{"fitComponentCentroidInsideIsland",SourceComponentCentroidInsideIslandCut},
	{"fitComponentIsolatedCentroid",SourceComponentCentroidDistanceCut},
	{"fitComponentPeakFluxToMaxRatio",SourceComponentPeakFluxCut},
	{"fitComponentPeakSignificance",SourceComponentPeakSignificanceCut},
	{"fitComponentType",SourceComponentTypeCut},
	{"fitComponentFlag",SourceComponentFlagCut}
};

SourceSelector::SourceSelector()
{

}

SourceSelector::~SourceSelector()
{

}

int SourceSelector::SelectSources(std::vector<Source*>& sources_sel,const std::vector<Source*>& sources,std::string cutFile)
{
	//Check args
	sources_sel.clear();
	if(sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Input source collection to be selected is empty, nothing to be done!");
		#endif
		return 0;
	}
	if(cutFile==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Input cut filename is empty string!");
		#endif
		return -1;
	}

	//Read cut file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading input cut file "<<cutFile<<" ...");
	#endif
	int status= CutParser::Instance().Parse(cutFile);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to parse input cut file "<<cutFile<<"!");
		#endif
		return -1;
	}

	//## Print parsed cuts
	#ifdef LOGGING_ENABLED
		INFO_LOG("Printing parsed cuts...");
	#endif
	CutParser::Instance().PrintCuts();

	typedef std::map<std::string, Cut*> CutMap;
	CutMap cuts= CutParser::Instance().GetCuts();
	
	//## Initialize cut rejection stats
	long int nSources= 0;
	long int nSources_sel= 0;
	long int nSourceComponents= 0;
	long int nSourceComponents_sel= 0;
	std::map<std::string,long int> nSourcesRejectedPerCut;
	std::map<std::string,long int> nSourceComponentsRejectedPerCut;
	for (CutMap::const_iterator it = cuts.begin(); it!=cuts.end(); it++){
		std::string cutName= it->first;
		nSourcesRejectedPerCut[cutName]= 0;
		nSourceComponentsRejectedPerCut[cutName]= 0;
	}
	
	//## Loop over source collection and apply cuts to each source
	for(size_t i=0;i<sources.size();i++)
	{
		//## Process mother source
		Source* aSource= new Source;
		*aSource= *(sources[i]);
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Selecting source no. "<<i+1<<"/"<<sources.size()<<" ...");
		#endif
		
		//- Increment source counter before selection		
		nSources++;		
		int nComponents= aSource->GetNFitComponents();
		nSourceComponents+= nComponents;

		//- Apply cuts
		bool passed= true;

		for (CutMap::const_iterator it = cuts.begin(); it!=cuts.end(); it++){
			std::string cutName= it->first;
			Cut* cut= it->second;
				
			int nComponents_before= aSource->GetNFitComponents();

			passed= ApplyCut(cutName,aSource,cut);
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source "<<aSource->GetName()<<": cut="<<cutName<<", passed? "<<passed);
			#endif

			int nComponents_after= aSource->GetNFitComponents();
			int nComponents_rejected= nComponents_before-nComponents_after;
			//nSourceComponentsRejectedPerCut[cutName]+= nComponents_rejected; 

			if(!passed) {
				nSourcesRejectedPerCut[cutName]++;
				nSourceComponentsRejectedPerCut[cutName]+= nComponents_before;
				break;
			}
			else{
				nSourceComponentsRejectedPerCut[cutName]+= nComponents_rejected; 
			}

		}//end loop cuts
		
		//## Process nested sources
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		std::vector<Source*> nestedSources_sel;

		for(size_t j=0;j<nestedSources.size();j++)
		{
			Source* aNestedSource= new Source;
			*aNestedSource= *(nestedSources[j]);

			//- Increment source counter before selection		
			nSources++;		
			int nComponents_nested= aNestedSource->GetNFitComponents();
			nSourceComponents+= nComponents_nested;
	
			//- Apply cuts
			bool passed_nested= true;

			for (CutMap::const_iterator it = cuts.begin(); it!=cuts.end(); it++){
				std::string cutName= it->first;
				Cut* cut= it->second;
			
				int nComponents_nested_before= aNestedSource->GetNFitComponents();
	
				passed_nested= ApplyCut(cutName,aNestedSource,cut);
				#ifdef LOGGING_ENABLED
					INFO_LOG("Source "<<aNestedSource->GetName()<<": cut="<<cutName<<", passed? "<<passed_nested);
				#endif

				int nComponents_nested_after= aNestedSource->GetNFitComponents();
				int nComponents_nested_rejected= nComponents_nested_before-nComponents_nested_after;
				//nSourceComponentsRejectedPerCut[cutName]+= nComponents_nested_rejected; 

				if(!passed_nested) {
					nSourcesRejectedPerCut[cutName]++;
					nSourceComponentsRejectedPerCut[cutName]+= nComponents_nested_before; 
					break;
				}
				else{
					nSourceComponentsRejectedPerCut[cutName]+= nComponents_nested_rejected; 
				}

			}//end loop cuts

			if(passed_nested){
				nestedSources_sel.push_back(aNestedSource);
			}

		}//end loop nested sources

		//Add selected sources
		if(passed)
		{
			//- Increment source counters for mother source selected
			nSources_sel++;
			nSourceComponents_sel+= aSource->GetNFitComponents();

			//- Increment source counters for nested sources 
			nSources_sel+= nestedSources_sel.size();//nested sources
			for(size_t j=0;j<nestedSources_sel.size();j++) nSourceComponents_sel+= nestedSources_sel[j]->GetNFitComponents();		

			//- Update nested sources in mother source
			if(nestedSources_sel.empty()) aSource->ClearNestedSources();
			else aSource->SetNestedSources(nestedSources_sel);

			//- Add selected source to collection
			sources_sel.push_back(aSource);
		}
		else
		{
			if(!nestedSources_sel.empty()){
				//- Increment source counters for nested sources only
				//  Add nested sources to collection as mother sources
				nSources_sel+= nestedSources_sel.size();
				for(size_t j=0;j<nestedSources_sel.size();j++) {
					nSourceComponents_sel+= nestedSources_sel[j]->GetNFitComponents();	
					sources_sel.push_back(nestedSources_sel[j]);
				}
			}
		}
		
	}//end loop sources


	//## Print cut stats
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<nSources_sel<<"/"<<nSources<<" sources selected");
		INFO_LOG("#"<<nSourceComponents_sel<<"/"<<nSourceComponents<<" source fitted components selected");
	#endif

	std::cout<<"== SOURCE CUT REJECTION STATS =="<<std::endl;
	for (CutMap::const_iterator it = cuts.begin(); it!=cuts.end(); it++){
		std::string cutName= it->first;
		std::cout.width(35); std::cout<<std::left<<cutName<<" - "<<(double)(nSourcesRejectedPerCut[cutName])/(double)(nSources)*100.<<"%"<<std::endl;
	}
	std::cout<<"================================"<<std::endl;

	std::cout<<std::endl;

	std::cout<<"== SOURCE COMPONENTS CUT REJECTION STATS =="<<std::endl;
	for (CutMap::const_iterator it = cuts.begin(); it!=cuts.end(); it++){
		std::string cutName= it->first;
		std::cout.width(35); std::cout<<std::left<<cutName<<" - "<<(double)(nSourceComponentsRejectedPerCut[cutName])/(double)(nSourceComponents)*100.<<"%"<<std::endl;
	}
	std::cout<<"================================"<<std::endl;
	


	return 0;

}//close SelectSources()


bool SourceSelector::ApplyCut(std::string cutName,Source* source,Cut* cut)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null Source ptr given to apply cut method, selection check returning false!");
		#endif
		return false;
	}

	//Check cut
	if(!cut){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null cut ptr given to apply cut method, selection check returning false!");
		#endif
		return false;
	}

	//Find cut name in registry
	CutFcnRegistry::iterator it= m_cutFcnRegistry.find(cutName);
	if(m_cutFcnRegistry.empty() || it==m_cutFcnRegistry.end()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Cut "<<cutName<<" parsed in file is not registered in cut registry, selection check returning true!");
		#endif
		return true; 
	}

	//Check ptr to cut function
	CutFcn fcn= it->second;
	if(!fcn) {
  	#ifdef LOGGING_ENABLED
			WARN_LOG("Cut "<<cutName<<" parsed in file is not registered in cut registry, selection check returning true!");
		#endif
		return true;
	}
                        
	
	//Apply cut
	bool passed= fcn(source,cut);

	return passed;

}//close ApplyCut()

//===========================================
//==         HAS FIT INFO CUT
//===========================================
bool SourceSelector::HasFitCut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	bool passed= cut->isPassed(hasFitInfo);
	return passed;

}//close HasFit()

//===========================================
//==         FIT STATUS CUT
//===========================================
bool SourceSelector::SourceFitStatusCut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;
	int fitStatus= source->GetFitStatus();
	bool passed= cut->isPassed(fitStatus);
	return passed;

}//close SourceFitStatusCut()

//===========================================
//==         FIT CHI2 CUT
//===========================================
bool SourceSelector::HasGoodFitChi2Cut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;
	SourceFitPars fitPars= source->GetFitPars();
	double chi2= fitPars.GetChi2();
	double ndf= fitPars.GetNDF();
	double redchi2= chi2/ndf;
	bool passed= cut->isPassed(redchi2);
	return passed;

}//close HasGoodFitChi2()

//===========================================
//==         FIT QUALITY FLAG CUT
//===========================================
bool SourceSelector::SourceFitQualityFlagCut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;
	int fitQuality= source->GetFitQuality();
	bool passed= cut->isPassed(fitQuality);
	return passed;

}//close SourceFitQualityFlagCut()

//===========================================
//==         FLUX CUT
//===========================================
bool SourceSelector::SourceFluxCut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;
	double fluxDensity= 0;
	int status= source->GetFluxDensity(fluxDensity);
	if(status<0) return false;
	bool passed= cut->isPassed(fluxDensity);
	return passed;

}//close SourceFluxCut()

//===========================================
//==    FLUX TO ISLAND FLUX RATIO CUT
//===========================================
bool SourceSelector::SourceFluxToIslandRatioCut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;
	double fluxDensity= 0;
	int status= source->GetFluxDensity(fluxDensity);
	if(status<0) return false;
	double S= source->GetS();
	double fluxRatio= fluxDensity/S;
	bool passed= cut->isPassed(fluxRatio);
	return passed;
	
}//close SourceFluxToIslandRatioCut()


//===========================================
//==         FIT COMPONENT FLUX CUT
//===========================================
bool SourceSelector::SourceComponentFluxCut(Source* source,Cut* cut)
{
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;
	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();
	int nComponents_sel= 0;
	std::vector<int> componentsToBeRemoved;

	for(int k=0;k<nComponents;k++)
	{
		double fluxDensity= fitPars.GetComponentFluxDensity(k);
		bool passed= cut->isPassed(fluxDensity);
		if(passed){
			nComponents_sel++;
		}
		else{
			componentsToBeRemoved.push_back(k);
		}
	}//end loop components
	
	//If no components are left, return false (not resetting fitPars here)
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentFluxCut()


//===========================================
//==   SOURCE COMPONENT CENTROID CUT
//===========================================
bool SourceSelector::SourceComponentCentroidInsideIslandCut(Source* source,Cut* cut)
{
	//Check if has fit 
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;

	//Check if has parameters & centroid
	if( (source->GetContours()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for this source, cannot perform check, returning passed!");
		#endif
		return true;
	}
	
	//Retrieve contour & fit pars
	Contour* contour= (source->GetContours())[0];
	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();
	int nComponents_sel= 0;
	std::vector<int> componentsToBeRemoved;

	//Check if component centroid is inside source island contour
	for(int k=0;k<nComponents;k++)
	{
		double x0= fitPars.GetParValue(k,"x0");	
		double y0= fitPars.GetParValue(k,"y0");
		bool isInsideContour= contour->IsPointInsideContour(x0,y0);
		bool passed= cut->isPassed(isInsideContour);
		if(passed){
			nComponents_sel++;
		}
		else{
			componentsToBeRemoved.push_back(k);
		}
	}//end loop components

	//If no components are left, return false (not resetting fitPars here)
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentCentroidInsideIslandCut()


//==============================================
//==   SOURCE COMPONENT CENTROID DISTANCE CUT
//==============================================
bool SourceSelector::SourceComponentCentroidDistanceCut(Source* source,Cut* cut)
{
	//Check if has fit 
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;

	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();

	//Check number of components
	if(nComponents<=0) return false;//cut not passed with no components (this should not occur at this stage)
	if(nComponents==1) return true;//cut always passed with just 1 one component

	//Sort components by descending peak flux
	std::vector<double> peakFluxes;
	std::vector<double> peakFluxes_sorted;
	std::vector<size_t> sorting_indexes;
	for(int k=0;k<nComponents;k++){
		double A= fitPars.GetParValue(k,"A");	
		peakFluxes.push_back(A);
	}
	CodeUtils::sort_descending(peakFluxes,peakFluxes_sorted,sorting_indexes);


	//Check if component centroid are separated by more than specified cut
	std::set<int> removedIds;

	for(size_t i=0;i<sorting_indexes.size();i++)
	{
		int componentId_i= sorting_indexes[i];
		double x0_i= fitPars.GetParValue(componentId_i,"x0");	
		double y0_i= fitPars.GetParValue(componentId_i,"y0");

		for(size_t j=i+1;j<sorting_indexes.size();j++)
		{
			int componentId_j= sorting_indexes[j];
			double x0_j= fitPars.GetParValue(componentId_j,"x0");	
			double y0_j= fitPars.GetParValue(componentId_j,"y0");
			double dx= fabs(x0_i-x0_j);
			double dy= fabs(y0_i-y0_j);
			bool passed_x= cut->isPassed(dx);
			bool passed_y= cut->isPassed(dy);
			bool passed= (passed_x && passed_y);
			if(!passed){//remove fainter component is (i,j) are too close
				removedIds.insert(componentId_j);
			}
		}//end loop components
	}//end loop components

	

	//If no components are left, return false (not resetting fitPars here)
	int nComponents_removed= static_cast<int>(removedIds.size());
	int nComponents_sel= nComponents-nComponents_removed;
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	std::vector<int> componentsToBeRemoved(removedIds.begin(),removedIds.end());
	
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentCentroidDistanceCut()

//===========================================
//==   SOURCE COMPONENT PEAK FLUX CUT
//===========================================
bool SourceSelector::SourceComponentPeakFluxCut(Source* source,Cut* cut)
{
	//Check if has fit 
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;

	double Smax= source->GetSmax();
	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();
	int nComponents_sel= 0;
	std::vector<int> componentsToBeRemoved;

	//Check if component centroid is inside source island contour
	for(int k=0;k<nComponents;k++)
	{
		double A= fitPars.GetParValue(k,"A");	
		double peakFluxRatio= A/Smax;
		bool passed= cut->isPassed(peakFluxRatio);
		if(passed){
			nComponents_sel++;
		}
		else{
			componentsToBeRemoved.push_back(k);
		}
	}//end loop components

	//If no components are left, return false (not resetting fitPars here)
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentPeakFluxCut()


//==============================================
//==   SOURCE COMPONENT PEAK SIGNIFICANCE CUT
//==============================================
bool SourceSelector::SourceComponentPeakSignificanceCut(Source* source,Cut* cut)
{
	//Check if has fit 
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;

	double nPixels= static_cast<double>(source->NPix);
	double bkgSum= source->GetBkgSum();
	double bkgRMSSum= source->GetBkgRMSSum();
	double bkgMean= bkgSum/nPixels;
	double rmsMean= bkgRMSSum/nPixels;
	if(rmsMean==0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No bkg info stored, returning passed!");
		#endif
		return true;
	}

	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();
	int nComponents_sel= 0;
	std::vector<int> componentsToBeRemoved;

	//Check if component centroid is inside source island contour
	for(int k=0;k<nComponents;k++)
	{
		double A= fitPars.GetParValue(k,"A");	
		double Z= (A-bkgMean)/rmsMean;
		bool passed= cut->isPassed(Z);
		if(passed){
			nComponents_sel++;
		}
		else{
			componentsToBeRemoved.push_back(k);
		}
	}//end loop components

	//If no components are left, return false (not resetting fitPars here)
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentPeakSignificanceCut()


//==============================================
//==   SOURCE COMPONENT TYPE CUT
//==============================================
bool SourceSelector::SourceComponentTypeCut(Source* source,Cut* cut)
{
	//Check if has fit 
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;

	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();
	int nComponents_sel= 0;
	std::vector<int> componentsToBeRemoved;

	//Check if component centroid is inside source island contour
	for(int k=0;k<nComponents;k++)
	{
		int type;
		fitPars.GetComponentType(type,k);
		bool passed= cut->isPassed(type);
		if(passed){
			nComponents_sel++;
		}
		else{
			componentsToBeRemoved.push_back(k);
		}
	}//end loop components

	//If no components are left, return false (not resetting fitPars here)
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentTypeCut()

//==============================================
//==   SOURCE COMPONENT FLAG CUT
//==============================================
bool SourceSelector::SourceComponentFlagCut(Source* source,Cut* cut)
{
	//Check if has fit 
	bool hasFitInfo= source->HasFitInfo();
	//if(!hasFitInfo) return false;
	if(!hasFitInfo) return true;

	SourceFitPars fitPars= source->GetFitPars();
	int nComponents= source->GetNFitComponents();
	int nComponents_sel= 0;
	std::vector<int> componentsToBeRemoved;

	//Check if component centroid is inside source island contour
	for(int k=0;k<nComponents;k++)
	{
		int flag;
		fitPars.GetComponentFlag(flag,k);
		bool passed= cut->isPassed(flag);
		if(passed){
			nComponents_sel++;
		}
		else{
			componentsToBeRemoved.push_back(k);
		}
	}//end loop components

	//If no components are left, return false (not resetting fitPars here)
	if(nComponents_sel<=0){
		source->SetHasFitInfo(false);
		return false;
	}	

	//Remove components not passing the cut and set update fitPars
	if(!componentsToBeRemoved.empty()){
		fitPars.RemoveComponents(componentsToBeRemoved); 
		source->SetFitPars(fitPars);
	}

	return true;

}//close SourceComponentFlagCut()


//===========================================
//==         NPIXELS CUT
//===========================================
bool SourceSelector::NPixelsCut(Source* source,Cut* cut)
{
	long int nPixels= source->NPix;
	bool passed= cut->isPassed(nPixels);
	return passed;

}//close NPixelsCut()

//===========================================
//==         MIN BOUNDING BOX CUT
//===========================================
bool SourceSelector::MinBoundingBoxCut(Source* source,Cut* cut)
{
	//Check for line-like source
	if( (source->GetContours()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for this source, cannot perform check, returning passed!");
		#endif
		return true;
	}

	double BoundingBoxMin= ((source->GetContours())[0])->BoundingBoxMin;
	bool passed= cut->isPassed(BoundingBoxMin);
	return passed;

}//close MinBoundingBoxCut()


//===========================================
//==       GOOD SOURCE FLAG CUT
//===========================================
bool SourceSelector::GoodSourceFlagCut(Source* source,Cut* cut)
{
	bool isGoodSource= source->IsGoodSource();
	bool passed= cut->isPassed(isGoodSource);
	return passed;

}//close GoodSourceFlagCut()


//===========================================
//==    CONTOUR CIRCULARITY RATIO CUT
//===========================================
bool SourceSelector::CircularityRatioCut(Source* source,Cut* cut)
{
	//Check source has parameters
	if(!source->HasParameters()) {	
		#ifdef LOGGING_ENABLED
			WARN_LOG("No parameters are available for this source (did you compute them?), returning passed!");
		#endif
		return true;
	}
	
	//Check for contours
	if( (source->GetContours()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for this source, cannot perform check, returning passed!");
		#endif
		return true;
	}

	double circRatio= ((source->GetContours())[0])->CircularityRatio;
	bool passed= cut->isPassed(circRatio);
	return passed;

}//close CircularityRatioCut()


//===========================================
//==    CONTOUR ELONGATION CUT
//===========================================
bool SourceSelector::ElongationCut(Source* source,Cut* cut)
{
	//Check source has parameters
	if(!source->HasParameters()) {	
		#ifdef LOGGING_ENABLED
			WARN_LOG("No parameters are available for this source (did you compute them?), returning passed!");
		#endif
		return true;
	}
	
	//Check for contours
	if( (source->GetContours()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for this source, cannot perform check, returning passed!");
		#endif
		return true;
	}

	double elongation= ((source->GetContours())[0])->Elongation;
	bool passed= cut->isPassed(elongation);
	return passed;

}//close ElongationCut()

//===========================================
//==    CONTOUR ELLIPSE AREA RATIO CUT
//===========================================
bool SourceSelector::EllipseAreaRatioCut(Source* source,Cut* cut)
{
	//Check source has parameters
	if(!source->HasParameters()) {	
		#ifdef LOGGING_ENABLED
			WARN_LOG("No parameters are available for source "<<source->GetName()<<" (did you compute them?), returning passed!");
		#endif
		return true;
	}
	
	//Check for contours
	if( (source->GetContours()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for source "<<source->GetName()<<", cannot perform check, returning passed!");
		#endif
		return true;
	}

	//Check for ellipse fit
	Contour* contour= (source->GetContours())[0];
	if(!contour->HasEllipseFit){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for source "<<source->GetName()<<", cannot perform check, returning passed!");
		#endif
		return true;
	}

	double ellipseAreaRatio= contour->EllipseAreaRatio;
	bool passed= cut->isPassed(ellipseAreaRatio);
	return passed;

}//close EllipseAreaRatioCut()

//===========================================
//==    BEAM AREA RATIO CUT
//===========================================
bool SourceSelector::BeamAreaRatioCut(Source* source,Cut* cut)
{
	//Check number of beams contained in source
	double nPixels= static_cast<double>(source->NPix);
	double beamArea= source->GetBeamFluxIntegral();
	if(beamArea<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No beam area info stored for source "<<source->GetName()<<", cannot perform check, returning passed!");
		#endif
		return true;
	}
	
	//Check cut
	double nBeams= nPixels/beamArea;
	bool passed= cut->isPassed(nBeams);
	return passed;

}//close BeamAreaRatioCut()


//===========================================
//==    SOURCE TYPE CUT
//===========================================
bool SourceSelector::SourceTypeCut(Source* source,Cut* cut)
{
	//Check cut
	int type= source->Type;
	bool passed= cut->isPassed(type);
	return passed;

}//close SourceTypeCut()

//===========================================
//==    SOURCE FLAG CUT
//===========================================
bool SourceSelector::SourceFlagCut(Source* source,Cut* cut)
{
	//Check cut
	int flag= source->Flag;
	bool passed= cut->isPassed(flag);
	return passed;

}//close SourceFlagCut()

//===========================================
//==    SOURCE SIM TYPE CUT
//===========================================
bool SourceSelector::SourceSimTypeCut(Source* source,Cut* cut)
{
	//Check cut
	int type= source->SimType;
	bool passed= cut->isPassed(type);
	return passed;

}//close SourceSimTypeCut()

}//close namespace

