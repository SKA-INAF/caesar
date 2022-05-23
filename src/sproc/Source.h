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
* @file Source.h
* @class Source
* @brief Source class
*
* Class representing an image source
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SOURCE_h
#define _SOURCE_h 1

#include <SourceFitter.h>
#include <Blob.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <SpectralIndexData.h>
#include <AstroObject.h>

#include <TObject.h>
#include <TMatrixD.h>

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

class Contour;
class WCS;

//======================================
//==      STRUCT: SourceOverlapMatchPars
//======================================
struct SourceOverlapMatchPars {

	//Standard constructor
	SourceOverlapMatchPars(){
		ResetPars();
	}

	//Param constructor
	SourceOverlapMatchPars(long int _index,float _fraction, float _fraction_rec):
		index(_index), overlapFraction(_fraction), overlapFraction_rec(_fraction_rec)
	{}
	
	//Reset pars
	void ResetPars(){
		index=-1; 
		overlapFraction=0; 
		overlapFraction_rec=0; 
		overlappingSourceIndexes.clear();
		Sratio= 0;
		Sratio_rec= 0;
		dX= 0;
		dY= 0;
	}

	//Pars
	long int index;//index of match source in collection
	float overlapFraction;//overlap fraction with respect to true source (>0)
	float overlapFraction_rec;//overlap fraction with respect to rec source (>0)
	float Sratio;//ratio of integrated flux of overlap pixels over total flux of true source
	float Sratio_rec;//ratio of integrated flux of overlap pixels over total flux of rec source
	float dX;//difference (in pixels) between signal-weighted centroids of true and rec sources in x coordinate (rec-true)
	float dY;//difference (in pixels) between signal-weighted centroids of true and rec sources in y coordinate (rec-true)
	std::vector<long int> overlappingSourceIndexes;//list of source index in collection overlapping with this source

};//close SourceOverlapMatchPars struct


//======================================
//==      STRUCT: SourcePosMatchPars
//======================================
struct SourcePosMatchPars {

	//Standard constructor
	SourcePosMatchPars(){
		ResetPars();
	}
	//Param constructor
	SourcePosMatchPars(long int _index,float _posDiff,long int _fitComponentIndex=-1,int _nestedIndex=-1):
		index(_index), posDiff(_posDiff), fitComponentIndex(_fitComponentIndex), nestedIndex(_nestedIndex)
	{}
	
	//Reset pars
	void ResetPars(){
		index= -1; 
		fitComponentIndex= -1; 
		posDiff= 0;
		nestedIndex= -1;
	}

	//Comparison operator for sorting
	bool operator<(const SourcePosMatchPars& obj) const {
  	return posDiff < obj.posDiff;
  }

	//Pars
	long int index;//index of match source in collection
	float posDiff;//posDiff (>0)
	long int fitComponentIndex;
	long int nestedIndex;
	
};//close SourcePosMatchPars struct

//======================================
//==      CLASS: SOURCE
//======================================
class Source : public Blob {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Source();
		/** 
		\brief Parametric constructor
 		*/
		Source(std::string name);		
		/** 
		\brief Parametric constructor
 		*/
		Source(std::vector<Pixel*>const& pixels,std::string name="");			
		/**
		* \brief Copy constructor
		*/
		Source(const Source& source);
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Source();

		/**
		* \brief Assignment Operator
		*/
		Source& operator=(const Source &source);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& source) const;

		
	public:
		/**
		* \brief Set source morphology id
		*/
		void SetMorphId(SourceMorphology choice){MorphId=choice;}
		/**
		* \brief Set sourceness id
		*/
		void SetSourcenessId(Sourceness choice){SourcenessId=choice;}
		/**
		* \brief Set source sim type
		*/
		void SetSimType(SimSourceType choice){SimType=choice;}
		/**
		* \brief Set source sim max scale
		*/
		void SetSimMaxScale(float val){SimMaxScale=val;}
		/**
		* \brief Set beam flux integral
		*/
		void SetBeamFluxIntegral(double val){m_BeamFluxIntegral= val;}
		/**
		* \brief Get beam flux integral
		*/
		double GetBeamFluxIntegral(){return m_BeamFluxIntegral;}
		/**
		* \brief Is a "good" source
		*/
		bool IsGoodSource(){return m_IsGoodSource;}
		/**
		* \brief Set source as "good"
		*/
		void SetGoodSourceFlag(bool flag){m_IsGoodSource=flag;}

		/**
		* \brief Set source depth level (0=mother, 1=nested)
		*/
		void SetDepthLevel(int level){m_DepthLevel=level;}		
		/**
		* \brief Get source depth level (0=mother, 1=nested)
		*/
		int GetDepthLevel(){return m_DepthLevel;}

		
		/**
		* \brief Add nested sources
		*/
		void AddNestedSource(Source* aNestedSource){
			if(!aNestedSource) return;
			int nNestedSources= (int)m_NestedSources.size();
			int nestedId= nNestedSources+1;
			TString nestedName= Form("%s_N%d",this->GetName(),nestedId);
			aNestedSource->Id= nestedId;
			aNestedSource->MorphId= aNestedSource->MorphId;
			aNestedSource->SetName(std::string(nestedName));
			aNestedSource->m_DepthLevel= this->m_DepthLevel+1;
			m_NestedSources.push_back(aNestedSource);
			m_HasNestedSources= true;
		}	
		/**
		* \brief Has nested sources?
		*/
		bool HasNestedSources(){return (m_HasNestedSources && m_NestedSources.size()>0);}
		/**
		* \brief Set has nested sources
		*/
		void SetHasNestedSources(bool val){m_HasNestedSources=val;}

		/**
		* \brief Get nested sources
		*/
		std::vector<Source*>& GetNestedSources(){return m_NestedSources;}

		/**
		* \brief Clear nested sources
		*/
		int ClearNestedSources()
		{
			if(m_NestedSources.empty()) return 0;
			for(size_t i=0;i<m_NestedSources.size();i++){
				if(m_NestedSources[i]){
					delete m_NestedSources[i];	
					m_NestedSources[i]= 0;
				}
			}
			m_NestedSources.clear();
			m_HasNestedSources= false;
			return 0;
		}
			

		/**
		* \brief Set nested sources
		*/
		int SetNestedSources(std::vector<Source*>& sources,bool clear_existing=true)
		{
			//Check input list
			if(sources.empty()){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Given nested collection to be set is empty, will remove all existing sources!");
				#endif
			}
	
			//Release memory of existing collection?
			if(clear_existing){
				for(size_t i=0;i<m_NestedSources.size();i++){
					if(m_NestedSources[i]){
						delete m_NestedSources[i];	
						m_NestedSources[i]= 0;
					}
				}
			}//close if
			m_NestedSources.clear();

			//Add new collection
			if(!sources.empty()){
				m_NestedSources.insert(m_NestedSources.end(),sources.begin(),sources.end());
				m_HasNestedSources= true;
			}
			else{
				m_HasNestedSources= false;
			}

			return 0;

		}//close SetNestedSources()

		/**
		* \brief Get nested source number
		*/
		int GetNestedSourceNumber(){return m_NestedSources.size();}

		/**
		* \brief Get nested source
		*/
		Source* GetNestedSource(int index){
			if(index<0 || index>=(int)m_NestedSources.size() || m_NestedSources.size()==0) return 0;
			return m_NestedSources[index];
		}

		/**
		* \brief Find nested sources
		*/
		int FindNestedSources(std::vector<Source*>& nestedSources,double nestedBlobMinScale=1,double nestedBlobMaxScale=3,double nestedBlobScaleStep=1,double nestedBlobPeakZThr=5,double nestedBlobPeakZMergeThr=2.5,int minPixels=5,double nestedBlobThrFactor=0,double nestedBlobKernFactor=6,double minNestedMotherDist=2,double maxMatchingPixFraction=0.5);
	
		/**
		* \brief Are nested part of composite source?
		*/
		bool AreNestedComponentsOfCompositeSource(){return m_NestedAsCompositeSourceComponents;}
		/**
		* \brief Set has nested sources
		*/
		void SetAreNestedComponentsOfCompositeSource(bool val){m_NestedAsCompositeSourceComponents=val;}


		/**
		* \brief Draw contours
		*/
		void Draw(bool drawBoundingBox=false,bool drawFittedEllipse=false,bool drawNested=false,int lineColor=kBlack,int lineStyle=kSolid);
		/**
		* \brief Draw source
		*/
		int Draw(int pixMargin=0,ImgType imgType=eFluxMap,bool drawImage=true,bool drawContours=true,bool drawNested=true,bool drawFitComponents=true,int lineColor=kBlack,int lineStyle=kSolid,bool useWCS=false,int coordSyst=0);

		/**
		* \brief Get DS9 region info
		*/
		const std::string GetDS9Region(bool dumpNestedSourceInfo=false,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1);
		/**
		* \brief Get DS9 ellipse info
		*/
		const std::string GetDS9EllipseRegion(bool dumpNestedSourceInfo=false);
		/**
		* \brief Get DS9 fitted ellipse info
		*/
		const std::string GetDS9FittedEllipseRegion(bool useFWHM=true,bool dumpNestedSourceInfo=false,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1,bool useWCSSimpleConversion=true);

		/**
		* \brief Get DS9 region color according to source type
		*/
		std::string GetDS9RegionColor()
		{
			std::string colorStr= "white";
			if(MorphId==eExtended) colorStr= "green";
			else if(MorphId==eCompactPlusExtended) colorStr= "magenta";
			else if(MorphId==ePointLike) colorStr= "red";
			else if(MorphId==eCompact) colorStr= "blue";
			else if(MorphId==eDiffuse) colorStr= "yellow";
			else colorStr= "white";
			return colorStr;
		}//close GetDS9RegionColor()

		/**
		* \brief Get DS9 region tag according to source type
		*/
		std::string GetDS9RegionTag()
		{
			/*
			std::string tagStr= "unknown-type";
			if(MorphId==eExtended) tagStr= "extended";
			else if(MorphId==eCompactPlusExtended) tagStr= "extended-compact";
			else if(MorphId==ePointLike) tagStr= "point-like";
			else if(MorphId==eCompact) tagStr= "compact";
			else tagStr= "unknown-type";
			return tagStr;
			*/

			std::string label= GetSourceMorphLabel(MorphId);
			return label;

		}//close GetDS9RegionTag()

		//================================================
		//==         UTILS
		//================================================
		/**
		* \brief Dump source info
		*/
		void Print(){
			cout<<"*** SOURCE NO. "<<Id<<" (NAME: "<<this->GetName()<<") ***"<<endl;
			cout<<"N= "<<NPix<<", MorphId="<<MorphId<<", Pos("<<X0<<","<<Y0<<"), BoundingBox(["<<m_Xmin<<","<<m_Xmax<<"], ["<<m_Ymin<<","<<m_Ymax<<"])"<<endl;
			cout<<"S="<<m_S<<", Smin/Smax="<<m_Smin<<"/"<<m_Smax<<", Mean="<<Mean<<", RMS="<<RMS<<", Median="<<Median<<", MedianRMS="<<MedianRMS<<endl;
			cout<<"****************************"<<endl;
		}

		/**
		* \brief Is source inside given source
		*/
		bool IsInsideSource(Source* aSource){
			if(!aSource) return false;
			bool isInsideX= (m_Xmin>=aSource->m_Xmin && m_Xmax<=aSource->m_Xmax);
			bool isInsideY= (m_Ymin>=aSource->m_Ymin && m_Ymax<=aSource->m_Ymax);
			bool isInside= (isInsideX && isInsideY);			
			return isInside;
		}	

		/**
		* \brief Check if source share boundary with given box
		*/
		bool IsAtBoxEdge(float xmin,float xmax,float ymin,float ymax){
			//Check if at box edge	
			bool isAtBoxEdgeX= (this->m_Xmin==xmin || this->m_Xmax==xmax);
			bool isAtBoxEdgeY= (this->m_Ymin==ymin || this->m_Ymax==ymax);
			bool isAtBoxEdge= (isAtBoxEdgeX || isAtBoxEdgeY);
			return isAtBoxEdge;
		}
		
		/**
		* \brief Is source inside given box
		*/
		bool HasBoxOverlap(float xmin,float xmax,float ymin,float ymax){
			//Check if overlapping
			if (this->m_Xmax < xmin) return false; // A is left of B
  		if (this->m_Xmin > xmax) return false; // A is right of B
  		if (this->m_Ymax < ymin) return false; // A is above B
  		if (this->m_Ymin > ymax) return false; // A is below B
			return true;//boxes overlap
		}

		/**
		* \brief Check if this source bounding box overlaps with another given source
		*/
		bool CheckBoxOverlapping(Source*);
	
		/**
		* \brief Check if this source is adjacent to another given source
		*/
		bool IsAdjacentSource(Source* aSource);

		/**
		* \brief Merge this source with given source
		*/
		int MergeSource(Source* aSource,bool copyPixels=false,bool checkIfAdjacent=true,bool computeStatPars=true,bool computeMorphPars=true,bool sumMatchingPixels=false);

		/**
		* \brief Get collection of matching pixels between this and another source
		*/
		long int GetNMatchingPixels(std::vector<Pixel*>& matching_pixels,Source* aSource,bool sorted=false);

		/**
		* \brief Get number of matching pixels between this and another source
		*/
		long int GetNMatchingPixels(Source* aSource,bool sorted=false){
			std::vector<Pixel*> matching_pixels;
			return GetNMatchingPixels(matching_pixels,aSource,sorted);
		}

		/**
		* \brief Find source match in a collection by overlapping area
		*/
		bool FindSourceMatchByOverlapArea(SourceOverlapMatchPars& pars, const std::vector<Source*>& sources, float overlapThr);

		/**
		* \brief Find source match in a collection by position
		*/
		bool FindSourceMatchByPos(std::vector<SourcePosMatchPars>& pars, const std::vector<Source*>& sources, float posThr);

		
		/**
		* \brief Get distance in pixels between source centroids
		*/
		float GetCentroidDistance(Source* aSource);

		/**
		* \brief Fit source with a multi-component gaussian model
		*/
		int Fit(SourceFitOptions& fitOptions);
		/**
		* \brief Fit source with a multi-component gaussian model using provided start fit parameters
		*/
		int Fit(SourceFitOptions& fitOptions,SourceFitPars& initfitPars);

		/**
		* \brief Set true source info
		*/
		void SetTrueInfo(double S_true,double X0_true,double Y0_true){
			m_S_true= S_true;
			m_X0_true= X0_true;
			m_Y0_true= Y0_true;
			m_HasTrueInfo= true;	
		}

		/**
		* \brief Has true source info
		*/
		bool HasTrueInfo(){return m_HasTrueInfo;}
		/**
		* \brief Get true source flux
		*/
		double GetTrueFlux(){return m_S_true;}
		/**
		* \brief Get true source position
		*/
		void GetTruePos(double& x,double& y){x= m_X0_true; y=m_Y0_true;}

		/**
		* \brief Has fit info
		*/
		bool HasFitInfo(){return m_HasFitInfo;}

		/**
		* \brief Set Has fit info (for serialization scopes)
		*/
		void SetHasFitInfo(bool flag){m_HasFitInfo= flag;}

		/**
		* \brief Get fit pars
		*/
		SourceFitPars& GetFitPars(){return m_fitPars;}

		/**
		* \brief Set fit pars (for serialization scopes)
		*/
		void SetFitPars(SourceFitPars& fitPars){
			m_fitPars= fitPars;
			m_HasFitInfo= true;
		}

		/**
		* \brief Set fit component flag
		*/
		int SetFitComponentSourcenessId(int componentId,int flag)
		{
			if(!m_HasFitInfo) return -1;
			return m_fitPars.SetComponentSourcenessId(componentId, flag);
		}

		/**
		* \brief Get fit component flag
		*/
		int GetFitComponentSourcenessId(int& flag, int componentId)
		{
			if(!m_HasFitInfo) return -1;
			return m_fitPars.GetComponentSourcenessId(flag, componentId);
		}

		/**
		* \brief Get fit status
		*/
		int GetFitStatus(){return m_fitStatus;}

		/**
		* \brief Set fit status (for serialization scopes)
		*/
		void SetFitStatus(int fitStatus){m_fitStatus=fitStatus;}

		/**
		* \brief Get integrated flux density
		*/
		int GetFluxDensity(double& fluxDensity){
			fluxDensity= 0;
			if(!m_HasFitInfo) return -1;
			fluxDensity= m_fitPars.GetFluxDensity();
			return 0;
		}

		/**
		* \brief Get integrated flux density error
		*/
		int GetFluxDensityErr(double& fluxDensityErr){
			fluxDensityErr= 0;
			if(!m_HasFitInfo) return -1;
			fluxDensityErr= m_fitPars.GetFluxDensityErr();
			return 0;
		}

		/**
		* \brief Get fit quality flag
		*/
		int GetFitQuality(){
			if(!m_HasFitInfo) return -1;
			return m_fitPars.GetFitQuality();
		}

		/**
		* \brief Get integrated flux density error on components according to Condon (1997) formula 14
		*/
		int GetCondonComponentFluxDensityErr(std::vector<double>& fluxDensityErrList){
			fluxDensityErrList.clear();
			if(!m_HasFitInfo) return -1;
			if(!m_imgMetaData) return -1;
			int nComponents= m_fitPars.GetNComponents();
			double rmsAvg= m_bkgRMSSum/NPix;
			double dx= fabs(m_imgMetaData->dX);
 			double dy= fabs(m_imgMetaData->dY);
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("rmsAvg="<<rmsAvg<<", dx="<<dx<<", dy="<<dy);
			#endif
 
			for(int k=0;k<nComponents;k++){
				double A= m_fitPars.GetParValue(k,"A");
				double sigmaX= m_fitPars.GetParValue(k,"sigmaX");
				double sigmaY= m_fitPars.GetParValue(k,"sigmaY");
				double SNR= A/rmsAvg;
				double I= m_fitPars.GetComponentFluxDensity(k);
				double IRelErr= sqrt(2./(TMath::Pi()*sigmaX*sigmaY))/SNR;
				//double IVar= 2.*I*I*dx*dy*rmsAvg*rmsAvg/(TMath::Pi()*sigmaX*sigmaY*A*A);
				double IErr= IRelErr*I;
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("A="<<A<<", sigmaX="<<sigmaX<<", sigmaY="<<sigmaY<<", SNR="<<SNR<<", I="<<I<<", IRelErr="<<IRelErr<<", IErr="<<IErr);
				#endif
				fluxDensityErrList.push_back(IErr);
			}//end loop components
			return 0;
		}
	
		/**
		* \brief Get fit ellipses
		*/
		int GetFitEllipses(std::vector<TEllipse*>& fitEllipses,bool useFWHM=true,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1,double pixOffset=0,bool useWCSSimpleConversion=true);

		
		/**
		* \brief Get number of fit components
		*/
		int GetNFitComponents(){
			if(!m_HasFitInfo) return 0;
			return m_fitPars.GetNComponents();
		}

		/**
		* \brief Get number of selected fit components
		*/
		int GetNSelFitComponents(){
			if(!m_HasFitInfo) return 0;
			return m_fitPars.GetNSelComponents();
		}

		/**
		* \brief Find component peaks
		*/
		int FindComponentPeaks(std::vector<ImgPeak>& peaks,double peakZThr=0,int maxPeaks=-1,int peakShiftTolerance=2,std::vector<int> kernels= {3,5,7},int peakKernelMultiplicityThr=1,bool invertSearch=false);

		/**
		* \brief Find blended source components
		*/
		int FindBlendedComponents(std::vector<Source*>& deblendedComponents,std::vector<ImgPeak>& deblendedPeaks,double peakZThr=0,int maxPeaks=-1,double sigmaMin=3,double sigmaMax=3,double sigmaStep=1,int minBlobSize=5,double thrFactor=0,int kernelFactor=1,int pixMargin=10);

		/**
		* \brief Return source position in WCS coordinates
		*/
		int GetWCSPos(double& xwcs,double& ywcs,WCS* wcs=0,int coordSystem=eJ2000)
		{
			return GetWCSCoords(xwcs,ywcs,X0,Y0,wcs,coordSystem);
		}

		/**
		* \brief Return signal-weighted source position in WCS coordinates
		*/
		int GetWCSWeightedPos(double& xwcs,double& ywcs,WCS* wcs=0,int coordSystem=eJ2000)
		{
			return GetWCSCoords(xwcs,ywcs,m_Sx,m_Sy,wcs,coordSystem);
		}

		/**
		* \brief Return bounding box in WCS coordinates
		*/
		int GetWCSSourceRange(double& xmin_wcs,double& xmax_wcs,double& ymin_wcs,double& ymax_wcs,WCS* wcs=0,int coordSystem=eJ2000)
		{
			xmin_wcs= -999;
			xmax_wcs= -999;
			ymin_wcs= -999;
			ymax_wcs= -999;
			int status= GetWCSCoords(xmin_wcs,ymin_wcs,m_Xmin,m_Ymin,wcs,coordSystem);
			if(status<0) return -1;
			status= GetWCSCoords(xmax_wcs,ymax_wcs,m_Xmax,m_Ymax,wcs,coordSystem);
			if(status<0) return -1;	
			return 0;
		}
	
		/**
		* \brief Return source name following IAU convention
		*/
		std::string GetIAUName(bool useWeightedPos=false,WCS* wcs=0,int coordSystem=eJ2000);
		
		/**
		* \brief Return spectral axis info
		*/
		int GetSpectralAxisInfo(double& val,double& dval,std::string& units);

		/**
		* \brief Has spectral index data
		*/
		bool HasSpectralIndexData(){return m_hasSpectralIndexData;}
		/**
		* \brief Set has spectral index data flag
		*/
		void SetHasSpectralIndexData(bool val){m_hasSpectralIndexData= val;}
		/**
		* \brief Get spectral index data
		*/
		SpectralIndexData& GetSpectralIndexData(){return m_spectralIndexData;}

		/**
		* \brief Set spectral index data 
		*/
		void SetSpectralIndexData(SpectralIndexData& data){
			m_spectralIndexData= data;
			if(data.hasSpectralIndex) m_hasSpectralIndexData= true;
		}
		/**
		* \brief Has component spectral index data
		*/
		bool HasComponentSpectralIndexData(){return m_hasComponentSpectralIndexData;}
		/**
		* \brief Set has spectral index data flag
		*/
		void SetHasComponentSpectralIndexData(bool val){m_hasComponentSpectralIndexData= val;}
		/**
		* \brief Get component spectral index data
		*/
		std::vector<SpectralIndexData>& GetComponentSpectralIndexData(){return m_componentSpectralIndexData;}
		/**
		* \brief Set component spectral index data 
		*/
		void SetComponentSpectralIndexData(std::vector<SpectralIndexData>& data){
			m_componentSpectralIndexData= data;
			if(!data.empty()) m_hasComponentSpectralIndexData= true;
		}
		/**
		* \brief Has astro object data
		*/
		bool HasAstroObjects(){return m_hasAstroObjectData;}
		/**
		* \brief Has astro object data
		*/
		void SetHasAstroObjects(bool choice){m_hasAstroObjectData= choice;}
		/**
		* \brief Get astro objects data
		*/
		std::vector<AstroObject>& GetAstroObjects(){return m_astroObjects;}
		/**
		* \brief Set astro objects data
		*/
		void SetAstroObjects(std::vector<AstroObject>& data){
			m_astroObjects= data;
			if(!data.empty()){
				m_hasAstroObjectData= true;
				ComputeObjClassId();
			}
		}
		/**
		* \brief Add astro objects data
		*/
		int AddAstroObject(AstroObject& astroObject);

		/**
		* \brief Has component astro object data
		*/
		bool HasComponentAstroObjects(){return m_hasComponentAstroObjectData;}
		/**
		* \brief Has astro object data
		*/
		void SetHasComponentAstroObjects(bool choice){m_hasComponentAstroObjectData= choice;}
		/**
		* \brief Get astro objects data
		*/
		std::vector<std::vector<AstroObject>>& GetComponentAstroObjects(){return m_componentAstroObjects;}
		/**
		* \brief Set astro objects data
		*/
		void SetComponentAstroObjects(std::vector<std::vector<AstroObject>>& data){
			m_componentAstroObjects= data;
			m_hasComponentAstroObjectData= true;
			ComputeComponentObjClassId();
		}
		/**
		* \brief Add component astro objects data
		*/
		int AddComponentAstroObject(int componentIndex,AstroObject& astroObject);

		/**
		* \brief Compute object class id from astro object data (if available)
		*/
		int ComputeObjClassId();
		/**
		* \brief Compute component object class id from astro object data (if available)
		*/
		int ComputeComponentObjClassId();
		/**
		* \brief Compute if source island is resolved according to XXL survey criterion
		*/
		bool IsResolved_XXLSurveyMethod(double p0=1.08, double p1=2.03)
		{	
			double beamArea= this->GetBeamFluxIntegral();	
			double S= m_S/beamArea;
			double Speak= m_Smax;
			double bkgRMSSum= m_bkgRMSSum;
			double nPix= NPix;
			double avgBkgRMS= bkgRMSSum/nPix;
			if(m_HasFitInfo){
				double fluxDensity= m_fitPars.GetFluxDensity();
				S= fluxDensity/beamArea;
			}
			double snr= Speak/avgBkgRMS;
			double fluxToPeakRatio= S/Speak;
			double sourceResolveThr= p0 + p1/snr;
			bool resolved= (fluxToPeakRatio>sourceResolveThr);

			return resolved;
		}

		/**
		* \brief Compute if source fit component is resolved according to XXL survey criterion
		*/
		bool IsComponentResolved_XXLSurveyMethod(int componentId, double p0=1.08, double p1=2.03)
		{	
			//Check if has fit info
			if(!m_HasFitInfo) return false;

			//Compute condition
			double beamArea= this->GetBeamFluxIntegral();	
			double Speak_comp= m_fitPars.GetComponentPeakFlux(componentId);
			double fluxDensity_comp= m_fitPars.GetComponentFluxDensity(componentId);
			double S_comp= fluxDensity_comp/beamArea;
			double bkgRMSSum= m_bkgRMSSum;
			double nPix= NPix;
			double avgBkgRMS= bkgRMSSum/nPix;
			double snr_comp= Speak_comp/avgBkgRMS;
			double fluxToPeakRatio_comp= S_comp/Speak_comp;
			double sourceResolveThr_comp= p0 + p1/snr_comp;
			bool resolved_comp= (fluxToPeakRatio_comp>sourceResolveThr_comp);
		
			return resolved_comp;

		}

		/**
		* \brief Get user tags
		*/
		std::vector<std::string>& GetTags(){return m_tags;}
		/**
		* \brief Set tags
		*/
		void SetTags(std::vector<std::string>& data){
			m_tags= data;
		}

	protected:
		/**
		* \brief Find source match by position
		*/
		bool FindSourceMatchByPos(std::vector<SourcePosMatchPars>& pars, long int source_index,long int nested_source_index,Source* source, float matchPosThr);
		
		/**
		* \brief Return centroid in WCS coordinates
		*/
		int GetWCSCoords(double& xwcs,double& ywcs,double x,double y,WCS* wcs=0,int coordSystem=eJ2000);

		/**
		* \brief Compute object class id from astro object data (if available)
		*/
		int ComputeObjClassId(int& id,int& subid,std::string& id_str,bool& confirmed,std::vector<AstroObject>& data);
		

	private:
	
		/**
		* \brief Initialize class members
		*/
		void Init();

	public:

		//Source flags
		int MorphId;
		int SourcenessId;
		float SourcenessScore;
 
		int SimType;
		float SimMaxScale;//in arcsec
		
		//Object type
		int ObjLocationId;
		std::string ObjClassStrId;
		int ObjClassId;
		int ObjClassSubId;
		bool ObjConfirmed;
		float ObjClassScore;

		//Component object type
		std::vector<int> componentObjLocationIds;
		std::vector<std::string> componentObjClassStrIds;
		std::vector<int> componentObjClassIds;
		std::vector<int> componentObjClassSubIds;
		std::vector<bool> componentObjConfirmed;

		

	private:
		double m_BeamFluxIntegral;

		//Is good source?
		bool m_IsGoodSource;

		//Nested source info
		int m_DepthLevel;
		bool m_HasNestedSources;
		//Source* m_NestedSource;
		std::vector<Source*> m_NestedSources;
		bool m_NestedAsCompositeSourceComponents;

		//True source info
		bool m_HasTrueInfo;
		double m_S_true;
		double m_X0_true;
		double m_Y0_true;

		//Fit info
		bool m_HasFitInfo;
		SourceFitPars m_fitPars;
		int m_fitStatus;

		//Spectral index data	
		bool m_hasSpectralIndexData;
		SpectralIndexData m_spectralIndexData;
		bool m_hasComponentSpectralIndexData;
		std::vector<SpectralIndexData> m_componentSpectralIndexData;

		//Astro object cross-match data
		bool m_hasAstroObjectData;
		std::vector<AstroObject> m_astroObjects;
		bool m_hasComponentAstroObjectData;
		std::vector<std::vector<AstroObject>> m_componentAstroObjects;

		//User tags
		std::vector<std::string> m_tags;
		
		ClassDef(Source,11)

	public:
		
};//close Source()

typedef std::vector<Source*> SourceCollection;

struct SourceCompareByPeakFlux {
	bool operator()(const Source* lhs, const Source* rhs) const { 
		return lhs->GetSmax() < rhs->GetSmax();
	}
};//close SourceCompareByPeakFlux()

struct SourceCompareByLargerPeakFlux {
	bool operator()(const Source* lhs, const Source* rhs) const { 
		return lhs->GetSmax() > rhs->GetSmax();
	}
};//close SourceCompareByLargerPeakFlux()


#ifdef __MAKECINT__
#pragma link C++ class Source+;
#pragma link C++ class vector<Source>+;
#pragma link C++ class vector<Source*>+;
#pragma link C++ enum SourceMorphology+;
#pragma link C++ enum Sourceness+;
#pragma link C++ enum SimSourceType+;
#endif

}//close namespace



#endif
