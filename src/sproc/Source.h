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

		/**
		* \brief Source type enumeration
		*/
		//enum SourceType {eUnknown=0,eCompact=1,ePointLike=2,eExtended=3,eCompactPlusExtended=4};

		/**
		* \brief Source flag enumeration
		*/
		//enum SourceFlag {eReal=1,eCandidate=2,eFake=3};

		/**
		* \brief Simulated source type enumeration
		*/
		//enum SimSourceType {eUnknownSimClass=0,eRingLike=1,eBubbleLike=2,eEllipseLike=3,eDiskLike=4,eBlobLike=5};
		

	public:
		/**
		* \brief Set source type
		*/
		void SetType(SourceType choice){Type=choice;}
		/**
		* \brief Set source flag
		*/
		void SetFlag(SourceFlag choice){Flag=choice;}
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
			aNestedSource->Type= aNestedSource->Type;
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
			if(sources.empty()){
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
		std::string GetDS9RegionColor(){
			std::string colorStr= "white";
			if(Type==eExtended) colorStr= "green";
			else if(Type==eCompactPlusExtended) colorStr= "magenta";
			else if(Type==ePointLike) colorStr= "red";
			else if(Type==eCompact) colorStr= "blue";
			else colorStr= "white";
			return colorStr;
		}//close GetDS9RegionColor()

		/**
		* \brief Get DS9 region tag according to source type
		*/
		std::string GetDS9RegionTag(){
			std::string tagStr= "unknown-type";
			if(Type==eExtended) tagStr= "extended";
			else if(Type==eCompactPlusExtended) tagStr= "extended-compact";
			else if(Type==ePointLike) tagStr= "point-like";
			else if(Type==eCompact) tagStr= "compact";
			else tagStr= "unknown-type";
			return tagStr;
		}//close GetDS9RegionTag()

		//================================================
		//==         UTILS
		//================================================
		/**
		* \brief Dump source info
		*/
		void Print(){
			cout<<"*** SOURCE NO. "<<Id<<" (NAME: "<<this->GetName()<<") ***"<<endl;
			cout<<"N= "<<NPix<<", Type="<<Type<<", Pos("<<X0<<","<<Y0<<"), BoundingBox(["<<m_Xmin<<","<<m_Xmax<<"], ["<<m_Ymin<<","<<m_Ymax<<"])"<<endl;
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
		int GetFitEllipses(std::vector<TEllipse*>& fitEllipses,bool useFWHM=true,bool convertToWCS=false,WCS* wcs=0,int coordSystem=-1,int pixOffset=0,bool useWCSSimpleConversion=true);
	
		
		/**
		* \brief Get number of fit components
		*/
		int GetNFitComponents(){
			if(!m_HasFitInfo) return 0;
			return m_fitPars.GetNComponents();
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

		

	protected:
		/**
		* \brief Find source match by position
		*/
		bool FindSourceMatchByPos(std::vector<SourcePosMatchPars>& pars, long int source_index,long int nested_source_index,Source* source, float matchPosThr);
		
		/**
		* \brief Return centroid in WCS coordinates
		*/
		int GetWCSCoords(double& xwcs,double& ywcs,double x,double y,WCS* wcs=0,int coordSystem=eJ2000);



	private:
	
		/**
		* \brief Initialize class members
		*/
		void Init();

	public:

		//Source flags
		int Type;
		int Flag;
		int SimType;
		float SimMaxScale;//in arcsec

	private:
		double m_BeamFluxIntegral;

		//Is good source?
		bool m_IsGoodSource;

		//Nested source info
		int m_DepthLevel;
		bool m_HasNestedSources;
		Source* m_NestedSource;
		std::vector<Source*> m_NestedSources;	

		//True source info
		bool m_HasTrueInfo;
		double m_S_true;
		double m_X0_true;
		double m_Y0_true;

		//Fit info
		bool m_HasFitInfo;
		SourceFitPars m_fitPars;
		int m_fitStatus;

		ClassDef(Source,3)

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
#pragma link C++ enum SourceType+;
#pragma link C++ enum SourceFlag+;
#pragma link C++ enum SimSourceType+;
#endif

}//close namespace



#endif
