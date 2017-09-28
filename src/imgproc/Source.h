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
		enum SourceType {eUnknown=0,eCompact=1,ePointLike=2,eExtended=3,eCompactPlusExtended=4};

		/**
		* \brief Source flag enumeration
		*/
		enum SourceFlag {eReal=1,eCandidate=2,eFake=3};

		/**
		* \brief Simulated source type enumeration
		*/
		enum SimSourceType {eUnknownSimClass=0,eRingLike=1,eBubbleLike=2,eEllipseLike=3,eDiskLike=4,eBlobLike=5};
		

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
		* \brief Get DS9 region info
		*/
		const std::string GetDS9Region(bool dumpNestedSourceInfo=false);
		/**
		* \brief Get DS9 ellipse info
		*/
		const std::string GetDS9EllipseRegion(bool dumpNestedSourceInfo=false);

		
		//================================================
		//==         UTILS
		//================================================
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
		int MergeSource(Source* aSource,bool copyPixels=false,bool checkIfAdjacent=true,bool computeStatPars=true,bool computeMorphPars=true);

		/**
		* \brief Get number of matching pixels between this and another source
		*/
		long int GetNMatchingPixels(Source* aSource);

		/**
		* \brief Get distance in pixels between source centroids
		*/
		float GetCentroidDistance(Source* aSource);

		/**
		* \brief Fit source with a multi-component gaussian model
		*/
		int Fit(BlobPars blobPars,int nMaxComponents=3);
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
		* \brief Get fit pars
		*/
		SourceFitPars& GetFitPars(){return m_fitPars;}
	
		/**
		* \brief Get number of fit components
		*/
		int GetNFitComponents(){
			if(!m_HasFitInfo) return 0;
			return m_fitPars.GetNComponents();
		}

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

		ClassDef(Source,1)

	public:
		
};//close Source()

typedef std::vector<Source*> SourceCollection;

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
