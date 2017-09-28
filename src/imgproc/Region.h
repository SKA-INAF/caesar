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
* @file Region.h
* @class Region
* @brief Region data class
*
* Superpixel Region data
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _REGION_h
#define _REGION_h 1

#include <Blob.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>

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

using namespace std;


namespace Caesar {


struct DistPars {
	DistPars(){
		dist2= 0;
		dist2_curv= 0;
		dist2_spatial= 0;
	}
	double dist2;
	double dist2_curv;
	double dist2_spatial;
};


class Region : public Blob {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Region();
		/** 
		\brief Parametric constructor
 		*/
		//Region(ImgRange img_range,std::string name="");		
		Region(std::string name);		
		/** 
		\brief Parametric constructor
 		*/
		//Region(std::vector<Pixel*>const& pixels,ImgRange img_range,std::string name="");
		Region(std::vector<Pixel*>const& pixels,std::string name="");
		/**
		* \brief Copy constructor
		*/
		Region(const Region& region);
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Region();
		/**
		* \brief Assignment Operator
		*/
		Region& operator=(const Region& region);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& region) const;

	public:
		/**
		* \brief Region tag enumeration
		*/
		enum RegionTag {eBkgTag=0,eSignalTag=1,eUntagged=2};
	
		/**
		* \brief Region pars
		*/
		struct RegionPars {
			TVectorD* pars;
			TVectorD* robustPars;
			TVectorD* spatialPars;
			RegionPars() {
				pars= 0;
				robustPars= 0;
				spatialPars= 0;
			}
			~RegionPars() { 
				if(pars) pars->Delete();
				if(robustPars) robustPars->Delete();
				if(spatialPars) spatialPars->Delete();	
			}
		};//close RegionPars()

	public:

		/**
		* \brief Get region parameters
		*/
		Region::RegionPars* GetParams(bool includeCurvPar=true);

		/**
		* \brief Get distance squared between this and given region
		*/	
		int GetDistance(DistPars& distPars,Region* aRegion,bool useRobustParams=false);
	
		/**
		* \brief Get asymmetric color & space symmetric distance between this and given region
		*/
		int GetAsymmDistance(DistPars& distPars,DistPars& distPars_neighbor,Region* aRegion,bool useRobustParams=false);


		/**
		* \brief Get color & space symmetric distance among two regions
		*/
		//int GetDistance(double& dist_color,double& dist_space,Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addCurvDist=true);

		
		/**
		* \brief Get color & space asymmetric distance among two regions
		*/
		//int GetAsymmDistance(double& dist,double& dist_neighbor,Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addSpatialDist=false,bool addCurvDist=true);

		/**
		* \brief Merge a region to this
		*/
		int AddRegion(Region* aRegion,bool addPixels=true,bool copyPixels=false);

		/**
		* \brief Add a sub-region id to list
		*/
		int AddSubRegionId(long int id){
			//Search if a sub region with same id was already added
			int pos= -1;
			if(Caesar::CodeUtils::FindItem(m_SubRegionIds,id,pos)){
				WARN_LOG("Sub-region id "<<id<<" was already added as sub-region for region id="<<Id<<" (CHECK!!!)");
				return -1;
			}
			m_SubRegionIds.push_back(id);
			return 0;
		}

		/**
		* \brief Get number of sub-regions present in list
		*/
		int GetNSubRegions(){return (int)(m_SubRegionIds.size());}

		/**
		* \brief Get the list of sub-region ids
		*/
		const std::vector<long int>& GetSubRegionIds() const {return m_SubRegionIds;}

		/**
		* \brief Get sub-region id at given index
		*/
		long int GetSubRegionId(int index){
			size_t nSubRegions= m_SubRegionIds.size();
			if(nSubRegions<=0 || index<0 || index>=(signed)(nSubRegions) ) return -1;
			return m_SubRegionIds[index];
		}

	private:
		/**
		* \brief Get color distance between this and a given region
		*/
		double GetColorDistanceSqr(Region* aRegion,bool useRobustParams=false);

		/**
		* \brief Get color curvature distance between this and a given region
		*/
		double GetColorCurvDistanceSqr(Region* aRegion);

		/**
		* \brief Get spatial distance between this and a given region centroid
		*/
		double GetSpatialDistanceSqr(Region* aRegion);


	public:
		int Tag;

	protected:
		std::vector<long int> m_SubRegionIds;

	ClassDef(Region,1)

};//close Region()


class RegionCollection : public TObject {
	
	public:

		/**
		* \brief Constructor
		*/
		RegionCollection(){
			regions.clear();
		};

		/**
		* \brief Destructor
		*/
		virtual ~RegionCollection(){
			for(unsigned int i=0;i<regions.size();i++){
				if(regions[i]) {
					delete regions[i];	
					regions[i]= 0;
				}
			}
			regions.clear();
		};

	private:

		/**
		* \brief Matcher predicated function to find regions in collection by id
		*/
		struct MatchId {
 			MatchId(const int& id) : m_id(id) {}
 			bool operator()(const Region* obj) const {
   			return obj->Id == m_id;
 			}
 			private:
   			const int& m_id;
		};
		

	public: 

		/**
		* \brief Add a region to collection
		*/
		void Add(Region* aRegion){regions.push_back(aRegion);}

		/**
		* \brief Get number of regions present in collection
		*/
		int GetN(){return static_cast<int>(regions.size());}

		/**
		* \brief Get region by id
		*/
		Region* FindRegionById(int id){
			int index = FindRegion(id);
			if(index<0) return 0;
			return regions[index];
		}
	
		/**
		* \brief Find region (return index with given id)
		*/
		int FindRegion(int id){
			if(GetN()<=0) return -1;
			std::vector<Region*>::iterator it = std::find_if(regions.begin(), regions.end(), MatchId(id));
			if (it==regions.end()) return 0;//not found in collection
			int index = it-regions.begin();
			return index;
		}
	
		/**
		* \brief Get region with given index
		*/
		Region* GetRegion(int index){
			if(index<0 || index>=GetN()) return 0;
			return regions[index];
		}

		/**
		* \brief Get region map regionId-->index
		*/
		std::map<int,int> GetRegionIdMap() const {
			std::map<int,int> regionIdMap;
			for(unsigned int k=0;k<regions.size();k++){
				int regionId= regions[k]->Id;
				regionIdMap.insert( std::pair<int,int>(regionId,k) );
			}
			return regionIdMap;
		}//close GetRegionMapping()

		/**
		* \brief Get region index map (index-->regionId)
		*/
		std::map<int,int> GetRegionIndexMap() const {
			std::map<int,int> regionIndexMap;
			for(unsigned int k=0;k<regions.size();k++){
				int regionId= regions[k]->Id;
				regionIndexMap.insert( std::pair<int,int>(k,regionId) );
			}
			return regionIndexMap;
		}//close GetRegionMapping()
	

	public:
		std::vector<Region*> regions;
		

	ClassDef(RegionCollection,1)

};//close RegionCollection


#ifdef __MAKECINT__
#pragma link C++ class Region+;
#pragma link C++ class RegionCollection+;
#endif

}//close namespace

#endif

