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
* @file SLICData.h
* @class SLICData
* @brief SLIC data class
*
* Superpixel data
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SLIC_DATA_h
#define _SLIC_DATA_h 1

#include <Region.h>
#include <Contour.h>
#include <Consts.h>
#include <CodeUtils.h>

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

class Image;

class SLICData : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SLICData();
		/**
		* \brief Copy constructor
		*/
		SLICData(const SLICData& data);
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SLICData();
		/**
		* \brief Assignment Operator
		*/
		SLICData& operator=(const SLICData& region);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& region) const;

	public:	

		
		/**
		* \brief Clear data
		*/
		void Clear();
		/**
		* \brief Clear images
		*/
		void ClearImages();
		/**
		* \brief Clear regions
		*/
		void ClearRegions();
		/**
		* \brief Set data (NB: Pointer ownership is taken by this class)
		*/
		int SetData(Image* img,Image* lapl_img,Image* edge_img=0);

		/** 
		\brief Generate a superpixel partition given the passed options
 		*/
		int SPGenerator(Image* img,int regionSize=10,double regParam=1,int minRegionSize=10,bool useLogScaleMapping=false,Image* edgeImg=0);
	
		/** 
		\brief Get segmented image
 		*/
		Image* GetSegmentedImage(Image* image,int selectedTag=-1,bool normalize=false,bool binarize=false);


		//=======================================
		//==       REGION METHODS
		//=======================================
		
		/**
		* \brief Get number of regions
		*/
		int GetNRegions() const {return static_cast<int>(regions.size());}
		/**
		* \brief Return regions
		*/
		std::vector<Region*> GetRegions() {return regions;}

		/**
		* \brief Set regions
		*/
		void SetRegions(std::vector<Region*>& list) {regions= list;}

		/**
		* \brief Get region
		*/
		Region* GetRegion(int index){
			if(regions.empty() || index<0 || index>=(signed)(regions.size())) return nullptr;
			return regions[index];
		}
		/**
		* \brief Add region
		*/
		int AddRegion(Region* aRegion){	
			if(!aRegion) return -1;
			regions.push_back(aRegion);
			return 0;
		}

		/**
		* \brief Get region id from index (NB: No check done)
		*/
		long int GetRegionId(long int index){
			return regions[index]->Id;
		}

		/**
		* \brief Get region size from index (NB: No check done)
		*/
		long int GetRegionSize(long int index){
			return regions[index]->NPix;
		}

		/**
		* \brief Compute region parameters
		*/
		int ComputeRegionParameters();

		/**
		* \brief Remove regions without pixels
		*/
		void RemoveEmptyRegions();

		/**
		* \brief Delete regions specified in the list (NB: No check done)
		*/
		void DeleteRegions(std::vector<size_t>const& delete_indexes){
			Caesar::CodeUtils::DeletePtrItems(regions, delete_indexes);
		}

		/**
		* \brief Get region map regionId-->index
		*/
		void GetRegionIdMap(std::map<long int,long int>& regionIdMap) const;
		
		
		/**
		* \brief Get region index map (index-->regionId)
		*/
		void GetRegionIndexMap(std::map<long int,long int>& regionIndexMap) const;



		//=======================================
		//==       PIXEL LABELS METHODS
		//=======================================
		/**
		* \brief Set pixel labels 
		*/
		void SetPixelLabels(std::vector<long int> list){pixel_labels= list;}

		/**
		* \brief Set pixel label (NB: No check done)
		*/
		int SetPixelLabel(long int gBin,long int label,bool check=true);
		int SetPixelLabel(long int ix,long int iy,long int label,bool check=true);

		/**
		* \brief Scale pixel label (NB: No check done)
		*/	
		void ScalePixelLabel(long int gBin,long int scale){
			pixel_labels[gBin]*= scale;
		}

		/**
		* \brief Get pixel label (NB: No check is done)
		*/
		long int GetPixelLabel(long int ix,long int iy){
			long int Nx= inputImg->GetNx();
			long int gBin= ix + iy*Nx;
			return pixel_labels[gBin];
		}

		/**
		* \brief Get pixel label (NB: No check is done)
		*/
		long int GetPixelLabel(long int gBin){
			return pixel_labels[gBin];
		}
		
		

		/**
		* \brief Get image curvature stats
		*/
		ImgStats* GetCurvStats(){
			if(!laplImg) return nullptr;
			return laplImg->GetPixelStats();
		}

		/**
		* \brief Get image stats
		*/
		ImgStats* GetStats(){
			if(!inputImg) return nullptr;
			return inputImg->GetPixelStats();
		}

		/**
		* \brief Get image edge stats
		*/
		ImgStats* GetEdgeStats(){
			if(!edgeImg) return nullptr;
			return edgeImg->GetPixelStats();
		}

		/**
		* \brief Get curvature pixel value (NB: No check is done!)
		*/
		double GetScurv(long int ix,long int iy){
			return laplImg->GetPixelValue(ix,iy);
		}
		/**
		* \brief Get edgeness pixel value 
		*/
		double GetSedge(long int ix,long int iy);
		/**
		* \brief Get pixel value (NB: No check is done!)
		*/
		double GetS(long int ix,long int iy){
			return inputImg->GetPixelValue(ix,iy);
		}

		/**
		* \brief Set region id (NB: No checks performed)
		*/
		void SetRegionId(int index,int id){
			regions[index]->Id= id;
		}
	
		/**
		* \brief Add pixel to region (NB: No checks performed)
		*/
		void AddPixelToRegion(int index,Pixel* pixel){
			regions[index]->AddPixel(pixel);
		}
		
	protected:
		/** 
		\brief Initialize the superpixel data structure
 		*/
		int Init(Image* img,bool useLogScaleMapping,Image* edgeImg);


	public:

		//- The image passed to SLIC generator (normalized)
		Image* inputImg;

		//- The image borders passed to SLIC generator
		Image* edgeImg;

		//- The image laplacian
		Image* laplImg; 

		//- Matrix with pixel region labels (size=image size)
		std::vector<long int> pixel_labels;

		//- The list of generated superpixels
		std::vector<Region*> regions;

	
	friend class SLIC;

		ClassDef(SLICData,1)

};//close SLICData


typedef std::map<int,std::vector<int>> SLICBoundaryPixMap;//key: neighbour id, value: list of shared pix ids
typedef SLICBoundaryPixMap::iterator SLICBoundaryPixMapIterator;
typedef std::vector< std::vector<int> > SLICConnectedRegions;


class SLICContourData : public TObject {

	public:
		/** 
		\brief SLIC contour data constructor
 		*/
		SLICContourData(){
			contour= 0;
			ResetList();
		}
	
		/** 
		\brief SLIC contour data destructor
 		*/
		virtual ~SLICContourData() { 
			ResetContour();		
			ResetList();
		}

	public:
		/** 
		\brief Reset SLIC contour
 		*/
		void ResetContour(){
			if(!contour) return;
			delete contour;
			contour= 0;	
		}

		/** 
		\brief Reset connected regions list and clear SLIC boundary data
 		*/
		void ResetList(){
			connectedRegionIds.clear();
			boundaryData.clear();
		}
	
	public:
		Contour* contour;
		std::vector<SLICBoundaryPixMap> boundaryData;
		SLICConnectedRegions connectedRegionIds; 

	ClassDef(SLICContourData,1)

};//close SLICContourData


class SLICNeighborData : public TObject {
	
	public:
		/** 
		\brief SLIC neighbor data constructor
 		*/
		SLICNeighborData(){
			Order= 0;
			Id= -1;
			Index= -1;
			Tag= -1;
			D= 0;
			D_n= 0;
			E= 0;
			E_n= 0;
			Dtot= 0;
			Dtot_n= 0;
			Dsym= 0;
		};

		/** 
		\brief SLIC neighbor data destructor
 		*/
		virtual ~SLICNeighborData(){};

	public:
		/** 
		\brief SLIC neighbor data equality operator for comparison
 		*/
		bool operator==(const SLICNeighborData& aNeighborData) const {
    	return ( (aNeighborData.Id==Id) && (aNeighborData.Index==Index) );
    }
	
		/** 
		\brief Print neighbor data info
 		*/
		void Print(){
			cout<<"** NEIGHBOR DATA #"<<Index<<" **"<<endl;
			cout<<"Id="<<Id<<", D="<<D<<" Dtot="<<Dtot<<" Dsym="<<Dsym<<endl;
			cout<<"*******************"<<endl;
		}

	public:
		int Order;
		int Id;
		int Index;
		int Tag;	
		double D;//dissimilarity Dij
		double D_n;//dissimilarity Dji
		double E;//edgeness
		double E_n;
		double Dtot;//total dissimilarity (normalized)
		double Dtot_n;
		double Dsym;//symmetric dissimilarity

	ClassDef(SLICNeighborData,1)

};//close SLICNeighborData
typedef std::vector< std::vector<SLICNeighborData> > SLICNeighbors;


class SLICNeighborCollection : public TObject {
	
	public:
		/** 
		\brief SLIC neighbor collection constructor
 		*/
		SLICNeighborCollection(){
			m_neighbors.clear();
		};

		/** 
		\brief SLIC neighbor collection destructor
 		*/
		virtual ~SLICNeighborCollection(){};

	private:
		/** 
		\brief SLIC neighbor data comparator by id 
 		*/
		struct MatchId {
 			MatchId(const int& id) : m_id(id) {}
 			bool operator()(const SLICNeighborData& obj) const {
   			return obj.Id == m_id;
 			}
 			private:
   			const int& m_id;
		};

		/** 
		\brief SLIC neighbor data comparator by index
 		*/
		struct MatchIndex {
 			MatchIndex(const int& index) : m_index(index) {}
 			bool operator()(const SLICNeighborData& obj) const {
   			return obj.Index == m_index;
 			}
 			private:
   			const int& m_index;
		};

		/** 
		\brief SLIC neighbor data comparator by dissimilarity measure 
 		*/
		static bool compareByDiss(const SLICNeighborData &a, const SLICNeighborData &b) {
   		return ( (a.D) < (b.D) );
		}

		/** 
		\brief SLIC neighbor data comparator by edgeness measure 
 		*/
		static bool compareByEdgeness(const SLICNeighborData &a, const SLICNeighborData &b) {
   		return ( (a.E) < (b.E) );
		}
		
		/** 
		\brief SLIC neighbor data comparator by total dissimilarity measure 
 		*/
		static bool compareByTotDiss(const SLICNeighborData &a, const SLICNeighborData &b) {
   		return ( (a.Dtot) < (b.Dtot) );
		}

	public: 

		/** 
		\brief Add a SLIC neighbor data to collection
 		*/
		void Add(SLICNeighborData nn){m_neighbors.push_back(nn);}

		/** 
		\brief Sort SLIC neighbor data in collection by dissimilarity measure
 		*/
		void SortByDiss(){
			std::sort(m_neighbors.begin(), m_neighbors.end(), compareByDiss);
		};

		/** 
		\brief Sort SLIC neighbor data in collection by edgeness measure
 		*/
		void SortByEdgeness(){
			std::sort(m_neighbors.begin(), m_neighbors.end(), compareByEdgeness);
		};

		/** 
		\brief Sort SLIC neighbor data in collection by total dissimilarity measure
 		*/
		void SortByTotDiss(){
			std::sort(m_neighbors.begin(), m_neighbors.end(), compareByTotDiss);
		};

		/** 
		\brief Get number of neighbor data present in collection
 		*/
		int GetN(){return (int)(m_neighbors.size());}

		/** 
		\brief Find a SLIC neighbor data in collection by id
 		*/
		int FindById(int id){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it = std::find_if(m_neighbors.begin(), m_neighbors.end(), MatchId(id));
			if (it==m_neighbors.end()) return -1;//not found in collection
			int pos = it-m_neighbors.begin();
			return pos;
		}

		/** 
		\brief Find a SLIC neighbor data in collection by index
 		*/
		int FindByIndex(int index){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it = std::find_if(m_neighbors.begin(), m_neighbors.end(), MatchIndex(index));
			if (it== m_neighbors.end()) return -1;//not found in collection
			size_t pos = it-m_neighbors.begin();
			return pos;
		}

		/** 
		\brief Find the closest neighbor data in collection wrt dissimilarity
 		*/
		int FindCloserByDiss(){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it= std::min_element(m_neighbors.begin(),m_neighbors.end(),compareByDiss);
			size_t pos = it-m_neighbors.begin();
			return pos;
		}

		/** 
		\brief Find the closest neighbor data in collection wrt total dissimilarity
 		*/
		int FindCloserByDissTot(){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it= std::min_element(m_neighbors.begin(),m_neighbors.end(),compareByTotDiss);
			size_t pos = it-m_neighbors.begin();
			return pos;
		}

		/** 
		\brief Find the closest neighbor data in collection in dissimilarity
 		*/
		std::vector<SLICNeighborData> GetNSortedByDiss(int N){
			std::vector<SLICNeighborData> topSorted;
			topSorted.clear();
			int collSize= GetN();
			if(collSize<=0 || N>collSize) return topSorted;
			
			//Partially sort Nth elements
			std::vector<SLICNeighborData> tmp;
			tmp.assign(m_neighbors.begin(),m_neighbors.end());
			partial_sort(tmp.begin(), tmp.begin()+N, tmp.end(), compareByDiss);

			//Copy first N sorted to a new collection
			topSorted.assign(m_neighbors.begin(), m_neighbors.begin()+N);
			
			return topSorted;
		}

		/** 
		\brief Returns the number of SLIC neighbor data sorted by total dissimilarity
 		*/
		std::vector<SLICNeighborData> GetNSortedByTotDiss(int N){
			std::vector<SLICNeighborData> topSorted;
			topSorted.clear();
			int collSize= GetN();
			if(collSize<=0 || N>collSize) return topSorted;
			
			//Partially sort Nth elements
			std::vector<SLICNeighborData> tmp;
			tmp.assign(m_neighbors.begin(),m_neighbors.end());
			partial_sort(tmp.begin(), tmp.begin()+N, tmp.end(), compareByTotDiss);

			//Copy first N sorted to a new collection
			topSorted.assign(tmp.begin(), tmp.begin()+N);
			
			return topSorted;
		}

		/** 
		\brief Check if neighbor data with given id is among the closest N in collection wrt dissimilarity measure
 		*/
		bool IsIdAmongNClosersByDiss(int id,int N){		
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchId(id));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}

		/** 
		\brief Check if neighbor data with given id is among the closest N in collection wrt total dissimilarity measure
 		*/
		bool IsIdAmongNClosersByTotDiss(int id,int N){		
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByTotDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchId(id));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}

		/** 
		\brief Check if neighbor data with given index is among the closest N in collection wrt dissimilarity measure
 		*/
		bool IsIndexAmongNClosersByDiss(int index,int N){
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchIndex(index));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}		

		/** 
		\brief Check if neighbor data with given index is among the closest N in collection wrt total dissimilarity measure
 		*/
		bool IsIndexAmongNClosersByTotDiss(int index,int N){
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByTotDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchIndex(index));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}		

		/** 
		\brief Print neighbor data collection info
 		*/
		void Print(){
			for(int i=0;i<GetN();i++) m_neighbors[i].Print();
		}

		/** 
		\brief Returns SLIC neighbor data collection
 		*/
		const std::vector<SLICNeighborData>& GetNeighbors() const { return m_neighbors; }

		/** 
		\brief Returns SLIC neighbor data with given index in collection
 		*/
		SLICNeighborData* GetNeighbor(int index){
			if(index<0 || index>=GetN()) return 0;
			return (m_neighbors.data()+index);
		}

		/** 
		\brief Set total dissimilarity info of neighbor data with the given index in collection
 		*/
		int SetDtot(int index,double Dtot,double Dtot_n);


	private:
		std::vector<SLICNeighborData> m_neighbors;
		
	ClassDef(SLICNeighborCollection,1)

};//close SLICNeighborCollection
typedef std::vector<SLICNeighborCollection> SLICNeighborCollections;




class SLICSimilarityData : public TObject {

	public:
	
		/** 
		\brief SLIC similarity data constructor
 		*/
		SLICSimilarityData()
		{
			AdjacencyMatrix= 0;
		}

		/** 
		\brief SLIC similarity data destructor
 		*/
		virtual ~SLICSimilarityData() 
		{ 
			//Clear matrix
			if(AdjacencyMatrix) AdjacencyMatrix->Delete();
		}

	public:

		TMatrixD* AdjacencyMatrix;
		double Dmin;
		double Dmax;
		double Dmedian;
		double Dmedianrms;
		double Emin;
		double Emax;
			
	ClassDef(SLICSimilarityData,1)
	
};//close SLICSimilarityData


#ifdef __MAKECINT__
#pragma link C++ class SLICData+;
#pragma link C++ enum SLICEdgeModel+;
#pragma link C++ class SLICContourData+;
#pragma link C++ class SLICNeighborData+;
#pragma link C++ class SLICNeighborCollection+;
#pragma link C++ class SLICSimilarityData+;
#endif

}//close namespace

#endif

