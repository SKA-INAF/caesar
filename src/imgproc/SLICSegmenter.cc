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
* @file SLICSegmenter.cc
* @class SLICSegmenter
* @brief SLICSegmenter
*
* In this class, an over-segmentation is created of an image, provided by the
* step-size (distance between initial cluster locations) and the colour
* distance parameter.
* @author S. Riggi
* @date 15/06/2015
*/

#include <SLICSegmenter.h>
#include <Region.h>
#include <SLICData.h>
#include <SLIC.h>
#include <StatsUtils.h>
#include <CodeUtils.h>
#include <Image.h>

#include <TObject.h>
#include <TMatrixD.h>

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


ClassImp(Caesar::SLICSegmenter)


namespace Caesar {


SLICSegmenter::SLICSegmenter() {
	
}//close constructor

SLICSegmenter::~SLICSegmenter(){
	
}//close destructor


int SLICSegmenter::CheckData(SLICData const& slicData){

	//Check given slic input data
	if(!slicData.inputImg){
		ERROR_LOG("Null ptr to input image in slicData!");
		return -1;
	}
	if(!slicData.edgeImg){
		ERROR_LOG("Null ptr to edge image in slicData!");
		return -1;
	}	
	//if(!slicData.pixelLabels){
	//	ERROR_LOG("Null ptr to pixel labels in slicData!");
	//	return -1;
	//}

	int nRegions= slicData.GetNRegions();
	if(nRegions<=1){
		WARN_LOG("No (or too few) regions (NR="<<nRegions<<" in slicData!");
		return -1;
	}

	return 0;

}//close Init()


int SLICSegmenter::FindSegmentation(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,int minMergedSP,double SPMergingRatio, double SPMergingMaxDissRatio,double SPMergingMaxDissRatio_2ndNeighbor,double SPMergingDissThreshold,bool SPMergingIncludeSpatialPars,bool SPMergingUseRobustPars,bool SPMergingUseCurvDist){

	//## Check input data
	INFO_LOG("Checking initial slic data given...");
	if(CheckData(slicData)<0){
		ERROR_LOG("Invalid slicData given (missing edge image information or empty region list?)!");
		return -1;
	}
	
	//## Create mapping between region id and index
	std::map<long int,long int> regionIdMap_initialSegm;
	slicData.GetRegionIdMap(regionIdMap_initialSegm);
	
	//## Copy given slic data
	SLICData* slicData_segmMultiStep= new SLICData;//for 1st stage segmentation
	*slicData_segmMultiStep= slicData;
	
	segmSlicData.Clear();//for final hierarchical merging
	segmSlicData= slicData;

	//----DEBUG----
	INFO_LOG("slic data edge image stats");
	(slicData.edgeImg)->PrintStats();

	INFO_LOG("slic data edge image copy stats");
	(segmSlicData.edgeImg)->PrintStats();
	//-------------
	
	//## Adaptively merge superpixels using max similarity measure
	INFO_LOG("Adaptively merge superpixels using max similarity measure...");
	int status= MultiStepSPMerger(
		slicData, *slicData_segmMultiStep,
		SPMergingRegularization, use2ndNeighborsInSPMerging,
		SPMergingIncludeSpatialPars, SPMergingUseRobustPars, SPMergingUseCurvDist
	);

	if(status<0){
		ERROR_LOG("Multi-step superpixel merger failed!");
		delete slicData_segmMultiStep;
		slicData_segmMultiStep= 0;
		return -1;
	}

	//## Check first stage segmentation
	int nSig= 0;
	int nBkg= 0;
	int nUntagged= 0;
	INFO_LOG("Counting tagged regions after adaptive merging stage...");
	SLIC::CountTaggedRegions(slicData_segmMultiStep->regions,nSig,nBkg,nUntagged);
	if(nSig<=0 || slicData_segmMultiStep->GetNRegions()<=1){
		WARN_LOG("No signal regions (or no regions) found after adaptive merging, stop processing.");
		delete slicData_segmMultiStep;
		slicData_segmMultiStep= 0;
		return 0;
	}


	//## Create a copy of original regions and change tag according to segmented regions
	INFO_LOG("Create a copy of initial segmented regions (N="<<segmSlicData.GetNRegions()<<") and change tag according to segmented regions (N="<<slicData_segmMultiStep->GetNRegions()<<") ...");
	for(int i=0;i<segmSlicData.GetNRegions();i++) {
		(segmSlicData.regions)[i]->Tag= Region::eBkgTag;
	}//end loop original regions

	for(int i=0;i<slicData_segmMultiStep->GetNRegions();i++) {
		Region* region_postSegm= slicData_segmMultiStep->GetRegion(i);
		if(!region_postSegm){
			ERROR_LOG("Failed to get region no. "<<i<<" from slic data (this should not occur!)");
			continue;
		}
		long int regionId= region_postSegm->Id;
		int regionTag= region_postSegm->Tag;	
		//int nSubRegions= (slicData_segmMultiStep->regions)[i]->GetNSubRegions();
		int nSubRegions= region_postSegm->GetNSubRegions();
		long int regionIndex_initialSegm= regionIdMap_initialSegm[regionId];//index of region id in initial segmentation region list

		//Get access to region in initial segmentation and set the new tag 
		Region* region_initialSegm= segmSlicData.GetRegion(regionIndex_initialSegm);
		if(!region_initialSegm){
			ERROR_LOG("Failed to get region with index "<<regionIndex_initialSegm<<" from initial slic data (this should not occur!)");
			continue;
		}
		region_initialSegm->Tag= regionTag;
		
		//Get access to merged regions in initial segmentation and set the new tag
		INFO_LOG("Region no. "<<i<<" (id="<<regionId<<", tag="<<regionTag<<" initial_index="<<regionIndex_initialSegm<<") has "<<nSubRegions<<" sub-regions...");
		for(int j=0;j<nSubRegions;j++){
			//long int subregionId= (slicData_segmMultiStep->regions)[i]->GetSubRegionId(j);	
			long int subregionId= region_postSegm->GetSubRegionId(j);
			long int subregionIndex_initialSegm= regionIdMap_initialSegm[subregionId]; 
			INFO_LOG("Sub-Region no. "<<j<<" (id="<<subregionId<<", initial_index="<<subregionIndex_initialSegm<<")");
		
			Region* subregion_initialSegm= segmSlicData.GetRegion(subregionIndex_initialSegm);
			if(!subregion_initialSegm){
				ERROR_LOG("Failed to get subregion with index "<<subregionIndex_initialSegm<<" from initial slic data (this should not occur!)");
				continue;
			}
			//(segmSlicData.regions)[subregionIndex_original]->Tag= regionTag;
			subregion_initialSegm->Tag= regionTag;
		}	
	}//end loop regions
	

	INFO_LOG("Remove slic data computed in adaptive merging stage...");
	if(slicData_segmMultiStep){
		delete slicData_segmMultiStep;
		slicData_segmMultiStep= 0;
	}
	
	//## Hierarchically merge signal regions
	INFO_LOG("Hierarchically merge signal regions (#"<<segmSlicData.GetNRegions()<<" regions present)...");
	if(segmSlicData.GetNRegions()>2 && nSig>0){
		bool includeSpatialPars= true;
		int nRegions_beforeMerging= segmSlicData.GetNRegions();
		int mergingStatus= SPHierarchicalMerger(
			segmSlicData, Region::eSignalTag, Region::eSignalTag,
			minMergedSP, SPMergingRatio, SPMergingRegularization, 
			use2ndNeighborsInSPMerging, SPMergingMaxDissRatio,SPMergingMaxDissRatio_2ndNeighbor,SPMergingDissThreshold,
			SPMergingIncludeSpatialPars, SPMergingUseRobustPars, SPMergingUseCurvDist
		);
		if(mergingStatus<0){
			ERROR_LOG("Merging of signal regions failed!");
			return -1;
		}	
		int nRegions_afterMerging= segmSlicData.GetNRegions();
		INFO_LOG("Merged signal regions (NR="<<nRegions_afterMerging-nRegions_beforeMerging<<")");
	}//close if
	

	
	//## Finally add all background regions together
	INFO_LOG("Finally add all background regions together...");
	int bkgLabel= -1;
	Region* bkgRegion= new Region;
	bkgRegion->Id= 0;
	//std::vector<Region*> signalPlusMergedBkgRegions;
	bool isLabelSet= false;
	bool addPixels= true;
	bool copyPixels= true;
	int nBkgRegions= 0;
	int nSignalRegions= 0;
	std::map<long int,long int> relabelMap;
	std::vector<size_t> regionsToBeDeleted; 

	for(int i=0;i<segmSlicData.GetNRegions();i++) {
		Region* thisRegion= segmSlicData.GetRegion(i);
		int tag= thisRegion->Tag;
		long int regionId= thisRegion->Id;
		relabelMap.insert( std::pair<long int,long int>(regionId,regionId) );

		if(tag==Region::eBkgTag) {
			if(!isLabelSet){
				bkgLabel= regionId;
				bkgRegion->Tag= Region::eBkgTag;
				isLabelSet= true;
			}
			else {
				regionsToBeDeleted.push_back(i);
			}

			relabelMap[regionId]= bkgLabel;
		
			bkgRegion->AddRegion(thisRegion,addPixels,copyPixels);		
			nBkgRegions++;
		}//close if bkg regions
		else if(tag==Region::eSignalTag){
			nSignalRegions++;
		}
		else{
			WARN_LOG("Untagged region present (this should not occur CHECK!!!!)");
		}
	}//end loop regions


	bkgRegion->Id= bkgLabel;
	
	//Remove merged regions and compute region pars	
	INFO_LOG("Remove merged regions and compute region pars	...");
	CodeUtils::DeleteItems(segmSlicData.regions, regionsToBeDeleted);
	for(int i=0;i<segmSlicData.GetNRegions();i++){
		segmSlicData.regions[i]->ComputeStats(true,false);
	}//end loop merged regions

	
	//Update pixel labels	
	INFO_LOG("Update pixel labels	...");
	/*
	for(size_t i=0;i<(segmSlicData.labels).size();i++) {
		for(size_t j=0;j<(segmSlicData.labels)[i].size();j++) {
			int oldLabel= (segmSlicData.labels)[i][j];
			int newLabel= relabelMap[oldLabel];
			(segmSlicData.labels)[i][j]= newLabel;
		}
	}
	*/
	for(size_t i=0;i<(segmSlicData.pixel_labels).size();i++) {
		long int label_old= (segmSlicData.pixel_labels)[i];
		long int label_new= relabelMap[label_old];
		(segmSlicData.pixel_labels)[i]= label_new;
	}

	return 0;

}//close FindSegmentation()


int SLICSegmenter::MultiStepSPMerger(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,bool SPMergingIncludeSpatialPars,bool SPMergingUseRobustPars,bool SPMergingUseCurvDist){
	
	//## Check if regions have been tagged
	//## Return if no signal-tagged regions is present
	int nSig= 0;
	int nBkg= 0;
	int nUntagged= 0;
	int status= SLIC::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
	if(status<0 || (nSig==0 && nBkg==0) ){
		ERROR_LOG("Regions are not tagged (hint: you must tag regions before running segmentation)!");
		return -1;
	}	
	if(nSig==0){//Tag all untagged regions as bkg and end!		
		WARN_LOG("No signal-tagged regions available, tagging all regions as background...");
		for(int i=0;i<segmSlicData.GetNRegions();i++) {
			int tag= (segmSlicData.regions)[i]->Tag;
			if(tag==Region::eUntagged) (segmSlicData.regions)[i]->Tag= Region::eBkgTag;
		}//end loop regions	
		return 0;
	}//close if
	 
	//## Run region merging
	bool stopAlgo= false;
	int nTotMergedRegions= 0;
	int stageCounter= 0;

	while(!stopAlgo){

		//## Update tagging stats
		nTotMergedRegions= 0;
		if(stageCounter>0){
			SLIC::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
		}
		stageCounter++;
		INFO_LOG("Stage no. "<<stageCounter<<" NR="<<segmSlicData.GetNRegions()<<": (nBkg,nSig,nUntagged)=("<<nBkg<<","<<nSig<<","<<nUntagged<<")");
		
		//## == STAGE 1==
		//## Merge non-tagged regions with bkg-tagged regions
		if(nBkg>0 && nUntagged>0) {//start 1st stage if there are untagged regions present	
			INFO_LOG("=================");
			INFO_LOG("==   STAGE 1   ==");
			INFO_LOG("=================");
			INFO_LOG("Start 1st stage...");		
			int nregions_preStage1= segmSlicData.GetNRegions();
			SPMaxSimilarityMerger(
				segmSlicData, Region::eBkgTag, Region::eUntagged,
				SPMergingRegularization, use2ndNeighborsInSPMerging,
				SPMergingIncludeSpatialPars, SPMergingUseRobustPars, SPMergingUseCurvDist
			);
			int nregions_postStage1= segmSlicData.GetNRegions();
			int nMergedRegions_1stStage= nregions_preStage1-nregions_postStage1;
			nTotMergedRegions+= nMergedRegions_1stStage;
		}//close if	

		//Update tagging stats
		SLIC::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
		INFO_LOG("After 1st stage: NR="<<segmSlicData.GetNRegions()<<": (nBkg,nSig,nUntagged)=("<<nBkg<<","<<nSig<<","<<nUntagged<<")");
		
		//## == STAGE 2 ==
		//## Merge non-tagged regions adaptively
		if(nUntagged>0){//start 2nd stage if there are untagged regions available
			INFO_LOG("=================");
			INFO_LOG("==   STAGE 2   ==");
			INFO_LOG("=================");
			INFO_LOG("Start 2nd stage...");
			int nregions_preStage2= segmSlicData.GetNRegions();
			SPMaxSimilarityMerger(	
				segmSlicData, Region::eUntagged, Region::eUntagged,
				SPMergingRegularization, use2ndNeighborsInSPMerging,
				SPMergingIncludeSpatialPars, SPMergingUseRobustPars, SPMergingUseCurvDist
			);
			int nregions_postStage2= segmSlicData.GetNRegions();
			int nMergedRegions_2ndStage= nregions_preStage2-nregions_postStage2;
			nTotMergedRegions+= nMergedRegions_2ndStage;	
		}//close if
		
		//Update tagging stats
		SLIC::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
		INFO_LOG("After 2nd stage: NR="<<segmSlicData.GetNRegions()<<": (nBkg,nSig,nUntagged)=("<<nBkg<<","<<nSig<<","<<nUntagged<<")");
	
		//## Check end condition
		if(nUntagged<=0) {
			INFO_LOG("NR="<<segmSlicData.GetNRegions()<<": No untagged regions left, algorithm end!");
			stopAlgo= true;
		}
		if(nTotMergedRegions==0){
			INFO_LOG("NR="<<segmSlicData.GetNRegions()<<": No regions merged in all stages, mark remaining as signal and end algorithm!");
			for(int i=0;i<segmSlicData.GetNRegions();i++) {
				int tag= (segmSlicData.regions)[i]->Tag;
				if(tag==Region::eUntagged) (segmSlicData.regions)[i]->Tag= Region::eSignalTag;
			}//end loop regions
			stopAlgo= true;
		}//close if

	}//end while loop 

	return 0;

}//close MultiStepSPMerger()



int SLICSegmenter::SPMaxSimilarityMerger(SLICData& segmSlicData,int mergerTag,int mergedTag,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,bool SPMergingIncludeSpatialPars,bool SPMergingUseRobustPars,bool SPMergingUseCurvDist){

	//## Check regions
	int nRegions= segmSlicData.GetNRegions();
	if(nRegions<=0) {
		WARN_LOG("No regions available, nothing to be merged!");
		return 0;
	}
	
	//## Compute region contour info (neighbors, ...)
	INFO_LOG("Finding region neighbors (NR="<<nRegions<<") ...");
	SLICContourData* contourData= SLIC::ComputeBoundaryContours(&segmSlicData);
	if(!contourData){
		ERROR_LOG("Failed to compute the slic contour data!");
		return -1;
	}
	
	//## Find neighbors list
	SLICNeighborCollections neighbors;
	int selectedTag= -1;//mergedTag
	bool normalizeParams= true;
	
	int status= SLIC::FindNeighbors(neighbors,&segmSlicData,contourData,use2ndNeighborsInSPMerging,selectedTag,SPMergingIncludeSpatialPars,normalizeParams,SPMergingUseRobustPars,SPMergingUseCurvDist);
	if(status<0){
		ERROR_LOG("Failed to find neighbors list!");
		return -1;
	}

	//Delete contour data		
	if(contourData){
		delete contourData;
		contourData= 0;
	}


	//## Compute similarity matrix	
	INFO_LOG("Compute region similarity matrix...");
	SLICSimilarityData* similarityData= SLIC::ComputeRegionSimilarity(&segmSlicData,neighbors,SPMergingRegularization);
	if(!similarityData){
		ERROR_LOG("Failed to compute region similarity matrix, cannot perform region merging!");
		return -1;
	}
	TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
	
	//## Compute page rank of segments and sort
	INFO_LOG("Compute page rank ...");
	std::vector<double> ranks;
	status= StatsUtils::ComputePageRank(ranks,AdjacencyMatrix->T());//pass transpose of adjacency matrix
	if(status<0 || ranks.empty()){
		WARN_LOG("PageRank failed, cannot perform region merging!");
		delete similarityData;
		similarityData= 0;
		return -1;
	}

	//## Delete similarity data
	DEBUG_LOG("Deleting similarity data...");
	if(similarityData){
		delete similarityData;
		similarityData= 0;
	}

	//## Sort ranks
	std::vector<size_t> sort_index;
	std::vector<double> ranks_sorted;
	CodeUtils::sort_descending(ranks,ranks_sorted,sort_index);
		
	
	//## Loop over sorted ranks and select regions to be merged
	INFO_LOG("Start region merging loop ...");
	std::vector< std::vector<int> > regionsToBeMerged;
	for(size_t i=0;i<(segmSlicData.regions).size();i++) regionsToBeMerged.push_back( std::vector<int>() );
	std::vector<size_t> regionsToBeDeleted;
	std::vector<size_t> regionsIdToBeDeleted;
	int nMergedRegions= 0;

	for(size_t i=0;i<sort_index.size();i++){
		size_t index= sort_index[i];//region index in regions list
		long int regionId= (segmSlicData.regions)[index]->Id;
		int regionTag= (segmSlicData.regions)[index]->Tag;
		
		//Skip if region tag is different from merger tag 
		if(regionTag!=mergerTag && mergerTag!=-1) {
			DEBUG_LOG("Skip region no. "<<i<<" (id="<<regionId<<", tag="<<regionTag<<")...");
			continue;
		}

		//Check if this seed was not already merged by a previous (best ranked) region
		INFO_LOG("Checking if seed region "<<index<<" (id="<<regionId<<") was not already merged by a previous (best ranked) region");
		std::vector<size_t>::iterator seedfinderIt= std::find(regionsToBeDeleted.begin(),regionsToBeDeleted.end(),index);
		if(seedfinderIt!=regionsToBeDeleted.end()){
			INFO_LOG("Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!");
			continue;
		}

		//Loop over untagged neighbors and find best merging
		int nGoodNeighbors= 0;
		std::vector<SLICNeighborData> thisNeighbors= neighbors[index].GetNeighbors();

		for(size_t j=0;j<thisNeighbors.size();j++){//loop over neighbors of region with ranked index

			int neighborIndex= thisNeighbors[j].Index;
			int neighborId= (segmSlicData.regions)[neighborIndex]->Id;
			int neighborTag= (segmSlicData.regions)[neighborIndex]->Tag;

			if(mergerTag==Region::eBkgTag && neighborTag==mergerTag){//for bkg-untagged merging look at neighbors with tag different than bkg
				continue;
			}
			if(mergerTag==Region::eUntagged && neighborTag!=mergedTag){//for untagged-untagged merging look at neighbors with 'untagged' tag
				continue;
			}
			nGoodNeighbors++;
			
			//double Delta_ij= (*DissimilarityMatrix)(index,neighborIndex);
			double Delta_ij= thisNeighbors[j].Dsym;
	
			//Loop over neighbors of this neighbors (2nd neighbors)
			int closerIndex= -1;
			double closerDiss= 1.e+99;	
			std::vector<SLICNeighborData> thisNeighbors_2nd= neighbors[neighborIndex].GetNeighbors();

			for(size_t k=0;k<thisNeighbors_2nd.size();k++){
			
				int neighborIndex_2nd= thisNeighbors_2nd[k].Index;
				double Delta_jk= thisNeighbors_2nd[k].Dsym; 

				if(Delta_jk<closerDiss && Delta_jk>0){	
					closerDiss= Delta_jk;	
					closerIndex= neighborIndex_2nd;
				}
			}//end loop 2nd neighbors

			//If diss Dij is the maximum among Djk(s) merge ij!
			if(closerIndex!=(int)index) {
				INFO_LOG("Skip merging of seed region "<<index<<" (id="<<regionId<<") with region "<<neighborIndex<<" (id="<<neighborId<<") as not the best possible merging for merged region...");
				continue;
			}


			//Check if this closer region was not already selected to be merged to another node previously	
			std::vector<size_t>::iterator finderIt= std::find(regionsToBeDeleted.begin(),regionsToBeDeleted.end(),neighborIndex);				
			if(finderIt!=regionsToBeDeleted.end()){	
				INFO_LOG("Closer neighbor no. "<<neighborIndex<<" (id="<<neighborId<<") was already selected for merging in a previous node, skip this merging!");
				continue;
			}

			//Check if this neighbor was not previously selected as a merger
			if(regionsToBeMerged[neighborIndex].size()>0){
				INFO_LOG("Closer neighbor no. "<<neighborIndex<<" (id="<<neighborId<<") was previously selected as a merger, skip this merging!");
				continue;
			}

			//Merging selected!
			INFO_LOG("Regions ("<<index<<"-"<<neighborIndex<<") selected for merging (Delta_ij="<<Delta_ij<<" closerDiss="<<closerDiss<<")");
			regionsToBeMerged[index].push_back(neighborIndex);
			regionsToBeDeleted.push_back(neighborIndex);
			regionsIdToBeDeleted.push_back(neighborId);		
			nMergedRegions++;

		}//end loop neighbors
				
	}//end loop ranks
		
	/*
	//## Delete similarity data
	DEBUG_LOG("Deleting similarity data...");
	if(similarityData){
		delete similarityData;
		similarityData= 0;
	}
	*/
		
	//## Merge regions and update region and label list
	if(nMergedRegions>0){
		INFO_LOG("Merge the selected regions...");
		std::map<long int,long int> newLabelMap;
		bool copyPixels= true;
		bool addPixels= true;

		for(size_t k=0;k<(segmSlicData.regions).size();k++){
			long int regionId= (segmSlicData.regions)[k]->Id;	
			int regionTag= (segmSlicData.regions)[k]->Tag;
			newLabelMap.insert( std::pair<long int,long int>(regionId,regionId) );//Initialize to same label

			for(size_t j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionIndex= regionsToBeMerged[k][j];
				long int mergedRegionId= (segmSlicData.regions)[mergedRegionIndex]->Id;
				int mergedRegionTag= (segmSlicData.regions)[mergedRegionIndex]->Tag;
				newLabelMap[mergedRegionId]= regionId;

				INFO_LOG("Region no. "<<k<<" (id="<<regionId<<", tag="<<regionTag<<"): merging region id="<<mergedRegionId<<", tag="<<mergedRegionTag<<" (index="<<mergedRegionIndex<<")");			
				(segmSlicData.regions)[k]->AddRegion((segmSlicData.regions)[mergedRegionIndex],addPixels,copyPixels);//NB: This also adds the sub-regions id
				//(segmSlicData.regions)[k]->AddSubRegionId(mergedRegionId);//No need to call this!!!
			
			}//end loop regions to be merged
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		INFO_LOG("Deleting regions aggregated in this step from the main list...");
		CodeUtils::DeleteItems((segmSlicData.regions), regionsToBeDeleted);
		
		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		INFO_LOG("Updating region parameters & contours...");
		bool forceRecomputing= true;
		for(unsigned int k=0;k<(segmSlicData.regions).size();k++){
			(segmSlicData.regions)[k]->ComputeStats(SPMergingUseRobustPars,forceRecomputing);
		}

		//## Update pixel labels
		INFO_LOG("Updating pixel labels...");
		for(size_t i=0;i<(segmSlicData.pixel_labels).size();i++) {
			long int label_old= (segmSlicData.pixel_labels)[i];
			long int label_new= newLabelMap[label_old];
			(segmSlicData.pixel_labels)[i]= label_new;
		}

	}//close if merge regions
		
	INFO_LOG(nMergedRegions<<"/"<<segmSlicData.GetNRegions()<<" regions merged at this stage...");

	return 0;

}//close SPMaxSimilarityMerger()


int SLICSegmenter::SPHierarchicalMerger(SLICData& slicData,int mergerTag,int mergedTag,int minMergedSP,double SPMergingRatio,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,double SPMergingMaxDissRatio,double SPMergingMaxDissRatio_2ndNeighbor,double SPMergingDissThreshold, bool SPMergingIncludeSpatialPars,bool SPMergingUseRobustPars,bool SPMergingUseCurvDist){

	//## Check regions
	int nRegions= slicData.GetNRegions();
	if(nRegions<=0) {
		WARN_LOG("No regions available, nothing to be merged!");
		return 0;
	}

	//----
	INFO_LOG("slic data edge image stats");
	(slicData.edgeImg)->PrintStats();
	//---

	//## Run a hierarchical clustering till all segments are merged in one
	int hierarchyLevel= 0;
	double DissMedian0= 0;
	int nMergedRegionsInHierarchyLevel= 0;
	int itemPos= -1;
	bool normalizeParams= true;
	
	INFO_LOG("Starting hierarchical merging...");

	while( (slicData.GetNRegions())>minMergedSP ){
		nMergedRegionsInHierarchyLevel= 0;
		
		//## Create the mapping of regionId and vector index
		//INFO_LOG("Create the mapping of regionId and region list index...");
		//std::map<long int,long int> regionIdMap;
		//slicData.GetRegionIdMap(regionIdMap);

		//## Compute region contour info (neighbors, ...)
		INFO_LOG("Finding region neighbors at hierarchy level "<<hierarchyLevel<<", NR="<<slicData.GetNRegions()<<" ...");
		SLICContourData* contourData= SLIC::ComputeBoundaryContours(&slicData);
		
		//## Find neighbors list
		SLICNeighborCollections neighbors;
		int selectedTag= mergedTag;
		int status= SLIC::FindNeighbors(neighbors,&slicData,contourData,use2ndNeighborsInSPMerging,selectedTag,SPMergingIncludeSpatialPars,normalizeParams,SPMergingUseRobustPars,SPMergingUseCurvDist);
		if(status<0){
			ERROR_LOG("Failed to find neighbors list!");
			return -1;
		}
		if(contourData){
			delete contourData;
			contourData= 0;
		}
	
		
		//## Compute similarity matrix	
		INFO_LOG("Compute region similarity at hierarchy level "<<hierarchyLevel<<" ...");
		SLICSimilarityData* similarityData= SLIC::ComputeRegionSimilarity(&slicData,neighbors,SPMergingRegularization);
		if(!similarityData){
			ERROR_LOG("Failed to compute region similarity matrix, cannot perform region hierarchical merging!");
			return -1;
		}
		TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
		double DissMedian= similarityData->Dmedian;
 		if(hierarchyLevel==0) {
			DissMedian0= DissMedian;
		}


		//## Compute page rank of segments and sort
		INFO_LOG("Compute page rank at hierarchy level "<<hierarchyLevel<<" ...");
		std::vector<double> ranks;
		status= StatsUtils::ComputePageRank(ranks,AdjacencyMatrix->T());//pass transpose of adjacency matrix
		if(status<0 || ranks.size()<=0){
			WARN_LOG("PageRank failed, cannot perform region merging!");
			return -1;
		}
		
		//## Delete similarity data
		DEBUG_LOG("Deleting similarity data...");
		if(similarityData){
			delete similarityData;
			similarityData= 0;
		}

		//## Sort ranks
		std::vector<size_t> sort_index;//sorting index
		std::vector<double> ranks_sorted;
		CodeUtils::sort_descending(ranks,ranks_sorted,sort_index);
		
		//## Count max number of mergeable regions
		int nMaxMergeableRegions= 0;
		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			int regionTag= slicData.regions[index]->Tag;
			if(regionTag!=mergerTag && mergerTag!=-1) continue;
			nMaxMergeableRegions++;
		}//end loop regions
		int maxRegionsToMerge= std::round(nMaxMergeableRegions*SPMergingRatio);
		INFO_LOG(maxRegionsToMerge<<" regions can be merged at maximum at this hierarchy level ("<<hierarchyLevel<<") given the number of mergeable regions found (N="<<nMaxMergeableRegions<<") and the chosen merging ratio option (ratio="<<SPMergingRatio<<")");

		//## Loop over sorted ranks and select regions to be merged
		int nMergedRegions= 0;
		int nMergeableRegions= 0;
		std::vector< std::vector<int> > regionsToBeMerged;
		for(int i=0;i<slicData.GetNRegions();i++) regionsToBeMerged.push_back( std::vector<int>() );
		std::vector<long int> regionsToBeDeleted;
		std::vector<long int> regionsIdToBeDeleted;

		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			//double thisRank= ranks[index];//region rank
			long int regionId= (slicData.regions)[index]->Id;
			int regionTag= (slicData.regions)[index]->Tag;
			//int mapIndex= regionIdMap[regionId];
			
			if(regionTag!=mergerTag && mergerTag!=-1) {
				INFO_LOG("Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as tag is different from merger region (tag="<<mergerTag<<")");
				continue;
			}

			//Stop merging above nmerge threshold						
			if(nMergeableRegions+1>maxRegionsToMerge) {
				INFO_LOG("Maximum number of regions mergeable reached for this hierarchy level ("<<nMergeableRegions+1<<">"<<maxRegionsToMerge<<")");
				break;
			}
			
			//Check if this seed was not already merged by a previous (best ranked) region
			//NB: If it is found int the regionsIdToBeDeleted collection skip the merging
			if(CodeUtils::FindItem(regionsIdToBeDeleted,regionId,itemPos)){
				INFO_LOG("Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!");
				continue;
			}


			//Chech if the seed region has any neighbors (according to the selected merged tag)
			int NN= (int)neighbors[index].GetN();
			if(NN<=0){
				INFO_LOG("Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as no neighbors are available!");
				continue;
			}
			/*
			int NN_1st= (int)neighborIndexList_1st[index].size();
			int NN_2nd= (int)neighborIndexList_2nd[index].size();	
			if(NN_1st<=0 && NN_2nd<=0){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as no neighbors are available!"<<endl;
				continue;
			}
			*/

			//Get closest neighbor
			/*
			if(DissSortIndexes[index].size()<=0) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip seed ranked region (id="<<regionId<<") as no neighbors are available!"<<endl;
				continue;
			}
			int closerNeighborIndex= DissSortIndexes[index][0];
			int closerNeighborId= (slicData.regions)[closerNeighborIndex]->Id;
			int closerNeighborTag= (slicData.regions)[closerNeighborIndex]->Tag;
			*/

			int closerNeighborPosInCollection= neighbors[index].FindCloserByDissTot();//NB: this is not the region index in original vector!!!
			SLICNeighborData* closerNeighbor= neighbors[index].GetNeighbor(closerNeighborPosInCollection);
			long int closerNeighborId= closerNeighbor->Id;
			int closerNeighborTag= closerNeighbor->Tag; 
			int closerNeighborOrder= closerNeighbor->Order;
			int closerNeighborIndex= closerNeighbor->Index;//regionIdMap[closerNeighborId];
			
			if(closerNeighborId==regionId) {//skip same region
				INFO_LOG("Skip neighbor (id="<<closerNeighborId<<") as equal to seed region!");
				continue;
			}
			if(closerNeighborTag!=mergedTag && mergedTag!=-1) {	
				INFO_LOG("Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as neighbor (id="<<closerNeighborId<<",tag="<<closerNeighborTag<<") tag is different from merged region (tag="<<mergedTag<<")");
				continue;
			}
			nMergeableRegions++;
			
			
			double Delta_ij= closerNeighbor->Dtot;
			double Delta_ji= closerNeighbor->Dtot_n;
			double AbsDelta_ij= closerNeighbor->D;
			
			//Find current region index among best neighbors of closer
			int N_half= ceil(neighbors[closerNeighborIndex].GetN()/2.);
			if(N_half>=2){
				bool isSeedAmongNeighborCloser= neighbors[closerNeighborIndex].IsIdAmongNClosersByTotDiss(regionId,N_half);
				if(!isSeedAmongNeighborCloser){
					INFO_LOG("Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<" since this is not among the closest neighbors of the selected neighbor!");
					continue;
				}
			}

			/*
			double Delta_ij= (*DissimilarityMatrix)(index,closerNeighborIndex);//include edgeness
			double Delta_ji= (*DissimilarityMatrix)(closerNeighborIndex,index);//include edgeness
			double AbsDelta_ij= (*AbsDissimilarityMatrix)(index,closerNeighborIndex);//without edgeness
			double AbsDelta_ji= (*AbsDissimilarityMatrix)(closerNeighborIndex,index);//without edgeness
			int closerNeighborness= (*NeighborMatrix)(index,closerNeighborIndex);
			
			//Get neighbors of the closer region
			std::vector<int> closerNeighborNNIndexes;
			closerNeighborNNIndexes.assign(
				DissSortIndexes[closerNeighborIndex].begin(),
				DissSortIndexes[closerNeighborIndex].begin() + ceil(DissSortIndexes[closerNeighborIndex].size()/2.)
			);

			//Find current region index among best neighbors of closer
			std::vector<int>::iterator it= std::find(closerNeighborNNIndexes.begin(),closerNeighborNNIndexes.end(),index);
			if(it==closerNeighborNNIndexes.end() ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<" since this is not among the closest neighbors of selected neighbor!"<<endl;
				continue;
			}
			*/
			
			//## Check similarity difference
			//## 1st neighbors are merged (flood-fill)
			//## 2nd neighbors are merged if their similarities are not too different
			/*
			if(closerNeighborness==1 && (Delta_ji>Delta_ij*SPMergingMaxDissRatio) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (1st neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			else if(closerNeighborness==2 && (Delta_ji>Delta_ij*SPMergingMaxDissRatio_2ndNeighbor) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			*/

			//Check if this closer region was not already selected as a seed merger previously
			if( regionsToBeMerged[closerNeighborIndex].size()>0 ){	
				INFO_LOG("Closer neighbor (id="<<closerNeighborId<<") selected for merging was before selected as a primary merger, skip this merging!");
				continue;
			}

			//Check if this closer region was not already selected to be merged to another node previously	
			if(CodeUtils::FindItem(regionsIdToBeDeleted,closerNeighborId,itemPos)){
				INFO_LOG("Closer neighbor (id="<<closerNeighborId<<") was already selected for merging in a previous node, skip this merging!");
				continue;
			}

			
			//## Apply dissimilarity threshold			
			if(use2ndNeighborsInSPMerging && closerNeighborOrder==2){
				double DissThreshold= SPMergingDissThreshold*DissMedian0;
			
				if(Delta_ji>Delta_ij*SPMergingMaxDissRatio_2ndNeighbor){
					INFO_LOG("Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<SPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij);
					continue;
				}
				if(AbsDelta_ij>DissThreshold){
					INFO_LOG("Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (AbsDiss="<<AbsDelta_ij<<", Diss="<<Delta_ij<<">"<<DissThreshold<<")");
					continue;
				}	
			}//close if


			/*
			if(closerNeighborness==2 && (Delta_ji>Delta_ij*SPMergingMaxDissRatio_2ndNeighbor) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}

			if(closerNeighborness==2 && AbsDelta_ij>DissThreshold){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (AbsDiss="<<AbsDelta_ij<<", Diss="<<Delta_ij<<">"<<DissThreshold<<")"<<endl;
				continue;
			}
			*/

			INFO_LOG("Region (id="<<regionId<<") merging neighbor (id="<<closerNeighborId<<"): order="<<closerNeighborOrder<<", Delta_ij="<<Delta_ij<<", Delta_ji="<<Delta_ji<<" ratio="<<Delta_ji/Delta_ij);

			
			//int regionIndex_A= regionIdMap_top[regionId];
			//int regionIndex_B= regionIdMap_top[closerNeighborId];
			
			regionsToBeMerged[index].push_back(closerNeighborIndex);
			regionsToBeDeleted.push_back(closerNeighborIndex);
			regionsIdToBeDeleted.push_back(closerNeighborId);		
			//mergedRegionList[regionIndex_A].push_back(regionIndex_B);
	
			nMergedRegions++;
		}//end loop ranks

		
		/*
		//## Delete similarity data
		if(similarityData){
			delete similarityData;
			similarityData= 0;
		}
		*/
		
		
		//## Merge regions
		INFO_LOG("Merge the selected regions at hierarchy level "<<hierarchyLevel<<" ...");
		std::map<long int,long int> newLabelMap;
		bool copyPixels= true;
		bool addPixels= true;

		for(int k=0;k<slicData.GetNRegions();k++){
			long int regionId= (slicData.regions)[k]->Id;
			newLabelMap.insert( std::pair<long int,long int>(regionId,regionId) );

			for(size_t j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionIndex= regionsToBeMerged[k][j];
				long int mergedRegionId= (slicData.regions)[mergedRegionIndex]->Id;
				newLabelMap[mergedRegionId]= regionId;
				INFO_LOG("Region no. "<<k<<" (id="<<regionId<<") : merging region id="<<mergedRegionId<<" (index="<<mergedRegionIndex<<")");		
				(slicData.regions)[k]->AddRegion((slicData.regions)[mergedRegionIndex],addPixels,copyPixels);
			}//end loop regions to be merged
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		INFO_LOG("Deleting regions aggregated in this step from the main list...");
		CodeUtils::DeleteItems((slicData.regions), regionsToBeDeleted);
		//for(size_t k=0;k<regionsIdToBeDeleted.size();k++) regionIdMap.erase(regionsIdToBeDeleted[k]);

		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		INFO_LOG("Updating region parameters & contours...");
		for(int k=0;k<slicData.GetNRegions();k++){
			(slicData.regions)[k]->ComputeStats(false,true);
			//int regionId= (slicData.regions)[k]->Id;
			//regionIdMap[regionId]= k;
		}//end loop regions
	
		//## Update pixel labels
		INFO_LOG("Updating pixel labels...");
		for(size_t i=0;i<(slicData.pixel_labels).size();i++) {
			long int label_old= (slicData.pixel_labels)[i];
			long int label_new= newLabelMap[label_old];
			(slicData.pixel_labels)[i]= label_new;
		}

		nMergedRegionsInHierarchyLevel= nMergedRegions;
		hierarchyLevel++;

		if(nMergedRegionsInHierarchyLevel==0){
			WARN_LOG("No regions merged in this stage, exit to avoid stuck!");
			break;
		}
		INFO_LOG(nMergedRegionsInHierarchyLevel<<"/"<<slicData.GetNRegions()<<" regions aggregated at this level hierarchy...");
		

	}//end while loop


	INFO_LOG("#"<<hierarchyLevel<<" hierarchy levels aggregated (N="<<slicData.GetNRegions()<<" regions left)");


	/*
	//## Final tag of composite regions according to a majority criterion, that is region tagged as bkg if the majority of sub-regions are non-significative
	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Tagging significant regions..."<<endl;
	for(int i=0;i<regions.size();i++){
		int regionId= regions[i]->fId;
		int regionTag= regions[i]->fTag;
		int regionIndex_top= regionIdMap_top[regionId];//region index at the top of the hierarchy
		std::vector<int> subregionIds= regions[i]->fSubRegionIds;

		bool isSignificant= inputRegions[regionIndex_top]->fIsSignificative;
		bool isSalient= inputRegions[regionIndex_top]->fIsSalient;
		double salientFraction= 0;
		int nSubRegions= 1+subregionIds.size();
		if(isSalient) salientFraction++;
		bool isAnySignificant= false;
		if(isSignificant) isAnySignificant= true;

		for(unsigned int j=0;j<subregionIds.size();j++){
			int subRegionId= subregionIds[j];
			int subRegionIndex_top= regionIdMap_top[subRegionId];
			//bool subIsSignificant= fRegions[subRegionIndex_top]->fIsSignificative;
			//bool subIsSalient= fRegions[subRegionIndex_top]->fIsSalient;	
			bool subIsSignificant= inputRegions[subRegionIndex_top]->fIsSignificative;
			bool subIsSalient= inputRegions[subRegionIndex_top]->fIsSalient;
			if(subIsSalient) salientFraction++;
			if(subIsSignificant) isAnySignificant= true;
		}
		salientFraction/= (double)nSubRegions;
			
		if( (salientFraction>fSignificantSPRatio) || isAnySignificant ){
			cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<" tagged as significant as (significanceRatio="<<salientFraction<<">"<<fSignificantSPRatio<<")..."<<endl;
			regions[i]->fIsSignificative= true;
		}
		else {
			regions[i]->fIsSignificative= false;
		}	
	}//end loop regions
		
	//## Tag background?
	if(fUsePixelRatioCut){
		//Set region with maximum number of pixels
		for(int i=0;i<regions.size();i++){
			int regionId= regions[i]->fId;	
			int npix= regions[i]->fNPix;
			double pixelRatio= (double)npix/(double)N;
			if(pixelRatio>fPixelRatioCut) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<" has too many pixels (ratio="<<pixelRatio<<">"<<fPixelRatioCut<<") and will be tagged as bkg..."<<endl;
				regions[i]->fIsSignificative= false;
			}
		}
	}//close if
	
	//## Copy final results to global data (fRegions, fPixelClusterIds)
	inputRegions.clear();
	inputRegions.assign(regions.begin(),regions.end());
	for(unsigned int i=0;i<labels.size();i++) {
		for(unsigned int j=0;j<labels[i].size();j++) {
			inputLabels[i][j]= labels[i][j];
		}
	}
	*/	
	
	return 0;

}//close SPHierarchicalMerger()

}//close namespace


