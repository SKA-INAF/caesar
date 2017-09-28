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
* @file SourceFinder.cc
* @class SourceFinder
* @brief Source finder class
*
* Class to perform source finding 
* @author S. Riggi
* @date 20/01/2015
*/

#include <SFinder.h>
#include <BlobFinder.h>
#include <Image.h>
#include <Source.h>
#include <Contour.h>
#include <ConfigParser.h>
#include <BkgData.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <SysUtils.h>
#include <Logger.h>
#include <Consts.h>

#include <SLIC.h>
#include <SLICSegmenter.h>
#include <ChanVeseSegmenter.h>
#include <LRACSegmenter.h>

#include <TaskData.h>
#include <Serializer.h>
#include <Graph.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TCanvas.h>

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
#include <chrono>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

using namespace std;

ClassImp(Caesar::SFinder)

namespace Caesar {

#define MASTER_ID 0
#define MAX_NTASKS_PER_WORKER 10000
 
SFinder::SFinder() 
{
	
}//close costructor


SFinder::~SFinder(){
	
	//Clearup
	INFO_LOG("[PROC "<<m_procId<<"] - Clearup source finder allocated data...");
	Clear();

}//close destructor


void SFinder::Clear()
{

	//## Close open ROOT file
	//## NB: When objects are written to file their memory is released, so don't delete them!
	if(m_procId==MASTER_ID){
		if( m_OutputFile && m_OutputFile->IsOpen() ) {
			DEBUG_LOG("[PROC "<<m_procId<<"] - Closing output ROOT file...");
			m_OutputFile->Close();
		}
	}

	//## Delete source tree
	if(m_procId==MASTER_ID){
		if(m_SourceTree && !m_saveSources) {
			DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting source tree...");
			m_SourceTree->Delete();
		}
	}

	//## Delete images & objects not written to file
	if(m_InputImg && !m_saveInputMap){	
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting input image...");
		delete m_InputImg;
		m_InputImg= 0;
	}
	if(m_BkgData && !m_saveBkgMap && !m_saveNoiseMap){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting bkg data...");
		delete m_BkgData;
		m_BkgData= 0;
	}
	if(m_ResidualImg && !m_saveResidualMap) {
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting residual image...");
		delete m_ResidualImg;
		m_ResidualImg= 0;
	}
	if(m_SignificanceMap && !m_saveSignificanceMap) {
		delete m_SignificanceMap;
		m_SignificanceMap= 0;
	}			
	if(m_EdgeImg && !m_saveEdgenessMap){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting edgeness image...");
		delete m_EdgeImg;
		m_EdgeImg= 0;
	}
	if(m_LaplImg && !m_saveCurvatureMap){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting Laplacian image...");
		delete m_LaplImg;
		m_LaplImg= 0;
	}
	if(m_SegmImg && !m_saveSegmentedMap){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting Segmented image...");
		delete m_SegmImg;
		m_SegmImg= 0;
	}
	if(m_SaliencyImg && !m_saveSaliencyMap){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting Saliency image...");
		delete m_SaliencyImg;
		m_SaliencyImg= 0;
	}
	

	//## Delete TApplication??	
	//if(m_Application) {
	//	m_Application->Delete();
	//}

	//## Delete task data
	DEBUG_LOG("[PROC "<<m_procId<<"] - Deleting task data...");
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++) {
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++) {
			if(m_taskDataPerWorkers[i][j]){
				delete m_taskDataPerWorkers[i][j];
				m_taskDataPerWorkers[i][j]= 0;
			}
		}
		m_taskDataPerWorkers[i].clear();
	}
	m_taskDataPerWorkers.clear();


	//## Free comm & groups	
	#ifdef MPI_ENABLED
		if(m_mpiEnabled){
			DEBUG_LOG("[PROC "<<m_procId<<"] - Freeing worker & world MPI groups...");
			MPI_Group_free(&m_WorkerGroup);
			if(m_WorkerComm!=MPI_COMM_NULL) MPI_Comm_free(&m_WorkerComm);
			MPI_Group_free(&m_WorldGroup);
		}
	#endif


}//close Clear()


void SFinder::InitOptions()
{

	//Check is MPI run is enabled at build & runtime
	m_mpiEnabled= SysUtils::IsMPIInitialized();
	INFO_LOG("MPI enabled? "<<m_mpiEnabled);

	//MPI vars (initialized after to actual value)
	//NB: When MPI is not used this should define only 1 process and 1 master
	m_nProc= 1;
	m_procId= 0;

	#ifdef MPI_ENABLED
	if(m_mpiEnabled){
		MPI_Comm_size(MPI_COMM_WORLD, &m_nProc);
		MPI_Comm_rank(MPI_COMM_WORLD, &m_procId);
	}
	#endif
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<m_nProc<<" processors for this run...");
	
	//Input file options
	m_InputFileName= "";
	m_InputImgName= "";
	m_InputFileExtension= "";
	m_InputImg= 0;

	//Beam info
	m_beamFWHM= 6.5;
	m_pixSize= 1;
	m_fluxCorrectionFactor= 1;
	
	//Output options
	m_OutputFile= 0;	
	m_OutputFileName= "";
	m_Application= 0;
	m_IsInteractiveRun= false;
	m_saveToFile= true;
	m_saveConfig= true;
	m_saveSources= true;	
	m_saveResidualMap= true;
	m_saveInputMap= true;
	m_saveSaliencyMap= false;
	m_saveEdgenessMap= false;
	m_saveCurvatureMap= false;
	m_saveSegmentedMap= true;
	m_SourceTree= 0;
	m_saveDS9Region= true;
	m_DS9CatalogFileName= "";
	m_DS9RegionFormat= 1;
	m_PerfTree= 0;
		
	//Source 
	m_Source= 0;
	m_SourceCollection.clear();	
	m_SourcesMergedAtEdges.clear();
	//m_CompactSources.clear();
	//m_ExtendedSources.clear();

	//Read options
	m_ReadTile= false;
	m_TileMinX= 0;
	m_TileMaxX= 0;
	m_TileMinY= 0;
	m_TileMaxY= 0;
	m_mergeSourcesAtEdge= true;
	
	//Bkg data
	m_BkgData= 0;
	m_SignificanceMap= 0;

	//Residual img
	m_ResidualImg= 0;

	//Saliency img
	m_SaliencyImg= 0;

	//Performance stats
	totTime= 0;
	initTime= 0;
	readImageTime= 0;
	compactSourceTime= 0;
	sourceSelectionTime= 0;
	imgResidualTime= 0;
	extendedSourceTime= 0;
	sourceFitTime= 0;	
	saveTime= 0;

	//Extended source finder
	m_EdgeImg= 0;
	m_LaplImg= 0;
	m_SegmImg= 0;

	
	//Task info
	m_TaskInfoTree= 0;


}//close InitOptions()


int SFinder::Init(){

	//## Init options
	InitOptions();

	//## Configure from parser
	if(Configure()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to configure options from parser!");
		return -1;
	}

	//## Create TApplication if interactive run is selected
	if(!m_Application && m_IsInteractiveRun){
		m_Application= new TApplication("Application", 0, 0);
	}	

	//## Create output file
	//## NB: Done only by processor 0 in MPI run
	if(m_saveToFile && m_procId==MASTER_ID){
		INFO_LOG("[PROC "<<m_procId<<"] - Opening ROOT output file "<<m_OutputFileName<<" ...");
		if(!m_OutputFile) m_OutputFile= new TFile(m_OutputFileName.c_str(),"RECREATE");	
		m_OutputFile->cd();
	
		//Init source tree
		if(m_saveSources){
			m_Source= 0;
			if(!m_SourceTree) {
				INFO_LOG("[PROC "<<m_procId<<"] - Creating ROOT source tree ...");	
				m_SourceTree= new TTree("SourceInfo","SourceInfo");
			}
			m_SourceTree->Branch("Source",&m_Source);
			m_SourceCollection.clear();
		}

		//Init task info tree
		if(!m_TaskInfoTree) {
			INFO_LOG("[PROC "<<m_procId<<"] - Creating ROOT task data tree ...");	
			m_TaskInfoTree= new TTree("TaskInfo","TaskInfo");
		}
		m_TaskInfoTree->Branch("xmin",&m_xmin);
		m_TaskInfoTree->Branch("xmax",&m_xmax);
		m_TaskInfoTree->Branch("ymin",&m_ymin);
		m_TaskInfoTree->Branch("ymax",&m_ymax);
	

		//Init time performance tree
		if(!m_PerfTree) {
			INFO_LOG("[PROC "<<m_procId<<"] - Creating ROOT run performance stats tree ...");	
			m_PerfTree= new TTree("PerformanceInfo","PerformanceInfo");
		}
		m_PerfTree->Branch("tot",&totTime,"tot/D");
		m_PerfTree->Branch("init",&initTime,"init/D");
		m_PerfTree->Branch("read",&readImageTime,"read/D");
		m_PerfTree->Branch("sfinder",&compactSourceTime,"sfinder/D");
		m_PerfTree->Branch("sselector",&sourceSelectionTime,"sselector/D");
		m_PerfTree->Branch("imgres",&imgResidualTime,"imgres/D");
		m_PerfTree->Branch("extsfinder",&extendedSourceTime,"extsfinder/D");
		m_PerfTree->Branch("sfit",&sourceFitTime,"sfit/D");
		//m_PerfTree->Branch("save",&saveTime,"save/D");

	}//close if saveToFile


	//## Init and fill task data
	//## NB: Done by all processors in MPI run
	INFO_LOG("[PROC "<<m_procId<<"] - Initializing and filling task data ...");	
	for(int i=0;i<m_nProc;i++){
		m_taskDataPerWorkers.push_back( std::vector<TaskData*>() );
	}
	if(PrepareWorkerTasks()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Preparation of tasks per worker failed!");
		return -1;
	}

	return 0;

}//close Init()

int SFinder::Configure(){

	//Get image read options
	if(GET_OPTION_VALUE(inputFile,m_InputFileName)<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to get inputFile option!");
		return -1;
	}
	GET_OPTION_VALUE(inputImage,m_InputImgName);
	GET_OPTION_VALUE(fitsHDUId,m_fitsHDUId);
	GET_OPTION_VALUE(readTileImage,m_ReadTile);
	if(m_ReadTile){
		GET_OPTION_VALUE(tileMinX,m_TileMinX);
		GET_OPTION_VALUE(tileMaxX,m_TileMaxX);
		GET_OPTION_VALUE(tileMinY,m_TileMinY);
		GET_OPTION_VALUE(tileMaxY,m_TileMaxY);
	}
	else {
		m_TileMinX= 0;
		m_TileMaxX= 0;
		m_TileMinY= 0;
		m_TileMaxY= 0;
	}

	//Get distributed image options
	GET_OPTION_VALUE(splitInTiles,m_splitInTiles);
	GET_OPTION_VALUE(tileSizeX,m_TileSizeX);
	GET_OPTION_VALUE(tileSizeY,m_TileSizeY);
	GET_OPTION_VALUE(useTileOverlap,m_UseTileOverlap);
	GET_OPTION_VALUE(tileStepSizeX,m_TileStepSizeX);
	GET_OPTION_VALUE(tileStepSizeY,m_TileStepSizeY);
	
	if(m_splitInTiles && (m_TileSizeX<=0 || m_TileSizeY<=0) ) {
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid tileSizeX/tileSizeY options!");
		return -1;	
	}
	if(m_splitInTiles && (m_TileStepSizeX<=0 || m_TileStepSizeY<=0 || m_TileStepSizeX>1 || m_TileStepSizeY>1)){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid tileStepSizeX/tileStepSizeY options!");
		return -1;
	}
	if(!m_UseTileOverlap){
		m_TileStepSizeX= 1;
		m_TileStepSizeY= 1;	
	}

	GET_OPTION_VALUE(mergeSourcesAtEdge,m_mergeSourcesAtEdge);
	
	//Get beam options
	GET_OPTION_VALUE(beamFWHM,m_beamFWHM);
	GET_OPTION_VALUE(pixSize,m_pixSize);
	GET_OPTION_VALUE(beamTheta,m_beamTheta);
	m_fluxCorrectionFactor= AstroUtils::GetBeamAreaInPixels(m_beamFWHM,m_beamFWHM,m_pixSize,m_pixSize);
	INFO_LOG("[PROC "<<m_procId<<"] - User-supplied beam info (fwhm="<<m_beamFWHM<<", dx="<<m_pixSize<<", theta="<<m_beamTheta<<", fluxCorrFactor="<<m_fluxCorrectionFactor<<")");

	//Get output file options
	GET_OPTION_VALUE(outputFile,m_OutputFileName);
	GET_OPTION_VALUE(saveToFile,m_saveToFile);
	GET_OPTION_VALUE(saveConfig,m_saveConfig);
	GET_OPTION_VALUE(saveDS9Region,m_saveDS9Region);
	GET_OPTION_VALUE(ds9RegionFile,m_DS9CatalogFileName);
	GET_OPTION_VALUE(DS9RegionFormat,m_DS9RegionFormat);
	GET_OPTION_VALUE(saveSources,m_saveSources);
	GET_OPTION_VALUE(isInteractiveRun,m_IsInteractiveRun);
	GET_OPTION_VALUE(saveResidualMap,m_saveResidualMap);
	GET_OPTION_VALUE(saveInputMap,m_saveInputMap);
	GET_OPTION_VALUE(saveSignificanceMap,m_saveSignificanceMap);
	GET_OPTION_VALUE(saveBkgMap,m_saveBkgMap);
	GET_OPTION_VALUE(saveNoiseMap,m_saveNoiseMap);
	GET_OPTION_VALUE(saveSaliencyMap,m_saveSaliencyMap);
	GET_OPTION_VALUE(saveEdgenessMap,m_saveEdgenessMap);
	GET_OPTION_VALUE(saveCurvatureMap,m_saveCurvatureMap);
	GET_OPTION_VALUE(saveSegmentedMap,m_saveSegmentedMap);
		
	//Get bkg options
	GET_OPTION_VALUE(useLocalBkg,m_UseLocalBkg);
	GET_OPTION_VALUE(use2ndPassInLocalBkg,m_Use2ndPassInLocalBkg);
	GET_OPTION_VALUE(skipOutliersInLocalBkg,m_SkipOutliersInLocalBkg);
	GET_OPTION_VALUE(localBkgMethod,m_LocalBkgMethod);
	GET_OPTION_VALUE(bkgEstimator,m_BkgEstimator);
	GET_OPTION_VALUE(boxSizeX,m_BoxSizeX);
	GET_OPTION_VALUE(boxSizeY,m_BoxSizeY);
	GET_OPTION_VALUE(gridSizeX,m_GridSizeX);
	GET_OPTION_VALUE(gridSizeY,m_GridSizeY);
	GET_OPTION_VALUE(useBeamInfoInBkg,m_UseBeamInfoInBkg);
	
	//Get source search options
	GET_OPTION_VALUE(searchCompactSources,m_SearchCompactSources);
	GET_OPTION_VALUE(minNPix,m_NMinPix);
	GET_OPTION_VALUE(seedBrightThr,m_SeedBrightThr);	
	GET_OPTION_VALUE(seedThr,m_SeedThr);
	GET_OPTION_VALUE(mergeThr,m_MergeThr);
	GET_OPTION_VALUE(compactSourceSearchNIters,m_compactSourceSearchNIters);
	GET_OPTION_VALUE(seedThrStep,m_seedThrStep);
	GET_OPTION_VALUE(mergeBelowSeed,m_MergeBelowSeed);
	GET_OPTION_VALUE(searchNegativeExcess,m_SearchNegativeExcess);
	
	//Get nested source search options
	GET_OPTION_VALUE(searchNestedSources,m_SearchNestedSources);
	GET_OPTION_VALUE(nestedBlobThrFactor,m_NestedBlobThrFactor);
	GET_OPTION_VALUE(minNestedMotherDist,m_minNestedMotherDist);
	GET_OPTION_VALUE(maxMatchingPixFraction,m_maxMatchingPixFraction);
	

	//Get source selection options
	GET_OPTION_VALUE(applySourceSelection,m_ApplySourceSelection);
	GET_OPTION_VALUE(sourceMinBoundingBox,m_SourceMinBoundingBox);	
	GET_OPTION_VALUE(useMinBoundingBoxCut,m_useMinBoundingBoxCut);
	GET_OPTION_VALUE(psCircRatioThr,m_psCircRatioThr);
	GET_OPTION_VALUE(useCircRatioCut,m_useCircRatioCut);
	GET_OPTION_VALUE(psElongThr,m_psElongThr);
	GET_OPTION_VALUE(useElongCut,m_useElongCut);
	GET_OPTION_VALUE(psEllipseAreaRatioMinThr,m_psEllipseAreaRatioMinThr);
	GET_OPTION_VALUE(psEllipseAreaRatioMaxThr,m_psEllipseAreaRatioMaxThr);
	GET_OPTION_VALUE(useEllipseAreaRatioCut,m_useEllipseAreaRatioCut);
	GET_OPTION_VALUE(psMaxNPix,m_psMaxNPix);
	GET_OPTION_VALUE(useMaxNPixCut,m_useMaxNPixCut);
	
	//Get source residual options
	GET_OPTION_VALUE(dilateNestedSources,m_DilateNestedSources);
	GET_OPTION_VALUE(dilateKernelSize,m_DilateKernelSize);
	GET_OPTION_VALUE(dilatedSourceType,m_DilatedSourceType);
	GET_OPTION_VALUE(dilateSourceModel,m_DilateSourceModel);
	GET_OPTION_VALUE(dilateRandomize,m_DilateRandomize);
	
	//Get source fitting options
	GET_OPTION_VALUE(fitSources,m_fitSources);
	GET_OPTION_VALUE(fitMaxNComponents,m_fitMaxNComponents);
	//GET_OPTION_VALUE(deblendCurvThr,m_deblendCurvThr);
	//GET_OPTION_VALUE(deblendComponentMinNPix,m_deblendComponentMinNPix);
	
	//Get smoothing options
	GET_OPTION_VALUE(usePreSmoothing,m_UsePreSmoothing);
	GET_OPTION_VALUE(smoothFilter,m_SmoothFilter);
	GET_OPTION_VALUE(gausFilterKernSize,m_GausFilterKernSize);
	GET_OPTION_VALUE(gausFilterSigma,m_GausFilterSigma);
	GET_OPTION_VALUE(guidedFilterRadius,m_GuidedFilterRadius);
	GET_OPTION_VALUE(guidedFilterColorEps,m_GuidedFilterColorEps);
	
	//Get saliency options
	GET_OPTION_VALUE(saliencyThrFactor,m_SaliencyThrFactor);
	GET_OPTION_VALUE(saliencyBkgThrFactor,m_SaliencyBkgThrFactor);
	GET_OPTION_VALUE(saliencyImgThrFactor,m_SaliencyImgThrFactor);
	GET_OPTION_VALUE(saliencyResoMin,m_SaliencyResoMin);
	GET_OPTION_VALUE(saliencyResoMax,m_SaliencyResoMax);
	GET_OPTION_VALUE(saliencyResoStep,m_SaliencyResoStep);
	GET_OPTION_VALUE(saliencyUseRobustPars,m_SaliencyUseRobustPars);
	GET_OPTION_VALUE(saliencyUseBkgMap,m_SaliencyUseBkgMap);
	GET_OPTION_VALUE(saliencyUseNoiseMap,m_SaliencyUseNoiseMap);
	GET_OPTION_VALUE(saliencyUseCurvInDiss,m_SaliencyUseCurvInDiss);
	GET_OPTION_VALUE(saliencyNNFactor,m_SaliencyNNFactor);
	GET_OPTION_VALUE(saliencySpatialRegFactor,m_SaliencySpatialRegFactor);
	GET_OPTION_VALUE(saliencyMultiResoCombThrFactor,m_SaliencyMultiResoCombThrFactor);
	GET_OPTION_VALUE(saliencyDissExpFalloffPar,m_SaliencyDissExpFalloffPar);
	GET_OPTION_VALUE(saliencySpatialDistRegPar,m_SaliencySpatialDistRegPar);
		
	//Get extended source options
	GET_OPTION_VALUE(searchExtendedSources,m_SearchExtendedSources);
	GET_OPTION_VALUE(extendedSearchMethod,m_ExtendedSearchMethod);
	GET_OPTION_VALUE(wtScaleExtended,m_wtScaleExtended);
	GET_OPTION_VALUE(useResidualInExtendedSearch,m_UseResidualInExtendedSearch);
	GET_OPTION_VALUE(activeContourMethod,m_activeContourMethod);
			
	//Get superpixel options
	GET_OPTION_VALUE(spSize,m_spSize);
	GET_OPTION_VALUE(spBeta,m_spBeta);
	GET_OPTION_VALUE(spMinArea,m_spMinArea);
	GET_OPTION_VALUE(spUseLogContrast,m_spUseLogContrast);

	//Chan-Vese options
	GET_OPTION_VALUE(cvNIters,m_cvNIters);
	GET_OPTION_VALUE(cvTimeStepPar,m_cvTimeStepPar);
	GET_OPTION_VALUE(cvWindowSizePar,m_cvWindowSizePar);
	GET_OPTION_VALUE(cvLambda1Par,m_cvLambda1Par);
	GET_OPTION_VALUE(cvLambda2Par,m_cvLambda2Par);
	GET_OPTION_VALUE(cvMuPar,m_cvMuPar);
	GET_OPTION_VALUE(cvNuPar,m_cvNuPar);
	GET_OPTION_VALUE(cvPPar,m_cvPPar);

	//LRAC algorithm options
	GET_OPTION_VALUE(lracNIters,m_lracNIters);	
	GET_OPTION_VALUE(lracLambdaPar,m_lracLambdaPar);
	GET_OPTION_VALUE(lracRadiusPar,m_lracRadiusPar);
	GET_OPTION_VALUE(lracEpsPar,m_lracEpsPar);

	//Hierarchical clustering options
	GET_OPTION_VALUE(spMergingEdgeModel,m_spMergingEdgeModel);
	GET_OPTION_VALUE(spMergingRegPar,m_spMergingRegPar);
	GET_OPTION_VALUE(spMergingNSegmentsToStop,m_spMergingNSegmentsToStop);
	GET_OPTION_VALUE(spMergingRatio,m_spMergingRatio);
	GET_OPTION_VALUE(spMergingMaxDissRatio,m_spMergingMaxDissRatio);
	GET_OPTION_VALUE(spMergingMaxDissRatio2ndNeighbours,m_spMergingMaxDissRatio2ndNeighbours);
	GET_OPTION_VALUE(spMergingDissThreshold,m_spMergingDissThreshold);
	GET_OPTION_VALUE(spMergingIncludeSpatialPars,m_spMergingIncludeSpatialPars);
	GET_OPTION_VALUE(spMergingUseRobustPars,m_spMergingUseRobustPars);
	GET_OPTION_VALUE(spMergingAddCurvDist,m_spMergingAddCurvDist);
	
	return 0;

}//close Configure()


int SFinder::RunTask(TaskData* taskData,bool storeData){

	//Check task data
	if(!taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to task data given!");
		return -1;
	}

	//Get task data
	long int ix_min= taskData->ix_min;
	long int ix_max= taskData->ix_max;
	long int iy_min= taskData->iy_min;
	long int iy_max= taskData->iy_max;

	//Define task data
	Image* taskImg= 0;
	ImgBkgData* bkgData= 0;
	Image* significanceMap= 0;
	Image* segmentedImg= 0;
	bool stopTask= false;
	int status= 0;

	//==================================
	//==   Read task input image
	//==================================
	INFO_LOG("[PROC "<<m_procId<<"] - Reading input image ["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"]...");
	auto t0_read = chrono::steady_clock::now();	
	
	FileInfo info;
	taskImg= ReadImage(info,m_InputFileName,m_InputImgName,ix_min,ix_max,iy_min,iy_max);
	
	// Set image physical boundary in task data
	if(taskImg){
		double xmin= taskImg->GetXmin();
		double xmax= taskImg->GetXmax();
		double ymin= taskImg->GetYmin();
		double ymax= taskImg->GetYmax();
		//taskData->x_min= xmin;
		//taskData->x_max= xmax;
		//taskData->y_min= ymin;
		//taskData->y_max= ymax;
	}	
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Reading of input image failed, skip to next task...");
		stopTask= true;
		status= -1;
	}
	auto t1_read = chrono::steady_clock::now();	
	readImageTime+= chrono::duration <double, milli> (t1_read-t0_read).count();
	

	//============================
	//== Find image stats & bkg
	//============================
	if(!stopTask){
		INFO_LOG("[PROC "<<m_procId<<"] - Computing image bkg...");
		bkgData= ComputeStatsAndBkg(taskImg);
		if(!bkgData){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute bkg for input image!");
			stopTask= true;
			status= -1;
		}
	}

	//============================
	//== Find compact sources
	//============================	
	if(!stopTask && m_SearchCompactSources ){ 	
		INFO_LOG("[PROC "<<m_procId<<"] - Searching compact sources...");
		auto t0_sfinder = chrono::steady_clock::now();	

		//significanceMap= FindCompactSources(taskImg,bkgData,taskData);
		significanceMap= FindCompactSourcesRobust(taskImg,bkgData,taskData,m_compactSourceSearchNIters);
		if(!significanceMap){
			ERROR_LOG("[PROC "<<m_procId<<"] - Compact source search failed!");
			//stopTask= true;
			status= -1;
		}
		auto t1_sfinder = chrono::steady_clock::now();	
		compactSourceTime+= chrono::duration <double, milli> (t1_sfinder-t0_sfinder).count();
	}

	//============================
	//== Find extended sources
	//============================
	if(!stopTask && m_SearchExtendedSources){
		INFO_LOG("[PROC "<<m_procId<<"] - Searching extended sources...");
		auto t0_extsfinder = chrono::steady_clock::now();	

		segmentedImg= FindExtendedSources(taskImg,bkgData,taskData,storeData);
		if(!segmentedImg){
			ERROR_LOG("[PROC "<<m_procId<<"] - Extended source search failed!");
			status= -1;
		}
		auto t1_extsfinder = chrono::steady_clock::now();	
		extendedSourceTime+= chrono::duration <double, milli> (t1_extsfinder-t0_extsfinder).count();

	}//close if search extended sources

	//============================
	//== Find edge sources
	//============================
	if(!stopTask){
		INFO_LOG("[PROC "<<m_procId<<"] - Finding sources at tile edges ...");
		if(FindSourcesAtEdge()<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Finding sources at tile edges failed!");
			status= -1;
		}
	}

	//============================
	//== Store & Clear data
	//============================
	//## Store data?
	if(storeData){
		if(taskImg) m_InputImg= taskImg;
		if(bkgData) m_BkgData= bkgData;
		if(significanceMap) m_SignificanceMap= significanceMap;
		if(segmentedImg) m_SegmImg= segmentedImg;
	}
	else{//clear all data
		if(taskImg){
			delete taskImg;
			taskImg= 0;
		}
		if(bkgData){
			delete bkgData;
			bkgData= 0;	
		}
		if(significanceMap){
			delete significanceMap;
			significanceMap= 0;
		}
		if(segmentedImg){
			delete segmentedImg;
			segmentedImg= 0;
		}
	}//close else

	return status;

}//close RunTask()



int SFinder::Run(){

	//Start timer
	auto t0 = chrono::steady_clock::now();

	//================================================
	//== Init options & data (done by all processors)
	//================================================
	INFO_LOG("[PROC "<<m_procId<<"] - Initializing source finder...");
	auto t0_init = chrono::steady_clock::now();
	if(Init()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Initialization failed!");
		return -1;
	}
	auto t1_init = chrono::steady_clock::now();
	initTime= chrono::duration <double, milli> (t1_init-t0_init).count();

	
	//================================================
	//== DISTRIBUTE TASKS TO WORKERS 
	//================================================
	//## Start loop on tasks per worker
	int status= 0;
	bool storeData= true;
	//if(m_mpiEnabled) storeData= false;
	if(m_taskDataPerWorkers.size()>1) storeData= false;
	size_t nTasks= m_taskDataPerWorkers[m_procId].size(); 
	INFO_LOG("[PROC "<<m_procId<<"] - Start processing of #"<<nTasks<<" tasks...");
	
	for(size_t j=0;j<nTasks;j++){

		//Run task
		if(RunTask(m_taskDataPerWorkers[m_procId][j],storeData)<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to run task no. "<<j<<", skip to next!");
			status= -1;
			continue;
		}

	}//end loop tasks per worker
	
	if(status<0){
		WARN_LOG("One or more errors occurred in source finding tasks...");
	}

	
	//## Update task info (tile physical range) from workers
	//## The updated list of task data is available in master processor
	//## (check if it is better to replace with MPI_Gather and put it available in all workers)
	#ifdef MPI_ENABLED
	if(m_mpiEnabled){
		INFO_LOG("Gathering task data from all workers...");
		if(GatherTaskDataFromWorkers()<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Gathering task data from workers failed!");
			return -1;
		}
	}
	#endif

	
	/*
	//## Find edge sources
	#ifdef MPI_ENABLED
		if(m_mpiEnabled) MPI_Barrier(MPI_COMM_WORLD);
	#endif
	if(m_procId==MASTER_ID) {
		INFO_LOG("[PROC "<<m_procId<<"] - Finding sources at tile edges ...");
		if(FindSourcesAtEdge()<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Finding sources at tile edges failed!");
			return -1;
		}
	}
	*/
	
	//## Merge sources at edge
	if(m_mergeSourcesAtEdge){
		#ifdef MPI_ENABLED
			if(m_mpiEnabled) MPI_Barrier(MPI_COMM_WORLD);
		#endif
		if(m_procId==MASTER_ID) {
			INFO_LOG("[PROC "<<m_procId<<"] - Merging sources find at tile edges by all workers...");
			if(MergeSourcesAtEdge()<0){
				ERROR_LOG("[PROC "<<m_procId<<"] - Merging sources at tile edges failed!");
				return -1;
			}
		}
	}//close if merge sources at edges

	//## Merge sources found in each task in unique collection
	if(m_procId==MASTER_ID && MergeTaskData()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Merging sources found in each task in unique collection failed!");
		return -1;
	}


	//============================
	//== Fit sources
	//============================
	if(m_fitSources && m_procId==MASTER_ID) {
		auto t0_sfit = chrono::steady_clock::now();	
		if(FitSources(m_SourceCollection)<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to fit sources!");
			return -1;
		}
		auto t1_sfit = chrono::steady_clock::now();	
		sourceFitTime= chrono::duration <double, milli> (t1_sfit-t0_sfit).count();
	}

	//Stop timer
	auto t1 = chrono::steady_clock::now();	
	totTime= chrono::duration <double, milli> (t1-t0).count();

	//============================
	//== Save to file
	//============================
	if(m_saveToFile && m_procId==MASTER_ID) {
		INFO_LOG("[PROC "<<m_procId<<"] - Saving data to files...");	
		auto t0_save = chrono::steady_clock::now();
		Save();	
		auto t1_save = chrono::steady_clock::now();	
		saveTime= chrono::duration <double, milli> (t1_save-t0_save).count();
	}
	#ifdef MPI_ENABLED
	if(m_mpiEnabled) MPI_Barrier(MPI_COMM_WORLD);
	#endif	

	//===============================
	//== Print performance stats
	//===============================
	if(m_procId==MASTER_ID) PrintPerformanceStats();


	return 0;

}//close Run()


/*
int SFinder::Run(){

	//Start timer
	auto t0 = chrono::steady_clock::now();	

	//===========================
	//== Init options & data
	//===========================
	INFO_LOG("Initializing source finder...");
	auto t0_init = chrono::steady_clock::now();	
	if(Init()<0){
		ERROR_LOG("Initialization failed!");
		return -1;
	}
	auto t1_init = chrono::steady_clock::now();	
	initTime= chrono::duration <double, milli> (t1_init-t0_init).count();
	
	//===========================
	//== Read image
	//===========================
	INFO_LOG("Reading input image...");
	auto t0_read = chrono::steady_clock::now();

	FileInfo info;
	if(m_ReadTile) m_InputImg= ReadImage(info,m_InputFileName,m_InputImgName,m_TileMinX,m_TileMaxX,m_TileMinY,m_TileMaxY);	
	else m_InputImg= ReadImage(info,m_InputFileName,m_InputImgName);
	
	//if(ReadImage()<0){
	if(!m_InputImg){
		ERROR_LOG("Reading of input image failed!");
		return -1;
	}
	m_InputFileExtension= info.extension;

	auto t1_read = chrono::steady_clock::now();	
	readImageTime= chrono::duration <double, milli> (t1_read-t0_read).count();

	//============================
	//== Find compact sources
	//============================
	INFO_LOG("Searching compact sources...");
	auto t0_sfinder = chrono::steady_clock::now();	
	if(m_SearchCompactSources && FindCompactSources()<0){
		ERROR_LOG("Compact source search failed!");
		return -1;
	}
	auto t1_sfinder = chrono::steady_clock::now();	
	compactSourceTime= chrono::duration <double, milli> (t1_sfinder-t0_sfinder).count();


	//============================
	//== Find extended sources
	//============================
	if(m_SearchExtendedSources){
		INFO_LOG("Searching extended sources...");

		// Find residual map
		INFO_LOG("Computing residual image ...");
		auto t0_res = chrono::steady_clock::now();	
		if(FindResidualMap()<0){
			ERROR_LOG("Residual map computation failed!");
			return -1;
		}
		auto t1_res = chrono::steady_clock::now();	
		imgResidualTime= chrono::duration <double, milli> (t1_res-t0_res).count();

		
		//Find extended sources
		auto t0_extsfinder = chrono::steady_clock::now();
		if(FindExtendedSources(m_ResidualImg)<0){
			ERROR_LOG("Extended source search failed!");
			return -1;
		}
		auto t1_extsfinder = chrono::steady_clock::now();	
		extendedSourceTime= chrono::duration <double, milli> (t1_extsfinder-t0_extsfinder).count();

	}//close if search extended sources

	
	//============================
	//== Fit sources
	//============================
	if(m_fitSources) {
		auto t0_sfit = chrono::steady_clock::now();	
		if(FitSources(m_SourceCollection)<0){
			ERROR_LOG("Failed to fit sources!");
			return -1;
		}
		auto t1_sfit = chrono::steady_clock::now();	
		sourceFitTime= chrono::duration <double, milli> (t1_sfit-t0_sfit).count();
	}

	//============================
	//== Draw images
	//============================
	if(m_IsInteractiveRun) {
		DrawSources(m_InputImg,m_SourceCollection);
	}

	//Stop timer
	auto t1 = chrono::steady_clock::now();	
	totTime= chrono::duration <double, milli> (t1-t0).count();

	//============================
	//== Save to file
	//============================
	if(m_saveToFile) {	
		auto t0_save = chrono::steady_clock::now();	
		Save();	
		auto t1_save = chrono::steady_clock::now();	
		saveTime= chrono::duration <double, milli> (t1_save-t0_save).count();
	}

	//===============================
	//== Print performance stats
	//===============================
	PrintPerformanceStats();

	//==========================================
	//== Run TApplication (interactive run)
	//==========================================
	if(m_Application && m_IsInteractiveRun) {
		m_Application->Run();
	}
	
	return 0;

}//close Run()
*/


int SFinder::FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr,Image* searchedImg){

	if(!inputImg) {
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to given input image!");
		return -1;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	//## Compute stats and bkg
	INFO_LOG("[PROC "<<m_procId<<"] - Computing image stats/bkg...");
	ImgBkgData* bkgData= ComputeStatsAndBkg(img);	
	if(!bkgData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	Image* significanceMap= img->GetSignificanceMap(bkgData,m_UseLocalBkg);
	if(!significanceMap){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute significance map!");
		return -1;
	}

	//## Find sources
	INFO_LOG("[PROC "<<m_procId<<"] - Finding sources...");	
	int status= inputImg->FindCompactSource(
		sources,
		significanceMap,bkgData,
		seedThr,mergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
		m_SearchNestedSources,m_NestedBlobThrFactor, m_minNestedMotherDist, m_maxMatchingPixFraction
	);

	//## Clear data
	if(bkgData){
		delete bkgData;
		bkgData= 0;
	}
	if(significanceMap) {
		delete significanceMap;	
		significanceMap= 0;
	}

	//## Check status
	if(status<0) {
		ERROR_LOG("[PROC "<<m_procId<<"] - Source finding failed!");	
		return -1;
	}
	int nSources= static_cast<int>(sources.size());
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" sources detected in input image...");	
	
	return 0;

}//close FindSources()

Image* SFinder::FindCompactSourcesRobust(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,int niter)
{
	//## Check img
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input img and/or bkg/task data!");
		return nullptr;
	}

	//## Copy input image (this will be masked with sources at each iteration)
	bool copyMetaData= true;
	bool resetStats= false;
	Image* img= inputImg->GetCloned("",copyMetaData,resetStats);
	if(!img){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to clone input image!");
		return nullptr;
	}

	//## Compute curvature map to be added to sources
	Image* curvMap= ComputeLaplacianImage(inputImg);
	if(!curvMap){
		WARN_LOG("[PROC "<<m_procId<<"] - Failed to compute curvature map!");
		return nullptr;
	}

	//## Iterate source finding subtracting sources found at each iteration
	std::vector<Source*> sources;
	Image* significanceMap= 0;
	int iterations_done= 0;
	double seedThr= m_SeedThr;
	
	for(int k=0;k<niter;k++){
	
		//## Compute stats & bkg at current iteration
		ImgBkgData* bkgData_iter= 0;
		if(k==0){//Copy bkg data
			seedThr= m_SeedThr;

			bkgData_iter= new ImgBkgData;
			*bkgData_iter= *bkgData;
		}
		else{//Compute bkg data
			double seedThr_iter= seedThr-m_seedThrStep;//reduce seedThr by step after the first iteration
			if(seedThr_iter>m_MergeThr) seedThr= seedThr_iter;//do not allow seed thr equal or smaller than mergeThr (if so stop seedThr adaptive decrease)

			INFO_LOG("[PROC "<<m_procId<<"] - Computing image stats & bkg at iter no. "<<k+1<<" ...");
			bkgData_iter= ComputeStatsAndBkg(img);
			if(!bkgData_iter){
				ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute bkg for input image at iter no. "<<k+1<<"!");
				delete img;
				img= 0;
				for(size_t i=0;i<sources.size();i++) {	
					if(sources[i]) {
						delete sources[i];
						sources[i]= 0;
					}
				}
				sources.clear();
				return nullptr;
			}
		}//close else

		//## Compute significance map
		Image* significanceMap_iter= img->GetSignificanceMap(bkgData_iter,m_UseLocalBkg);
		if(!significanceMap_iter){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute significance map at iter no. "<<k+1<<"!");
			delete img;
			img= 0;
			delete bkgData_iter;
			bkgData_iter= 0;
			for(size_t i=0;i<sources.size();i++) {	
				if(sources[i]) {
					delete sources[i];
					sources[i]= 0;
				}
			}
			sources.clear();
			delete curvMap;
			curvMap= 0;
			return nullptr;
		}

		//## Find compact sources
		INFO_LOG("[PROC "<<m_procId<<"] - Finding compact sources at iter no. "<<k+1<<" ...");
		std::vector<Source*> sources_iter;
		int status= img->FindCompactSource(
			sources_iter,
			significanceMap_iter,bkgData_iter,
			seedThr,m_MergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
			m_SearchNestedSources,m_NestedBlobThrFactor, m_minNestedMotherDist, m_maxMatchingPixFraction,
			curvMap
		);

		if(status<0) {
			ERROR_LOG("[PROC "<<m_procId<<"] - Compact source finding failed!");
			delete img;
			img= 0;
			delete bkgData_iter;
			bkgData_iter= 0;
			delete significanceMap_iter;
			significanceMap_iter= 0;
			for(size_t i=0;i<sources.size();i++) {	
				if(sources[i]) {
					delete sources[i];
					sources[i]= 0;
				}
			}
			delete curvMap;
			curvMap= 0;
			return nullptr;
		}

		//## If no sources have been found stop iteration loop
		if(sources_iter.empty()){
			INFO_LOG("[PROC "<<m_procId<<"] - No compact sources found at iter "<<k+1<<", stop iteration.");
			delete bkgData_iter;
			bkgData_iter= 0;
			significanceMap= significanceMap_iter;
			break;
		}
		INFO_LOG("[PROC "<<m_procId<<"] - #"<<sources_iter.size()<<" compact sources found at iter "<<k+1<<", appending them to list...");
		iterations_done++;

		//## Append sources found at this iteration to main collection
		sources.insert(sources.end(),sources_iter.begin(),sources_iter.end());

		//## Mask sources found (replace pixels with zeros) at this iteration
		INFO_LOG("[PROC "<<m_procId<<"] - Masking #"<<sources_iter.size()<<" sources found at iter "<<k+1<<"...");
		img->MaskSources(sources_iter,0.);		

		//## Clear iter data
		delete bkgData_iter;
		bkgData_iter= 0;
		delete significanceMap;
		significanceMap= 0;	

	}//end loop iterations

	//Clearup data
	if(img){
		delete img;
		img= 0;
	}
	if(curvMap){
		delete curvMap;
		curvMap= 0;
	}

	//## Tag found sources as compact 
	int nSources= static_cast<int>( sources.size() );
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" compact sources detected in input image after #"<<iterations_done<<" iterations ...");
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetType(Source::eCompact);
	}
	
	//## Apply source selection?
	if(m_ApplySourceSelection && nSources>0){
		if(SelectSources(sources)<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to select sources!");
			delete significanceMap;
			significanceMap= 0;
			return nullptr;
		}
		nSources= static_cast<int>(sources.size());
	}//close if source selection


	//## Set flux correction factor
	double fluxCorrection= 1;
	bool hasBeamData= false;
	if(inputImg->HasMetaData()){
		fluxCorrection= inputImg->GetMetaData()->GetBeamFluxIntegral();
		if(fluxCorrection>0 && std::isnormal(fluxCorrection)) hasBeamData= true;
	}

	if(!hasBeamData){
		INFO_LOG("[PROC "<<m_procId<<"] - Beam information are not available in image or invalid, using correction factor ("<<m_fluxCorrectionFactor<<") computed from user-supplied beam info ...");	
		fluxCorrection= m_fluxCorrectionFactor;
	}
	
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetName(Form("S%d",(signed)k));
		sources[k]->SetBeamFluxIntegral(fluxCorrection);
		sources[k]->Print();
	}//end loop sources
			
	//## Add sources to task data sources
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());			
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" compact sources added to task data ...");

	return significanceMap;

}//close FindCompactSourcesRobust()


Image* SFinder::FindCompactSources(Image* inputImg, ImgBkgData* bkgData, TaskData* taskData){

	//## Check img
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input img and/or bkg/task data!");
		return nullptr;
	}

	//## Compute significance map
	Image* significanceMap= inputImg->GetSignificanceMap(bkgData,m_UseLocalBkg);
	if(!significanceMap){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute significance map!");
		return nullptr;
	}

	//## Find sources
	INFO_LOG("[PROC "<<m_procId<<"] - Finding compact sources...");	
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(
		sources,
		significanceMap,bkgData,
		m_SeedThr,m_MergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
		m_SearchNestedSources,m_NestedBlobThrFactor,m_minNestedMotherDist,m_maxMatchingPixFraction
	);

	if(status<0) {
		ERROR_LOG("[PROC "<<m_procId<<"] - Compact source finding failed!");
		delete significanceMap;
		significanceMap= 0;
		return nullptr;
	}

	//## Tag found sources as compact 
	int nSources= static_cast<int>( sources.size() );
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" compact sources detected in input image ...");
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetType(Source::eCompact);
	}
	
	//## Apply source selection?
	//int nSources= static_cast<int>(taskData->sources.size());
	
	if(m_ApplySourceSelection && nSources>0){
		if(SelectSources(sources)<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to select sources!");
			delete significanceMap;
			significanceMap= 0;
			return nullptr;
		}
		nSources= static_cast<int>(sources.size());
	}//close if source selection


	//## Set flux correction factor
	double fluxCorrection= 1;
	bool hasBeamData= false;
	if(inputImg->HasMetaData()){
		fluxCorrection= inputImg->GetMetaData()->GetBeamFluxIntegral();
		if(fluxCorrection>0) hasBeamData= true;
	}

	if(!hasBeamData){
		INFO_LOG("Beam information are not available in image or invalid, using correction factor ("<<m_fluxCorrectionFactor<<") computed from user-supplied beam info ...");	
		fluxCorrection= m_fluxCorrectionFactor;
	}
	
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetName(Form("S%d",(signed)k));
		sources[k]->SetBeamFluxIntegral(fluxCorrection);
		sources[k]->Print();
	}	
			
	//## Add sources to task data sources
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());			
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" compact sources added to task data ...");

	return significanceMap;

}//close FindCompactSources()



Image* SFinder::FindResidualMap(Image* inputImg,ImgBkgData* bkgData,std::vector<Source*> const & sources){

	//Check input image
	if(!inputImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to given input image!");
		return nullptr;
	}

	//Compute residual map
	Image* residualImg= 0;
	if(m_UseResidualInExtendedSearch && sources.size()>0){
		INFO_LOG("[PROC "<<m_procId<<"] - Computing residual image (#"<<sources.size()<<" sources present)...");
		residualImg= inputImg->GetSourceResidual(
			sources,
			m_DilateKernelSize,m_DilateSourceModel,m_DilatedSourceType,m_DilateNestedSources,	
			bkgData,m_UseLocalBkg,
			m_DilateRandomize,m_SeedBrightThr
		);
	}//close if
	else{
		INFO_LOG("[PROC "<<m_procId<<"] - Setting residual image to input image...");
		residualImg= inputImg->GetCloned("",true,true);
	}

	if(!residualImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute residual map!");
		return nullptr;
	}
	
	return residualImg;

}//close FindResidualMap()



Image* SFinder::ComputeSmoothedImage(Image* inputImg,int model){

	//Check input image
	if(!inputImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to given input image!");
		return nullptr;
	}

	//Compute smoothed image
	Image* smoothedImg= 0;	
	if(model==eGausFilter){
		smoothedImg= inputImg->GetSmoothedImage(m_GausFilterKernSize,m_GausFilterKernSize,m_GausFilterSigma,m_GausFilterSigma);
	}
	else if(model==eGuidedFilter){
		smoothedImg= inputImg->GetGuidedFilterImage(m_GuidedFilterRadius,m_GuidedFilterColorEps);
	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid smoothing algo ("<<model<<") selected!");
		return nullptr;
	}

	if(!smoothedImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Image smoothing failed!");
		return nullptr;
	}

	return smoothedImg;

}//close ComputeSmoothedImage()


Image* SFinder::FindExtendedSources(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,bool storeData)
{
	//## Check input image
	if(!inputImg || !taskData || !bkgData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input image and/or bkg/task data given!");
		return nullptr;
	}

	//## Check segmentation method
	std::vector<int> segmMethods= {eHClust,eActiveContour,eWaveletTransform,eSaliencyThr};
	bool foundMethod= false;
	for(size_t i=0;i<segmMethods.size();i++){
		if(m_ExtendedSearchMethod==segmMethods[i]){
			foundMethod= true;
			break;
		}
	}
	if(!foundMethod){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid extended source method selected (method="<<m_ExtendedSearchMethod<<")!");
		return nullptr;
	}

	Image* searchImg= 0;
	
	//****************************
	//** Find residual map
	//****************************
	INFO_LOG("[PROC "<<m_procId<<"] - Computing residual image ...");
	Image* residualImg= FindResidualMap(inputImg,bkgData,taskData->sources);
	if(!residualImg){
		ERROR_LOG("Residual map computation failed!");
		return nullptr;
	}
	if(storeData) m_ResidualImg= residualImg;
	searchImg= residualImg;

	
	//Compute bkg & noise map for residual img
	INFO_LOG("[PROC "<<m_procId<<"] - Computing residual image stats & bkg...");
	ImgBkgData* residualBkgData= ComputeStatsAndBkg(residualImg);
	if(!residualBkgData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute bkg data for residual map!");
		if(residualImg && !storeData){
			delete residualImg;
			residualImg= 0;
		}
		return nullptr;
	}
	if(storeData) m_ResidualBkgData= residualBkgData;


	//****************************
	//** Smooth image
	//****************************
	Image* smoothedImg= 0;
	if(m_UsePreSmoothing){
		smoothedImg= ComputeSmoothedImage(residualImg,m_SmoothFilter);
		if(!smoothedImg){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute residual smoothed image!");
			if(residualImg && !storeData){
				delete residualImg;
				residualImg= 0;
			}
			if(residualBkgData && !storeData){
				delete residualBkgData;
				residualBkgData= 0;
			}
			return nullptr;
		}

		//Set extended source search image to smoothed residual image
		searchImg	= smoothedImg;
	}
	

	//****************************
	//** Run segmentation
	//****************************
	INFO_LOG("[PROC "<<m_procId<<"] - Run extended source segmentation algorithm...");
	Image* segmentedImg= 0;
	if(m_ExtendedSearchMethod==eHClust){
		segmentedImg= FindExtendedSources_HClust(inputImg,residualBkgData,taskData,searchImg,storeData);
	}
	else if(m_ExtendedSearchMethod==eActiveContour){
		segmentedImg= FindExtendedSources_AC(inputImg,residualBkgData,taskData,searchImg,storeData);
	}
	else if(m_ExtendedSearchMethod==eWaveletTransform){
		segmentedImg= FindExtendedSources_WT(inputImg,taskData,searchImg);
	}
	else if(m_ExtendedSearchMethod==eSaliencyThr){
		segmentedImg= FindExtendedSources_SalThr(inputImg,residualBkgData,taskData,searchImg,storeData);
	}

	//Check if segmentation succeeded
	if(!segmentedImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to run the segmentation algorithm!");	
		if(residualImg && !storeData){
			delete residualImg;
			residualImg= 0;
		}
		if(residualBkgData && !storeData){
			delete residualBkgData;
			residualBkgData= 0;
		}
		if(smoothedImg){
			delete smoothedImg;
			smoothedImg= 0;
		}
	}


	//Clear data
	if(residualImg && !storeData) {
		delete residualImg;
		residualImg= 0;
	}
	if(residualBkgData && !storeData) {
		delete residualBkgData;
		residualBkgData= 0;
	}
	if(smoothedImg && !storeData) {
		delete smoothedImg;		
		smoothedImg= 0;
	}
	
	return segmentedImg;

}//close FindExtendedSources()



Image* SFinder::FindExtendedSources_SalThr(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg,bool storeData){

	//Check input image
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input image and/or bkg/task data given!");
		return nullptr;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	Image* saliencyImg= img->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,bkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(saliencyImg){
		if(storeData) m_SaliencyImg= saliencyImg;	
	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute saliency map!");
		return nullptr;
	}

	//## Get saliency map optimal threshold
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= saliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	
	//==========================================
	//==    FIND SOURCES
	//==========================================
	//## Find compact blobs in saliency map by simple thresholding
	INFO_LOG("[PROC "<<m_procId<<"] - Finding blobs in saliency map with threshold="<<signalThr<<"...");	
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	int minNPix= m_NMinPix;
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(	
		sources, saliencyImg, 0,
		signalThr,signalThr,minNPix,
		findNegativeExcess,mergeBelowSeed,findNestedSources
	);

	//Delete saliency map if not needed
	if(saliencyImg && !storeData){
		delete saliencyImg;
		saliencyImg= 0;
	}

	if(status<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Compact source finding with saliency map failed!");
		return nullptr;
	}

	//## Tag found sources as extended 
	int nSources= static_cast<int>( sources.size() );
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" extended sources detected in input image by thresholding the saliency map...");
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetName(Form("Sext%d",(signed)k));
		sources[k]->SetType(Source::eExtended);
	}
	
	//## Add sources to extended sources?
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());	
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	
	//## Compute segmented map
	bool isBinary= true;
	bool invert= false;
	Image* segmMap= inputImg->GetSourceMask(sources,isBinary,invert);
	if(!segmMap){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute segmented map!");
		return nullptr;
	}
	
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" extended sources to the list...");

	return segmMap;

}//close FindExtendedSources_SalThr()



Image* SFinder::FindExtendedSources_HClust(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg,bool storeData){

	//Check input image
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input image and/or bkg/task data given!");
		return nullptr;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	Image* saliencyImg= inputImg->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,bkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(saliencyImg){
		if(storeData) m_SaliencyImg= saliencyImg;
 	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute saliency map!");
		return nullptr;
	}

	//## Threshold saliency map and get signal and bkg markers
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= saliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	double bkgThr= saliencyImg->FindMedianThreshold(m_SaliencyBkgThrFactor);
	if(TMath::IsNaN(signalThr) || fabs(signalThr)==TMath::Infinity()){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid numeric threshold returned as threshold computation failed!");
		return nullptr;
	}

	INFO_LOG("[PROC "<<m_procId<<"] - Computing binarized saliency maps (signalThr="<<signalThr<<", bkgThr="<<bkgThr);
	double fgValue= 1;
	Image* signalMarkerImg= saliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	Image* bkgMarkerImg= saliencyImg->GetBinarizedImage(bkgThr,fgValue,true);
	
	//## Delete saliency map if not needed
	if(saliencyImg && !storeData){
		delete saliencyImg;
		saliencyImg= 0;
	}

	//## Compute Laplacian filtered image
	INFO_LOG("[PROC "<<m_procId<<"] - Computing laplacian image...");
	Image* laplImg= ComputeLaplacianImage(img);
	if(m_LaplImg){
		if(storeData) m_LaplImg= laplImg;
	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute laplacian image, cannot perform extended source finding!");
		return nullptr;
	}

	//## Compute edge image	
	INFO_LOG("[PROC "<<m_procId<<"] - Computing edgeness image...");
	Image* edgeImg= ComputeEdgeImage(img,m_spMergingEdgeModel);
	if(edgeImg){
		if(storeData) m_EdgeImg= edgeImg;
	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute the edgeness image, cannot perform extended source finding!");
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		return nullptr;
	}

	//## Compute the Superpixel partition
	bool normalizeImage= true;
	SLICData* slicData_init= SLIC::SPGenerator(img,m_spSize,m_spBeta,m_spMinArea,normalizeImage,m_spUseLogContrast,laplImg,edgeImg);
	if(!slicData_init){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute the initial superpixel partition, cannot perform extended source finding!");	
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		if(edgeImg && !storeData){
			delete edgeImg;		
			edgeImg= 0;
		}
		return nullptr;
	}

	//## Tag the superpixel partition
	if(SLIC::TagRegions(slicData_init->regions,bkgMarkerImg,signalMarkerImg)<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to tag (signal vs bkg) the initial superpixel partition, cannot perform extended source finding!");
		delete slicData_init;
		slicData_init= 0;
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		if(edgeImg && !storeData){
			delete edgeImg;		
			edgeImg= 0;
		}
		if(bkgMarkerImg){
			delete bkgMarkerImg;
			bkgMarkerImg= 0;
		}	
		if(signalMarkerImg){
			delete signalMarkerImg;
			signalMarkerImg= 0;
		}
		return nullptr;
	}

	//## Delete signal/bkg marker images
	if(bkgMarkerImg){
		delete bkgMarkerImg;
		bkgMarkerImg= 0;
	}	
	if(signalMarkerImg){
		delete signalMarkerImg;
		signalMarkerImg= 0;
	}

	//==========================================
	//==    RUN SEGMENTATION
	//==========================================
	//## Run the segmentation
	INFO_LOG("[PROC "<<m_procId<<"] - Running the hierarchical clustering segmenter...");	
	SLICData slicData_segm;
	SLICSegmenter::FindSegmentation(
		*slicData_init, slicData_segm,
		m_spMergingRegPar, m_spMergingUse2ndNeighbours,
		m_spMergingNSegmentsToStop,m_spMergingRatio,
		m_spMergingMaxDissRatio,m_spMergingMaxDissRatio2ndNeighbours,m_spMergingDissThreshold,
		m_spMergingIncludeSpatialPars, m_spMergingUseRobustPars, m_spMergingAddCurvDist
	);

	//## Get segmentation results	
	INFO_LOG("[PROC "<<m_procId<<"] - Computing the segmented map from slic segmented data...");
	bool normalizeSegmImg= true;
	bool binarizeSegmImg= true;
	Image* segmentedImg= SLIC::GetSegmentedImage(inputImg,slicData_segm.regions,Region::eSignalTag,normalizeSegmImg,binarizeSegmImg);
	if(!segmentedImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute the segmented image from slic segmented data!");
		delete slicData_init;
		slicData_init= 0;
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		if(edgeImg && !storeData){
			delete edgeImg;		
			edgeImg= 0;
		}
		return nullptr;
	}

	//## Clear-up
	if(slicData_init){
		delete slicData_init;
		slicData_init= 0;
	}
	if(laplImg && !storeData){
		delete laplImg;		
		laplImg= 0;
	}
	if(edgeImg && !storeData){
		delete edgeImg;		
		edgeImg= 0;
	}


	//## Finding blobs in masked image
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(
		sources, segmentedImg,
		bkgData, fgValue, fgValue, 
		m_NMinPix, findNegativeExcess, mergeBelowSeed, findNestedSources
	);
	if(status<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Finding sources in hierarchical algorithm segmented mask failed!");
		return nullptr;
	}

	
	//## Remove sources of negative excess (THIS METHOD SHOULD BE IMPROVED)
	if(inputImg->HasStats()){
		ImgStats* stats= inputImg->GetPixelStats();
		double imgMedian= stats->median;

		std::vector<size_t> sourcesToBeRemoved;				
		for(size_t k=0;k<sources.size();k++){	
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);

	}//close if
	else {
		WARN_LOG("[PROC "<<m_procId<<"] - Input image has no stats computed (hint: you must have computed them before!), cannot remove negative excess from sources!");
	}	

	//## Tag sources as extended
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetName(Form("Sext%d",(signed)k));
		sources[k]->SetType(Source::eExtended);
	}
	
	//## Add sources to extended sources
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());		
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	

	return segmentedImg;

}//close FindExtendedSources_HClust()




Image* SFinder::ComputeLaplacianImage(Image* inputImg){

	//Check input image
	if(!inputImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to given input image!");
		return nullptr;
	}

	//Compute laplacian image
	INFO_LOG("[PROC "<<m_procId<<"] - Computing Laplacian image ...");
	Image* laplImg= inputImg->GetLaplacianImage(true);
	if(!laplImg){
		ERROR_LOG("Failed to compute Laplacian image!");
		return nullptr;
	}

	//Compute laplacian image stats
	INFO_LOG("[PROC "<<m_procId<<"] - Compute Laplacian image stats...");
	if(laplImg->ComputeStats(true,false,false)<0){	
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute Laplacian image stats, returning nullptr!");
		delete laplImg;
		laplImg= 0;
		return nullptr;
	}
	INFO_LOG("[PROC "<<m_procId<<"] - Laplacian image stats");
	laplImg->PrintStats();

	return laplImg;

}//close ComputeLaplacianImage()

Image* SFinder::ComputeEdgeImage(Image* inputImg,int edgeModel){

	//Check input image
	if(!inputImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to given input image!");
		return nullptr;
	}

	//Compute edge image according to desired model
	Image* edgeImg= 0;
	if(edgeModel == eKirschEdge){
		INFO_LOG("[PROC "<<m_procId<<"] - Computing edge image using a Kirsch model...");
		edgeImg= inputImg->GetKirschImage();	
	}
	else if(edgeModel == eChanVeseEdge){
		INFO_LOG("[PROC "<<m_procId<<"] - Computing edge image using a Chan-Vese contour model...");	
		bool returnContourImg= true;
		edgeImg= ChanVeseSegmenter::FindSegmentation (
			inputImg, 0, returnContourImg,
			m_cvTimeStepPar, m_cvWindowSizePar, m_cvLambda1Par, m_cvLambda2Par, m_cvMuPar, m_cvNuPar, m_cvPPar, m_cvNIters
		);
	}
	else {
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid edge model selected!");
		return nullptr;
	}

	//Check if edge image computing failed
	if(!edgeImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute edge image!");
		return nullptr;
	}

	//Compute edge image stats
	INFO_LOG("[PROC "<<m_procId<<"] - Compute edge image stats...");
	if(edgeImg->ComputeStats(true,false,false)<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute edge image stats, returning nullptr!");
		delete edgeImg;
		edgeImg= 0;
		return nullptr;
	}

	INFO_LOG("[PROC "<<m_procId<<"] - Edgeness image stats");
	edgeImg->PrintStats();
	
	return edgeImg;

}//close ComputeEdgeImage()

Image* SFinder::FindExtendedSources_AC(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg,bool storeData){

	//## Check input image
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input image and/or to bkg/task data given!");
		return nullptr;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	INFO_LOG("[PROC "<<m_procId<<"] - Searching extended sources with the active contour method...");

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================	
	//## Compute saliency
	INFO_LOG("[PROC "<<m_procId<<"] - Computing image saliency map...");
	Image* saliencyImg= img->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,bkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(saliencyImg){
		if(storeData) m_SaliencyImg= saliencyImg;
	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute saliency map!");
		return nullptr;
	}
	
	//## Get saliency map optimal threshold
	INFO_LOG("[PROC "<<m_procId<<"] - Computing saliency map optimal threshold...");
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= saliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	if(TMath::IsNaN(signalThr) || fabs(signalThr)==TMath::Infinity()){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid numeric threshold returned as threshold computation failed!");
		return nullptr;
	}	

	//## Get saliency binarized image
	INFO_LOG("[PROC "<<m_procId<<"] - Thresholding the saliency map @ thr="<<signalThr<<" and compute binarized map...");
	double fgValue= 1;
	Image* signalMarkerImg= saliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	if(!signalMarkerImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to get saliency binarized map!");
		if(saliencyImg && !storeData){
			delete saliencyImg;
			saliencyImg= 0;
		}
		return nullptr;
	}
		
	//Delete saliency if not needed
	if(saliencyImg && !storeData){
		delete saliencyImg;
		saliencyImg= 0;
	}

	//## If binarized mage is empty (e.g. only background) do not run contour algorithm
	if(signalMarkerImg->GetMaximum()<=0){
		WARN_LOG("[PROC "<<m_procId<<"] - No signal objects detected in saliency map (only background), will not run active contour (NB: no extended sources detected in this image!)");
		delete signalMarkerImg;
		signalMarkerImg= 0;
		return nullptr;
	}

	//==========================================
	//==    RUN ACTIVE CONTOUR SEGMENTATION
	//==========================================

	//## Compute segmented image
	Image* segmentedImg= 0;
	bool returnContourImg= false;
	if(m_activeContourMethod==eChanVeseAC){//Standard ChanVese algorithm
		segmentedImg= ChanVeseSegmenter::FindSegmentation (
			img, signalMarkerImg, returnContourImg,
			m_cvTimeStepPar,m_cvWindowSizePar,m_cvLambda1Par,m_cvLambda2Par,m_cvMuPar,m_cvNuPar,m_cvPPar,m_cvNIters
		);
	}
	else if(m_activeContourMethod==eLRAC){//LRAC algorithm (with Chan-vese energy)
		segmentedImg= LRACSegmenter::FindSegmentation (
			img, signalMarkerImg,
			m_lracNIters,m_lracLambdaPar,m_lracRadiusPar,m_lracEpsPar
		);
	}
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid active contour method specified ("<<m_activeContourMethod<<")!");
		delete signalMarkerImg;
		signalMarkerImg= 0;
		return nullptr;
	}

	//Delete signal mask
	if(signalMarkerImg){
		delete signalMarkerImg;
		signalMarkerImg= 0;
	}

	if(!segmentedImg){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute Active Contour image segmentation!");
		return nullptr;
	}
	
	//## Finding blobs in masked image
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(
		sources, segmentedImg,
		bkgData, fgValue, fgValue, 
		m_NMinPix, findNegativeExcess, mergeBelowSeed, findNestedSources
	);
	if(status<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Finding sources in active contour segmented mask failed!");
		return nullptr;
	}

	
	
	//## Remove sources of negative excess (because Chan-Vese detects them) (THIS METHOD SHOULD BE IMPROVED)
	if(inputImg->HasStats()){
		ImgStats* stats= inputImg->GetPixelStats();
		double imgMedian= stats->median;

		std::vector<size_t> sourcesToBeRemoved;				
		for(size_t k=0;k<sources.size();k++){	
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);

	}//close if
	else {
		WARN_LOG("[PROC "<<m_procId<<"] - Input image has no stats computed (hint: you must have computed them before!), cannot remove negative excess from sources!");
	}	


	//## Tag sources as extended
	for(size_t k=0;k<sources.size();k++) {
		sources[k]->SetName(Form("Sext%d",(signed)k));
		sources[k]->SetType(Source::eExtended);
	}

	//## Add sources to extended sources
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());		
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	
	return segmentedImg;

}//close FindExtendedSources_AC()



Image* SFinder::FindExtendedSources_WT(Image* inputImg,TaskData* taskData,Image* searchedImg){

	//## Check input image
	if(!inputImg || !taskData){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input image and/or task data given!");
		return nullptr;
	}
	Image* img= inputImg;
	if(searchedImg) img= searchedImg;
	
	//## Find extended sources in the scales of the residual image (with POINT-LIKE SOURCES removed)
	INFO_LOG("[PROC "<<m_procId<<"] - Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ...");
	std::vector<Image*> wt_extended= img->GetWaveletDecomposition(m_wtScaleExtended);
	
	std::vector<Source*> sources;
	int status= FindSources(
		sources,
		inputImg,
		m_SeedThr,m_MergeThr,
		wt_extended[m_wtScaleExtended]
	);

	//## Clear-up
	for(size_t i=0;i<wt_extended.size();i++){
		if(wt_extended[i]) {
			delete wt_extended[i];
			wt_extended[i]= 0;
		}
	}

	if(status<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Extended source finding failed!");
		return nullptr;
	}
	
	
	//## Tag sources as extended
	int nSources= static_cast<int>( sources.size() );		
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" found...");

	for(size_t i=0;i<sources.size();i++){
		sources[i]->SetName(Form("Sext%d",(signed)(i)));	
		sources[i]->SetType(Source::eExtended);
	}

	//## Add sources to extended sources
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());		
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	

	//## Compute segmented map
	bool isBinary= false;
	bool invert= false;
	Image* segmMap= inputImg->GetSourceMask(sources,isBinary,invert);
	if(!segmMap){
		ERROR_LOG("[PROC "<<m_procId<<"] - Failed to compute segmented map!");
		return nullptr;
	}

	return segmMap;

}//close FindExtendedSources_WT()



int SFinder::SelectSources(std::vector<Source*>& sources){

	//## Apply source selection?
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) return 0;
	
	int nSelSources= 0;
	std::vector<Source*> sources_sel;

	for(int i=0;i<nSources;i++){	
		std::string sourceName= sources[i]->GetName();
		int sourceId= sources[i]->Id;
		long int NPix= sources[i]->NPix;
		double X0= sources[i]->X0;
		double Y0= sources[i]->Y0;

		//Is bad source (i.e. line-like blob, etc...)?
		if(!IsGoodSource(sources[i])) {
			INFO_LOG("[PROC "<<m_procId<<"] - Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!");
			sources[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource(sources[i]) ){
			INFO_LOG("[PROC "<<m_procId<<"] - Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ...");
			sources[i]->SetType(Source::ePointLike);
		}

		//Tag nested sources
		std::vector<Source*> nestedSources= sources[i]->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->GetName();
			int nestedSourceId= nestedSources[j]->Id;
			long int nestedNPix= nestedSources[j]->NPix;
			double nestedX0= nestedSources[j]->X0;
			double nestedY0= nestedSources[j]->Y0;

			if(!IsGoodSource(nestedSources[j])) {
				INFO_LOG("[PROC "<<m_procId<<"] - Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!");
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				INFO_LOG("[PROC "<<m_procId<<"] - Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ...");
				nestedSources[j]->SetType(Source::ePointLike);
			}
		}//end loop nested sources
			
		//Add source to the list	
		sources_sel.push_back(sources[i]);
		nSelSources++;
	}//end loop sources

	
	INFO_LOG("[PROC "<<m_procId<<"] - Added "<<nSelSources<<" bright sources to the selected list...");

	//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	sources.clear();
	sources.insert(sources.end(),sources_sel.begin(),sources_sel.end());
	sources_sel.clear();

	return 0;

}//close SelectSources()

bool SFinder::IsGoodSource(Source* aSource){
	
	if(!aSource) return false;

	//## Check for pixels 	
	if(aSource->NPix<=0 || (aSource->GetPixels()).size()<=0) return false;

	//## Check for line-like source
	if( (aSource->GetContours()).size()<=0) {
		WARN_LOG("[PROC "<<m_procId<<"] - No contour stored for this source, cannot perform check!");
		return true;
	}

	double BoundingBoxMin= ((aSource->GetContours())[0])->BoundingBoxMin;
	if(m_useMinBoundingBoxCut && BoundingBoxMin<m_SourceMinBoundingBox) {
		INFO_LOG("[PROC "<<m_procId<<"] - BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<m_SourceMinBoundingBox<<")");
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close IsGoodSource()

bool SFinder::IsPointLikeSource(Source* aSource){

	if(!aSource) return false;
	if(!aSource->HasParameters()) {
		WARN_LOG("[PROC "<<m_procId<<"] - No parameters are available for this source (did you compute them?)...point-like check cannot be performed!");
		return true;
	}

	std::string sourceName= aSource->GetName();
	int sourceId= aSource->Id;

	//Loop over contours and check if all of them have circular features
	bool isPointLike= true;
	std::vector<Contour*> contours= aSource->GetContours();

	for(unsigned int i=0;i<contours.size();i++){
		Contour* thisContour= contours[i];

		//Test circularity ratio: 1= circle
		if(m_useCircRatioCut && thisContour->CircularityRatio<m_psCircRatioThr) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<"<<m_psCircRatioThr<<")");
			isPointLike= false;
			break;
		}

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(m_useElongCut && thisContour->Elongation>m_psElongThr) {
			DEBUG_LOG("[PROC "<<m_procId<<"] - Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">"<<m_psElongThr<<")");
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if( m_useEllipseAreaRatioCut && thisContour->HasEllipseFit && 
				(thisContour->EllipseAreaRatio<m_psEllipseAreaRatioMinThr || thisContour->EllipseAreaRatio>m_psEllipseAreaRatioMaxThr) 
		) 
		{
			DEBUG_LOG("[PROC "<<m_procId<<"] - Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<" outside range ["<<m_psEllipseAreaRatioMinThr<<","<<m_psEllipseAreaRatioMaxThr<<"])");
			isPointLike= false;
			break;	
		}

	}//end contour loop
	
	//Check number of pixels
	DEBUG_LOG("[PROC "<<m_procId<<"] - Source (name="<<sourceName<<","<<"id="<<sourceId<<") (NPix="<<aSource->NPix<<">"<<m_psMaxNPix<<")");
	if(m_useMaxNPixCut && aSource->NPix>m_psMaxNPix){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nMaxPix cut (NPix="<<aSource->NPix<<">"<<m_psMaxNPix<<")");
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close IsPointLikeSource()

int SFinder::FitSources(std::vector<Source*>& sources){

	//Check given source list
	if(sources.empty()){
		WARN_LOG("[PROC "<<m_procId<<"] - Empty source list, nothing to be fitted!");
		return 0;
	}

	//## Loop over image sources and perform fitting stage for non-extended sources
	INFO_LOG("[PROC "<<m_procId<<"] - Loop over image sources and perform fitting stage for non-extended sources...");
		
	BlobPars blobPars;
	blobPars.bmaj= static_cast<int>( ceil(m_beamBmaj/m_pixSizeX) );//converted in pixels
	blobPars.bmin= static_cast<int>( ceil(m_beamBmin/m_pixSizeY) );//converted in pixels
	blobPars.bpa= m_beamBpa;

	for(size_t i=0;i<sources.size();i++){

		//Skip non point-like sources
		int sourceType= sources[i]->Type;
		if(sourceType!=Source::ePointLike) continue;
	
		//Fit mother source
		if(sources[i]->Fit(blobPars,m_fitMaxNComponents)<0) {
			WARN_LOG("[PROC "<<m_procId<<"] - Failed to fit source no. "<<i<<"!");
			continue;
		}

		//Fit nested sources
		//std::vector<Source*> nestedSources= sources[i]->GetNestedSources();	
		//for(size_t j=0;j<nestedSources.size();j++){
			//if(nestedSources[j] && nestedSources[j]->Deblend(m_deblendCurvThr,deblendComponentMinNPix)<0){//TO BE IMPLEMENTED!!!
			//	WARN_LOG("Failed to deblend nested source no. "<<j<<" of source no. "<<i<<"...");
			//}
		//}//end loop nested sources

	}//end loop sources
	
	return 0;

}//close FitSources()



Image* SFinder::ReadImage(FileInfo& info,std::string filename,std::string imgname,long int ix_min,long int ix_max,long int iy_min,long int iy_max)
{
	//## Check file
	bool match_extension= false;
	if(!SysUtils::CheckFile(filename,info,match_extension,"")){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid input file name specified (filename="<<filename<<"), invalid file path?!");
		return nullptr;
	}
	
	//## Read image from file
	bool readTile= (ix_min!=-1 && ix_max!=-1 && iy_min!=-1 && iy_max!=-1);
	Image* img= 0;
	
	//=== ROOT reading ===
	if(info.extension==".root"){// Read image from ROOT file
		TFile* inputFile = new TFile(filename.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("[PROC "<<m_procId<<"] - Cannot open input file "<<filename<<"!");
			return nullptr;
		}
		
		if(readTile){
			//Read full image
			Image* fullImg= (Image*)inputFile->Get(imgname.c_str());
			if(!fullImg){
				ERROR_LOG("[PROC "<<m_procId<<"] - Cannot get image "<<imgname<<" from input file "<<filename<<"!");
				return nullptr;
			}
			
			//Read tile
			INFO_LOG("[PROC "<<m_procId<<"] - Reading image tile (file="<<filename<<", hdu="<<m_fitsHDUId<<", range[xmin,xmax]=["<<ix_min<<","<<ix_max<<"], [ymin,ymax]=["<<iy_min<<","<<iy_max<<"])");
			img= fullImg->GetTile(ix_min,ix_max,iy_min,iy_max);	
			if(!img){
				ERROR_LOG("[PROC "<<m_procId<<"] - Failed to read image tile [xmin,xmax]=["<<ix_min<<","<<ix_max<<"], [ymin,ymax]=["<<iy_min<<","<<iy_max<<"]");
				delete fullImg;
				fullImg= 0;
				return nullptr;
			}				
		}//close if read tile
		else{
			img= (Image*)inputFile->Get(imgname.c_str());
			if(!img){
				ERROR_LOG("[PROC "<<m_procId<<"] - Cannot get image "<<imgname<<" from input file "<<filename<<"!");
				return nullptr;
			}	
		}

	}//close if

	//=== FITS reading ===
	else if(info.extension==".fits"){// Read image from FITS file
		img= new Image;

		int status= 0;
		if(readTile) {
			INFO_LOG("[PROC "<<m_procId<<"] - Reading image tile (file="<<filename<<", hdu="<<m_fitsHDUId<<", range[xmin,xmax]=["<<ix_min<<","<<ix_max<<"], [ymin,ymax]=["<<iy_min<<","<<iy_max<<"])");
			status= img->ReadFITS(filename,m_fitsHDUId,ix_min,ix_max,iy_min,iy_max); 
		}
		else {
			INFO_LOG("[PROC "<<m_procId<<"] - Reading image (file="<<filename<<", hdu="<<m_fitsHDUId<<")");
			status= img->ReadFITS(filename,m_fitsHDUId);
		}

		if(status<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to read image from input file "<<filename<<"!");
			if(img) {
				delete img;
				img= 0;
			}
			return nullptr;
		}
	}//close else if FITS reading

	//== Invalid extension ==
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid file extension detected (ext="<<info.extension<<")!");
		return nullptr;
	}
	
	return img;

}//close ReadImage()




ImgBkgData* SFinder::ComputeStatsAndBkg(Image* img){

	//## Check input img
	if(!img){
		ERROR_LOG("[PROC "<<m_procId<<"] - Null ptr to input image given!");
		return nullptr;
	}

	//## Compute stats
	INFO_LOG("[PROC "<<m_procId<<"] - Computing image stats...");
	bool computeRobustStats= true;
	bool skipNegativePix= false;
	bool forceRecomputing= false;
	if(img->ComputeStats(computeRobustStats,skipNegativePix,forceRecomputing)<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Stats computing failed!");
		return nullptr;
	}
	img->LogStats("INFO");

	//## Set local bkg grid/box
	//## If MetaData & beam info are available, interpret grid & box options as multiple of beam
	//## If no info is available (or use of beam info is off) interpret grid & box options as fractions wrt image size
	double boxSizeX= m_BoxSizeX;
	double boxSizeY= m_BoxSizeY;

	int pixelWidthInBeam= 0;
	if(m_UseBeamInfoInBkg){

		//Retrieve info from image metadata
		bool hasBeamData= false;
		if(img->HasMetaData()){
			pixelWidthInBeam= img->GetMetaData()->GetBeamWidthInPixel();	
			if(pixelWidthInBeam>0) {
				hasBeamData= true; 
				m_beamBmaj= img->GetMetaData()->Bmaj;
				m_beamBmin= img->GetMetaData()->Bmin;
				m_beamBpa= img->GetMetaData()->Bpa;
				m_pixSizeX= img->GetMetaData()->dX; 
				m_pixSizeY= img->GetMetaData()->dY; 
			}
		}
		
		//If beam data are not present in metadata, use those provided in the config file
		if(!hasBeamData){
			WARN_LOG("[PROC "<<m_procId<<"] - Using user-provided beam info to set bkg box size (beam info are not available/valid in image)...");
			pixelWidthInBeam= AstroUtils::GetBeamWidthInPixels(m_beamFWHM,m_beamFWHM,m_pixSize,m_pixSize);
			m_beamBmaj= m_beamFWHM;
			m_beamBmin= m_beamFWHM;
			m_beamBpa= m_beamTheta;
			m_pixSizeX= m_pixSize;
			m_pixSizeY= m_pixSize;	
		}

		if(pixelWidthInBeam<=0){
			ERROR_LOG("Invalid pixel width in beam computed from user-supplied beam info (beamFWHM,pixSize)=("<<m_beamFWHM<<","<<m_pixSize<<")!");
			return nullptr;
		}

		INFO_LOG("[PROC "<<m_procId<<"] - Setting bkg boxes as ("<<m_BoxSizeX<<","<<m_BoxSizeY<<") x beam (beam=#"<<pixelWidthInBeam<<" pixels) ...");
		boxSizeX= pixelWidthInBeam*m_BoxSizeX;
		boxSizeY= pixelWidthInBeam*m_BoxSizeY;

	}//close if use beam info
	else{
		WARN_LOG("[PROC "<<m_procId<<"] - Using image fractions to set bkg box size (beam info option is turned off)...");
		double Nx= static_cast<double>(img->GetNx());
		double Ny= static_cast<double>(img->GetNy());
		boxSizeX= m_BoxSizeX*Nx;
		boxSizeY= m_BoxSizeY*Ny;
		INFO_LOG("[PROC "<<m_procId<<"] - Setting bkg boxes to ("<<boxSizeX<<","<<boxSizeY<<") pixels ...");	
	}

	double gridSizeX= m_GridSizeX*boxSizeX;
	double gridSizeY= m_GridSizeY*boxSizeY;
	INFO_LOG("[PROC "<<m_procId<<"] - Setting grid size to ("<<gridSizeX<<","<<gridSizeY<<") pixels ...");

	//## Compute Bkg
	ImgBkgData* bkgData= img->ComputeBkg (
		m_BkgEstimator,
		m_UseLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY,
		m_Use2ndPassInLocalBkg,
		m_SkipOutliersInLocalBkg,m_SeedThr,m_MergeThr,m_NMinPix
	);

	if(!bkgData) {
		ERROR_LOG("[PROC "<<m_procId<<"] - Bkg computing failed!");
		return nullptr;
	}
		
	return bkgData;

}//close ComputeStatsAndBkg()


int SFinder::DrawSources(Image* image,std::vector<Source*>& sources){

	/*
	INFO_LOG("Drawing sources...");
	if(!image) return -1;

	bool useCurrentCanvas= false;
	bool drawFull= false;
	image->Plot(sources,useCurrentCanvas,drawFull,eRAINBOW,true);
	*/

	return 0;

}//close DrawSources()


int SFinder::SaveDS9RegionFile(){

	//## Open file
	FILE* fout= fopen(m_DS9CatalogFileName.c_str(),"w");

	//## Saving DS9 file region
	DEBUG_LOG("[PROC "<<m_procId<<"] - Saving DS9 region header...");
	fprintf(fout,"global color=red font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"image\n");

	DEBUG_LOG("[PROC "<<m_procId<<"] - Saving "<<m_SourceCollection.size()<<" sources to file...");

	for(unsigned int k=0;k<m_SourceCollection.size();k++){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Dumping DS9 region info for source no. "<<k<<" ...");
		std::string regionInfo= "";
		if(m_DS9RegionFormat==ePolygonRegion) regionInfo= m_SourceCollection[k]->GetDS9Region(true);
		else if(m_DS9RegionFormat==eEllipseRegion) regionInfo= m_SourceCollection[k]->GetDS9EllipseRegion(true);
		else {
			WARN_LOG("[PROC "<<m_procId<<"] - Invalid DS9RegionType given ("<<m_DS9RegionFormat<<")");
			return -1;
		}

		fprintf(fout,"%s\n",regionInfo.c_str());
	  	
	}//end loop sources
		
	DEBUG_LOG("[PROC "<<m_procId<<"] - Closing DS9 file region...");
	fclose(fout);

	return 0;

}//close SaveDS9RegionFile()


int SFinder::Save(){

	INFO_LOG("[PROC "<<m_procId<<"] - Storing results to file & catalog...");

	//Save DS9 regions?
	if(m_saveDS9Region && SaveDS9RegionFile()<0){
		WARN_LOG("[PROC "<<m_procId<<"] - Failed to save sources to DS9 region file!");
	}

	//Check ROOT output file
	if(!m_OutputFile) {
		WARN_LOG("[PROC "<<m_procId<<"] - Null ptr to output file, nothing will be saved in ROOT file!");
		return -1;
	}
	m_OutputFile->cd();

	//Save source tree?
	if(m_saveSources){
		DEBUG_LOG("[PROC "<<m_procId<<"] - Filling source ROOT TTree...");
		for(size_t k=0;k<m_SourceCollection.size();k++){
			m_Source= m_SourceCollection[k];
			m_SourceTree->Fill();
		}
		DEBUG_LOG("[PROC "<<m_procId<<"] - Writing tree to file...");
		m_SourceTree->Write();
	}
	
	//Save config?
	if(m_saveConfig){
		TTree* ConfigTree= ConfigParser::Instance().GetConfigTree("ConfigInfo");
		if(ConfigTree) ConfigTree->Write();
	}	

	//Save performance stats	
	if(m_PerfTree){
		m_PerfTree->Fill();
		m_PerfTree->Write();
	}

	//Save input image to file?
	if(m_saveInputMap && m_InputImg){
		INFO_LOG("Saving input map to file...");
		m_InputImg->SetNameTitle("img","img");
		m_InputImg->Write();
	}
	
	//Save residual map?
	if(m_saveResidualMap && m_ResidualImg){
		INFO_LOG("Saving residual map to file...");
		m_ResidualImg->SetNameTitle("img_residual","img_residual");
		m_ResidualImg->Write();
	}
	
	//Save significance map?
	if(m_saveSignificanceMap && m_SignificanceMap){
		INFO_LOG("Saving significance map to file...");
		m_SignificanceMap->SetNameTitle("img_significance","img_significance");
		m_SignificanceMap->Write();
	}

	//Save bkg & noise maps
	if(m_saveBkgMap && m_BkgData && m_BkgData->BkgMap){
		INFO_LOG("Saving bkg map to file...");
		(m_BkgData->BkgMap)->SetNameTitle("img_bkg","img_bkg");
		(m_BkgData->BkgMap)->Write();
	}
	if(m_saveNoiseMap && m_BkgData && m_BkgData->NoiseMap){
		(m_BkgData->NoiseMap)->SetNameTitle("img_rms","img_rms");
		(m_BkgData->NoiseMap)->Write();
	}

	//Save saliency map
	if(m_saveSaliencyMap && m_SaliencyImg){
		m_SaliencyImg->SetNameTitle("img_saliency","img_saliency");
		m_SaliencyImg->Write();
	}

	//Save Laplacian
	if(m_saveCurvatureMap && m_LaplImg){
		m_LaplImg->SetNameTitle("img_lapl","img_lapl");
		m_LaplImg->Write();
	}

	//Save Edgeness
	if(m_saveEdgenessMap && m_EdgeImg){
		m_EdgeImg->SetNameTitle("img_edge","img_edge");
		m_EdgeImg->Write();
	}

	//Save segmented map
	if(m_saveSegmentedMap && m_SegmImg){
		m_SegmImg->SetNameTitle("img_segm","img_segm");
		m_SegmImg->Write();
	}

	//## Close ROOT output file
	INFO_LOG("[PROC "<<m_procId<<"] - Closing output file...");
	if(m_OutputFile && m_OutputFile->IsOpen()) m_OutputFile->Close();

	INFO_LOG("[PROC "<<m_procId<<"] - End save to file");

	return 0;

}//close Save()


void SFinder::PrintPerformanceStats()
{
	INFO_LOG("===========================");
	INFO_LOG("===   PERFORMANCE INFO  ===");
	INFO_LOG("===========================");
	INFO_LOG("tot (ms)= "<<totTime);
	INFO_LOG("init (ms)= "<<initTime<<" ["<<initTime/totTime*100.<<"%]");
	INFO_LOG("read image (ms)= "<<readImageTime<<" ["<<readImageTime/totTime*100.<<"%]");
	INFO_LOG("source finding (ms)= "<<compactSourceTime<<" ["<<compactSourceTime/totTime*100.<<"%]");
	INFO_LOG("source selection (ms)= "<<sourceSelectionTime<<" ["<<sourceSelectionTime/totTime*100.<<"%]");
	INFO_LOG("img residual (ms)= "<<imgResidualTime<<" ["<<imgResidualTime/totTime*100.<<"%]");
	INFO_LOG("ext source finding (ms)= "<<extendedSourceTime<<" ["<<extendedSourceTime/totTime*100.<<"%]");
	INFO_LOG("source fitting (ms)= "<<sourceFitTime<<" ["<<sourceFitTime/totTime*100.<<"%]");
	INFO_LOG("save (ms)= "<<saveTime<<" ["<<saveTime/totTime*100.<<"%]");
	INFO_LOG("===========================");

}//close PrintPerformanceStats()


int SFinder::PrepareWorkerTasks()
{
	//## Generate a uuid for this job
	std::string jobId= CodeUtils::GenerateUUID();
	DEBUG_LOG("[PROC "<<m_procId<<"] - Generated jobId: "<<jobId);
	
	//## Get input image size
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
		return -1;
	}
	
	long int Nx= -1;
	long int Ny= -1;
	m_ImgXmin= 0;
	m_ImgYmin= 0;

	if(info.extension==".root"){// ROOT format

		//Read image from file
		TFile* inputFile= new TFile(m_InputFileName.c_str(),"READ");	
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
			return -1;
		}

		Image* inputImg= (Image*)inputFile->Get(m_InputImgName.c_str());
		if(!inputImg) {
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
			return -1;	
		}

		if(m_ReadTile){//READ TILE
			Nx= m_TileMaxX-m_TileMinX+1;
			Ny= m_TileMaxY-m_TileMinY+1;
			m_ImgXmin= inputImg->GetX(m_TileMinX);
			m_ImgYmin= inputImg->GetY(m_TileMinY);
		}
		else{//READ FULL MAP
			Nx= inputImg->GetNx();
			Ny= inputImg->GetNy();
			m_ImgXmin= inputImg->GetXmin();
			m_ImgYmin= inputImg->GetYmin();
		}

		if(inputImg){
			delete inputImg;
			inputImg= 0;
		}
	}//close if

	else if(info.extension==".fits"){//FITS
		
		if(m_ReadTile){//READ TILE
			Nx= m_TileMaxX-m_TileMinX+1;
			Ny= m_TileMaxY-m_TileMinY+1;
			m_ImgXmin= m_TileMinX;
			m_ImgYmin= m_TileMinY;	
		}
		else{//READ FULL MAP
			if(SysUtils::GetFITSImageSize(m_InputFileName,Nx,Ny)<0){
				ERROR_LOG("[PROC "<<m_procId<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
		}
	}//close else if		
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid/unsupported file extension ("<<info.extension<<") detected!");
		return -1;
	}

	INFO_LOG("[PROC "<<m_procId<<"] - Image size: "<<Nx<<"x"<<Ny<<", Image coord origin("<<m_ImgXmin<<","<<m_ImgYmin<<")");
	
	//==========================================
	//==    IMAGE PARTITION IN TILES
	//=========================================
	std::vector<long int> ix_min;
	std::vector<long int> ix_max;
	std::vector<long int> iy_min;
	std::vector<long int> iy_max;

	long int tileSizeX= Nx;
	long int tileSizeY= Ny;
	double tileStepSizeX= 1;
	double tileStepSizeY= 1;
	long int tileOverlapX= 0;
	long int tileOverlapY= 0;
	if(m_splitInTiles){
		tileSizeX= m_TileSizeX;
		tileSizeY= m_TileSizeY;
		tileStepSizeX= m_TileStepSizeX;
		tileStepSizeY= m_TileStepSizeY;
		if(m_UseTileOverlap){
			tileOverlapX= m_TileStepSizeX;
			tileOverlapY= m_TileStepSizeY;
		}
	}

	INFO_LOG("[PROC "<<m_procId<<"] - Computing tile partition: tileSize("<<tileSizeX<<","<<tileSizeY<<"), tileOverlap("<<tileOverlapX<<","<<tileOverlapY<<")");
	if(MathUtils::Compute2DGrid(ix_min,ix_max,iy_min,iy_max,Nx,Ny,tileSizeX,tileSizeY,tileStepSizeX,tileStepSizeY)<0){
		WARN_LOG("[PROC "<<m_procId<<"] - Failed to compute a 2D partition from input image!");
		return -1;
	}
	int nExpectedTasks= ix_min.size()*iy_min.size();
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nExpectedTasks<<" expected number of tasks ("<<ix_min.size()<<"x"<<iy_min.size()<<")");

	//## Compute worker tasks (check max number of tasks per worker)
	INFO_LOG("[PROC "<<m_procId<<"] - Computing worker task list...");
	TaskData* aTaskData= 0;
	long int workerCounter= 0;

	for(size_t j=0;j<iy_min.size();j++){
		for(size_t i=0;i<ix_min.size();i++){
			//Assign worker
			INFO_LOG("[PROC "<<m_procId<<"] - Assign task ("<<i<<","<<j<<") to worker no. "<<workerCounter<<"...");
				
			aTaskData= new TaskData;
			///aTaskData->filename= m_InputFileName;
			//aTaskData->jobId= jobId;
			aTaskData->workerId= workerCounter;
			//aTaskData->taskId= workerCounter;
			aTaskData->SetTile(
				ix_min[i] + m_TileMinX, ix_max[i] + m_TileMinX,
				iy_min[j] + m_TileMinY, iy_max[j] + m_TileMinY
			);
			//aTaskData->ix_min= ix_min[i] + m_TileMinX;
			//aTaskData->ix_max= ix_max[i] + m_TileMinX;
			//aTaskData->iy_min= iy_min[j] + m_TileMinY;
			//aTaskData->iy_max= iy_max[j] + m_TileMinY;
			m_taskDataPerWorkers[workerCounter].push_back(aTaskData);

			if(workerCounter>=m_nProc-1) workerCounter= 0;
			else workerCounter++;
		}//end loop x
	}//end loop y

	
	//Fill neighbor task list
	std::vector<int> workerIds;

	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		//Add only processors with tasks
		workerIds.push_back(i);
		
		//Loop over tasks present in this worker
		int nTasksInWorker= static_cast<int>(m_taskDataPerWorkers[i].size()); 
		
		for(int j=0;j<nTasksInWorker;j++){
			TaskData* task= m_taskDataPerWorkers[i][j];

			/*
			long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
			long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
			long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
			long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			*/

			//Find first neighbors among tasks inside the same worker
			for(int k=j+1;k<nTasksInWorker;k++){
				if(j==k) continue;
				TaskData* task_N= m_taskDataPerWorkers[i][k];
				bool areNeighbors= task->IsTaskTileNeighbor(task_N);		
				if(areNeighbors){
					task->AddNeighborInfo(k,i);
					task_N->AddNeighborInfo(j,i);
				}
	
				/*
				long int next_ix_min= m_taskDataPerWorkers[i][k]->ix_min;
				long int next_ix_max= m_taskDataPerWorkers[i][k]->ix_max;
				long int next_iy_min= m_taskDataPerWorkers[i][k]->iy_min;
				long int next_iy_max= m_taskDataPerWorkers[i][k]->iy_max;
				
				bool isAdjacentInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1 || (ix_min==next_ix_min && ix_max==next_ix_max));
				bool isAdjacentInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1 || (iy_min==next_iy_min && iy_max==next_iy_max));
				bool isOverlappingInX= ( (next_ix_min>=ix_min && next_ix_min<=ix_max) || (next_ix_max>=ix_min && next_ix_max<=ix_max) );
				bool isOverlappingInY= ( (next_iy_min>=iy_min && next_iy_min<=iy_max) || (next_iy_max>=iy_min && next_iy_max<=iy_max) );
				bool isAdjacent= isAdjacentInX && isAdjacentInY;
				bool isOverlapping= isOverlappingInX && isOverlappingInY;
				
				std::stringstream ss;	
				ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", Task "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"], NextTask "<<k<<"["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"] ==> isAdjacentInX? "<<isAdjacentInX<<", isAdjacentInY? "<<isAdjacentInY;
				if(m_procId==MASTER_ID) INFO_LOG(ss.str());
				
				if(isAdjacent || isOverlapping) {
					(m_taskDataPerWorkers[i][j]->neighborTaskId).push_back(k);
					(m_taskDataPerWorkers[i][k]->neighborTaskId).push_back(j);
					(m_taskDataPerWorkers[i][j]->neighborWorkerId).push_back(i);
					(m_taskDataPerWorkers[i][k]->neighborWorkerId).push_back(i);
				}
				*/
			}//end loop next task in worker


			//Find neighbors across workers
			for(size_t s=i+1;s<m_taskDataPerWorkers.size();s++){
				for(size_t t=0;t<m_taskDataPerWorkers[s].size();t++){
					TaskData* task_N= m_taskDataPerWorkers[s][t];
					bool areNeighbors= task->IsTaskTileNeighbor(task_N);		
					if(areNeighbors){
						task->AddNeighborInfo(t,s);
						task_N->AddNeighborInfo(j,i);
					}

					/*
					long int next_ix_min= m_taskDataPerWorkers[s][t]->ix_min;
					long int next_ix_max= m_taskDataPerWorkers[s][t]->ix_max;
					long int next_iy_min= m_taskDataPerWorkers[s][t]->iy_min;
					long int next_iy_max= m_taskDataPerWorkers[s][t]->iy_max;
				
					bool isAdjacentInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1 || (ix_min==next_ix_min && ix_max==next_ix_max));
					bool isAdjacentInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1 || (iy_min==next_iy_min && iy_max==next_iy_max));
					bool isOverlappingInX= ( (next_ix_min>=ix_min && next_ix_min<=ix_max) || (next_ix_max>=ix_min && next_ix_max<=ix_max) );
					bool isOverlappingInY= ( (next_iy_min>=iy_min && next_iy_min<=iy_max) || (next_iy_max>=iy_min && next_iy_max<=iy_max) );
					bool isAdjacent= isAdjacentInX && isAdjacentInY;
					bool isOverlapping= isOverlappingInX && isOverlappingInY;

					std::stringstream ss;	
					ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", Task "<<j<<", NextWorker no. "<<s<<", NextTask "<<t<<"["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"] ==> isAdjacentInX? "<<isAdjacentInX<<", isAdjacentInY? "<<isAdjacentInY;
					if(m_procId==MASTER_ID) INFO_LOG(ss.str());

					if(isAdjacent || isOverlapping) {
						(m_taskDataPerWorkers[i][j]->neighborTaskId).push_back(t);
						(m_taskDataPerWorkers[s][t]->neighborTaskId).push_back(j);
						(m_taskDataPerWorkers[i][j]->neighborWorkerId).push_back(s);
						(m_taskDataPerWorkers[s][t]->neighborWorkerId).push_back(i);
					}
					*/

				}//end loop tasks in next worker
			}//end loop workers 

		}//end loop tasks
	}//end loop workers

	int nWorkers= static_cast<int>(workerIds.size());
	std::stringstream ss;
	ss<<"[PROC "<<m_procId<<"] - # "<<nWorkers<<" workers {";
	for(int i=0;i<nWorkers;i++){
		ss<<workerIds[i]<<",";
	}
	ss<<"}";
	INFO_LOG(ss.str());
	

	
	//## Create a worker group (if MPI run is performed)
	m_workerRanks= -1;
	m_nWorkers= 1;

	#ifdef MPI_ENABLED
	if(m_mpiEnabled){

		//Get main processor group
		MPI_Comm_group(MPI_COMM_WORLD, &m_WorldGroup);
	
		// Construct a group containing all of the workers (proc with tasks assigned)
		MPI_Group_incl(m_WorldGroup, nWorkers, workerIds.data() , &m_WorkerGroup);

		// Create a new communicator based on the group
		int commTag= 10;
		MPI_Comm_create_group(MPI_COMM_WORLD, m_WorkerGroup, commTag, &m_WorkerComm);

		m_workerRanks = -1;
		m_nWorkers = -1;
	
		// If this rank isn't in the new communicator, it will be
		// MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
		// MPI_Comm_size is erroneous
		if (m_WorkerComm!=MPI_COMM_NULL) {
    	MPI_Comm_rank(m_WorkerComm, &m_workerRanks);
    	MPI_Comm_size(m_WorkerComm, &m_nWorkers);
		}
		else {
			WARN_LOG("[PROC "<<m_procId<<"] - Worker MPI communicator is null (this processor has no tasks and was not inserted in the worker group)!");
		}

	}//close if	
	#endif
	
	INFO_LOG("[PROC "<<m_procId<<"] - WORLD RANK/SIZE: "<<m_procId<<"/"<<m_nProc<<" WORKER RANK/SIZE: "<<m_workerRanks<<"/"<<m_nWorkers);


	//Print
	if(m_procId==MASTER_ID){
		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
				std::stringstream ss;	
				ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", ";

				long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
				long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
				long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
				long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			
				ss<<"Task no. "<<j<<" ["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"] neighbors{";
			
				for(size_t k=0;k<m_taskDataPerWorkers[i][j]->neighborTaskId.size();k++){
					long int neighborWorkerId= m_taskDataPerWorkers[i][j]->neighborWorkerId[k];
					long int neighborTaskId= m_taskDataPerWorkers[i][j]->neighborTaskId[k];
					long int next_ix_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_min;
					long int next_ix_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_max;
					long int next_iy_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_min;
					long int next_iy_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_max;

					ss<<"("<<neighborWorkerId<<","<<neighborTaskId<<") ["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"], ";
				}	
				ss<<"}";
				INFO_LOG(ss.str());	
			}//end loop tasks
		}//end loop workers
	}//close if MASTER


	bool hasTooManyTasks= false;
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		long int nTasksPerWorker= static_cast<long int>(m_taskDataPerWorkers[i].size());
		if(nTasksPerWorker>MAX_NTASKS_PER_WORKER){
			hasTooManyTasks= true;
			break;
		}
	}

	if(hasTooManyTasks){
		WARN_LOG("[PROC "<<m_procId<<"] - Too many tasks per worker (thr="<<MAX_NTASKS_PER_WORKER<<")");
		return -1;
	}

	return 0;

}//close PrepareWorkerTasks()


#ifdef MPI_ENABLED
int SFinder::GatherTaskDataFromWorkers()
{
	//## Put a barrier and collect all sources from workers in the master processor
	MPI_Barrier(MPI_COMM_WORLD);

	//## Sum and average all the elapsed timers across workers 
	INFO_LOG("[PROC "<<m_procId<<"] - Summing up and averaging he elapsed cpu timers across workers...");
	double initTime_sum;
	double readImageTime_sum;
	double compactSourceTime_sum;
	double sourceSelectionTime_sum;
	double imgResidualTime_sum;
	double extendedSourceTime_sum;
	/*
	MPI_Reduce(&initTime, &initTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&readImageTime, &readImageTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&compactSourceTime, &compactSourceTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&imgResidualTime, &imgResidualTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&sourceSelectionTime, &sourceSelectionTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&extendedSourceTime, &extendedSourceTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	*/
	MPI_Reduce(&initTime, &initTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, MPI_COMM_WORLD);
	MPI_Reduce(&readImageTime, &readImageTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, MPI_COMM_WORLD);
	MPI_Reduce(&compactSourceTime, &compactSourceTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, MPI_COMM_WORLD);
	MPI_Reduce(&imgResidualTime, &imgResidualTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, MPI_COMM_WORLD);
	MPI_Reduce(&sourceSelectionTime, &sourceSelectionTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, MPI_COMM_WORLD);
	MPI_Reduce(&extendedSourceTime, &extendedSourceTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, MPI_COMM_WORLD);
	
	INFO_LOG("[PROC "<<m_procId<<"] - CPU times (ms): {init="<<initTime<<", read="<<readImageTime<<", sourcefind="<<compactSourceTime<<", residual="<<imgResidualTime<<", sourcesel="<<sourceSelectionTime<<", extsourcefind="<<extendedSourceTime<<"}");

	if (m_procId == MASTER_ID) {
		initTime= initTime_sum;
		readImageTime= readImageTime_sum;
		compactSourceTime= compactSourceTime_sum;
		sourceSelectionTime= sourceSelectionTime_sum;
		imgResidualTime= imgResidualTime_sum;
		extendedSourceTime= extendedSourceTime_sum;
		INFO_LOG("[PROC "<<m_procId<<"] - Cumulative cpu times (ms): {init="<<initTime<<", read="<<readImageTime<<", sourcefind="<<compactSourceTime<<", residual="<<imgResidualTime<<", sourcesel="<<sourceSelectionTime<<", extsourcefind="<<extendedSourceTime);
	}

	//## Merge all sources found by workers in a unique collection
	INFO_LOG("[PROC "<<m_procId<<"] - Gathering task data found by workers in master processor...");
	int MSG_TAG= 1;
	if (m_procId == MASTER_ID) {//Receive data from the workers

		for (int i=1; i<m_nProc; i++) {
			//Check if this processor has tasks assigned, otherwise skip!
			if(m_taskDataPerWorkers[i].size()==0){
				INFO_LOG("[PROC "<<m_procId<<"] - No tasks assigned to process "<<i<<", nothing to be collected, skip to next worker...");
				continue;
			}
  
			//## Probe for an incoming message from process zero
			INFO_LOG("[PROC "<<m_procId<<"] - Probing for message from process "<<i<<"...");
    	MPI_Status status;

			if(MPI_Probe(i, MSG_TAG, MPI_COMM_WORLD, &status)==MPI_SUCCESS){
				INFO_LOG("[PROC "<<m_procId<<"] - a message has been probed from process "<<i<<" (tag="<< status.MPI_TAG << ", source " << status.MPI_SOURCE<<")");

    		//## When probe returns, the status object has the size and other
    		//## attributes of the incoming message. Get the message size
				INFO_LOG("[PROC "<<m_procId<<"] - Getting size of message received from process "<<i<<" ... ");
    	
				int rcvMsgSize= 0;
    		MPI_Get_count(&status, MPI_CHAR, &rcvMsgSize);

				//## Allocate a buffer to hold the incoming numbers
				DEBUG_LOG("[PROC "<<m_procId<<"] - Allocating a message of size "<<rcvMsgSize);
				if(rcvMsgSize<=0){
					ERROR_LOG("[PROC "<<m_procId<<"] - rcvMsg size is negative/null!");
					continue;
				}
				char* recvBuffer= (char*)malloc(rcvMsgSize);
    		
    		//## Now receive the message with the allocated buffer
    		//MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, i, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    		INFO_LOG("[PROC "<<m_procId<<"] - Received a message of size "<<rcvMsgSize<<") from process "<<i);

				//## Update task data with received worker data	
				bool isTaskCollectionPreAllocated= true;
				if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[i],recvBuffer,rcvMsgSize,isTaskCollectionPreAllocated)<0 ){
					ERROR_LOG("[PROC "<<m_procId<<"] - Failed to decode recv message into task data list!");
    			if(recvBuffer) free(recvBuffer);
				}

				//## Free received buffer
    		if(recvBuffer) free(recvBuffer);
			}//close if
			else{
				ERROR_LOG("[PROC "<<m_procId<<"] - Message probing from process "<<i<<" failed!");
				return -1;
				//continue;
			}
		}//end loop workers
	}//close if

	else {//Send data to master
		//## First encode taskData in protobuf
		INFO_LOG("[PROC "<<m_procId<<"] - Encoding task data collection to buffer...");
		long int msg_size= 0;
		char* msg= Serializer::TaskDataCollectionToCharArray(msg_size,m_taskDataPerWorkers[m_procId]);
		if(!msg){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to encode task data to protobuf!");
			return -1;
		}

		//## Send buffer to master processor	
		INFO_LOG("[PROC "<<m_procId<<"] - Sending task data to master process (msg: "<<msg<<", size="<<msg_size<<")...");
		MPI_Send((void*)(msg),msg_size, MPI_CHAR, MASTER_ID, MSG_TAG, MPI_COMM_WORLD);

		//## Free buffer
		if(msg) free(msg);
	}//close else

	MPI_Barrier(MPI_COMM_WORLD);

	//## Print sources
	/*
	if (m_procId == 0) {
		INFO_LOG("[PROC "<<m_procId<<"] - Printing aggregated sources...");
		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
			
			
			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
				INFO_LOG("[PROC "<<m_procId<<"] - Sources NOT at edge found in process "<<i<<", task "<<j<<" ...");
				for(size_t k=0;k<(m_taskDataPerWorkers[i][j]->sources).size();k++){
					(m_taskDataPerWorkers[i][j]->sources)[k]->Print();
				}//end loop sources
	
				INFO_LOG("[PROC "<<m_procId<<"] - Sources at edge found in process "<<i<<", task "<<j<<" ...");
				for(size_t k=0;k<(m_taskDataPerWorkers[i][j]->sources_edge).size();k++){
					(m_taskDataPerWorkers[i][j]->sources_edge)[k]->Print();
				}//end loop sources

			}//end loop tasks in worker
		}//end process loop
	}//close if
	*/

	return 0;

}//close GatherTaskDataFromWorkers()
#endif


int SFinder::MergeTaskData()
{
	//## Update sources in list
	//## NB: Push to collection only sources NON at edge
	if (m_procId == 0) {

		//Print task data
		INFO_LOG("[PROC "<<m_procId<<"] - Printing task data...");
		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
			
			std::stringstream ss;
			ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", ";
			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
				long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
				long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
				long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
				long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
				m_xmin= ix_min;
				m_xmax= ix_max;
				m_ymin= iy_min;
				m_ymax= iy_max;

				//double x_min= m_taskDataPerWorkers[i][j]->x_min;
				//double x_max= m_taskDataPerWorkers[i][j]->x_max;
				//double y_min= m_taskDataPerWorkers[i][j]->y_min;
				//double y_max= m_taskDataPerWorkers[i][j]->y_max;
				//m_xmin= x_min;
				//m_xmax= x_max;
				//m_ymin= y_min;
				//m_ymax= y_max;

				
				if(m_TaskInfoTree) m_TaskInfoTree->Fill();

				//ss<<"Task no. "<<j<<", PixelRange["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"], PhysCoordRange["<<x_min<<","<<x_max<<"] ["<<y_min<<","<<y_max<<"], ";
				ss<<"Task no. "<<j<<", PixelRange["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"], ";
			}//end loop tasks
			INFO_LOG(ss.str());			
	
		}//end loop workers

		//Add sources
		int nCompactSources= 0;
		int nExtendedSources= 0;
		int nPointLikeSources= 0;
		int nUnknown= 0;
		int nEdgeSources= 0;
		int nSources= 0;
		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){

				//Process sources not at edge
				for(size_t k=0;k<(m_taskDataPerWorkers[i][j]->sources).size();k++){
					int sourceType= (m_taskDataPerWorkers[i][j]->sources)[k]->Type;
					if(sourceType==Source::eCompact) nCompactSources++;
					else if(sourceType==Source::eExtended) nExtendedSources++;
					else if(sourceType==Source::ePointLike) nPointLikeSources++;
					else nUnknown++;
					nSources++;
					m_SourceCollection.push_back( (m_taskDataPerWorkers[i][j]->sources)[k] );
				}//end loop sources per task

				//Add edge sources not merged (if merging was not selected as option)
				for(size_t k=0;k<(m_taskDataPerWorkers[i][j]->sources_edge).size();k++){
					int sourceType= (m_taskDataPerWorkers[i][j]->sources_edge)[k]->Type;
					if(sourceType==Source::eCompact) nCompactSources++;
					else if(sourceType==Source::eExtended) nExtendedSources++;
					else if(sourceType==Source::ePointLike) nPointLikeSources++;
					else nUnknown++;
					nSources++;
					nEdgeSources++;
					if(!m_mergeSourcesAtEdge) m_SourceCollection.push_back( (m_taskDataPerWorkers[i][j]->sources_edge)[k] );
				}//end loop sources per task		

			}//end loop tasks
		}//end loop workers

		//Add merged sources at edge	
		m_SourceCollection.insert(m_SourceCollection.end(),m_SourcesMergedAtEdges.begin(),m_SourcesMergedAtEdges.end());
	
		INFO_LOG("[PROC "<<m_procId<<"] - #"<<nSources<<" sources found in total (#"<<nCompactSources<<" compact, #"<<nPointLikeSources<<" point-like, #"<<nExtendedSources<<" extended, #"<<nUnknown<<" unknown/unclassified, #"<<nEdgeSources<<" edge sources, #"<<m_SourcesMergedAtEdges.size()<<" merged at edges), #"<<m_SourceCollection.size()<<" sources added to collection ...");
		
	}//close if

	return 0;

}//close MergeTaskData()


int SFinder::MergeSourcesAtEdge()
{
	//## Executed only by master processor
	if(m_procId != MASTER_ID) return 0;

	//## Loop over edge sources and merge them
	//## TBD: Merge all sources regardless of their tag (compact or extended)
	struct MergedSourceInfo {
		long int source_index;
		long int worker_index;
		long int task_index;
		MergedSourceInfo(long int sindex,long int windex,long int tindex)
			: source_index(sindex), worker_index(windex), task_index(tindex)
		{}
		bool operator==(const MergedSourceInfo& obj) const {
			bool areEqual= ( 
				(obj.source_index==source_index) && 	
				(obj.worker_index==worker_index) &&
				(obj.task_index==task_index)
			);
    	return areEqual;
    }
	};//close MergedSourceInfo

	//## Fill list of edge sources to be merged and fill corresponding Graph
	INFO_LOG("[PROC "<<m_procId<<"] - Fill list of edge sources to be merged and fill corresponding graph data struct...");
	std::vector<MergedSourceInfo> sourcesToBeMerged;
	Graph mergedSourceGraph;
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
			int nEdgeSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources_edge).size());
			for(int k=0;k<nEdgeSources;k++){
				sourcesToBeMerged.push_back(MergedSourceInfo(k,i,j));
				mergedSourceGraph.AddVertex();
			}
		}
	}
	
	//## Return if there are no edge sources
	if(sourcesToBeMerged.empty()){
		INFO_LOG("[PROC "<<m_procId<<"] - No edge sources to be merged, nothing to be done...");
		return 0;
	}

	//## Find adjacent edge sources	
	INFO_LOG("[PROC "<<m_procId<<"] - Finding adjacent edge sources (#"<<sourcesToBeMerged.size()<<" edge sources present, graph nvertex="<<mergedSourceGraph.GetNVertexes()<<") ...");

	int itemPos= -1;
	for(size_t i=0;i<sourcesToBeMerged.size()-1;i++){	
		long int sindex= sourcesToBeMerged[i].source_index;
		long int windex= sourcesToBeMerged[i].worker_index;
		long int tindex= sourcesToBeMerged[i].task_index; 
		Source* source= (m_taskDataPerWorkers[windex][tindex]->sources_edge)[sindex];
		
		//Loop neighbors
		for(size_t j=i+1;j<sourcesToBeMerged.size();j++){	
			long int sindex_neighbor= sourcesToBeMerged[j].source_index;
			long int windex_neighbor= sourcesToBeMerged[j].worker_index;
			long int tindex_neighbor= sourcesToBeMerged[j].task_index; 
			Source* source_neighbor= (m_taskDataPerWorkers[windex_neighbor][tindex_neighbor]->sources_edge)[sindex_neighbor];
			
			//Check if main worker tile is physically neighbor to this
			//If not they cannot be adjacent, so skip the following check
			if(windex!=windex_neighbor && !CodeUtils::FindItem(m_taskDataPerWorkers[windex][tindex]->neighborWorkerId,windex_neighbor,itemPos)){
				INFO_LOG("Worker id "<<windex_neighbor<<" is not physically neighbor to worker "<<windex<<", so skip the source adjacency check and go to next...");
				continue;
			}

			/*
			//Check if bouding boxes are overlapping
			//NB: If not skip the adjacency check
			bool areBoundingBoxesOverlapping= source->CheckBoxOverlapping(source_neighbor);
			if(!areBoundingBoxesOverlapping){
				DEBUG_LOG("Sources (i,j)=("<<i<<" {"<<sindex<<","<<windex<<","<<tindex<<"} , "<<j<<" {"<<sindex_neighbor<<","<<windex_neighbor<<","<<tindex_neighbor<<"}) have NON-overlapping bounding boxes, skip the adjacency check...");
			}
			*/

			//Check is sources are adjacent
			//NB: This is time-consuming (N1xN2 more or less)!!!
			bool areAdjacentSources= source->IsAdjacentSource(source_neighbor);
			if(!areAdjacentSources) continue;

			//If they are adjacent add linking in graph
			INFO_LOG("Sources (i,j)=("<<i<<" {"<<sindex<<","<<windex<<","<<tindex<<"} , "<<j<<" {"<<sindex_neighbor<<","<<windex_neighbor<<","<<tindex_neighbor<<"}) are adjacent and selected for merging...");
			mergedSourceGraph.AddEdge(i,j);

		}//end loop sources

	}//end loop sources


	//## Find all connected components in graph corresponding to 
	//## edge sources to be merged
	INFO_LOG("Find all connected components in graph corresponding to edge sources to be merged...");
	std::vector<std::vector<int>> connected_source_indexes;
	mergedSourceGraph.GetConnectedComponents(connected_source_indexes);
	INFO_LOG("#"<<connected_source_indexes.size()<<"/"<<sourcesToBeMerged.size()<<" edge sources will be left after merging...");
		
	//## Now merge the sources
	std::vector<int> sourcesToBeRemoved;
	bool copyPixels= true;//do not create memory for new pixels
	bool checkIfAdjacent= false;//already done before
	bool computeStatPars= false;//do not compute stats& pars at each merging
	bool computeMorphPars= false;
	bool computeRobustStats= true;
	bool forceRecomputing= false;//no need to re-compute moments (already updated in AddPixel())

	m_SourcesMergedAtEdges.clear();

	INFO_LOG("Merging edge sources and adding them to collection...");
	for(size_t i=0;i<connected_source_indexes.size();i++){
		if(connected_source_indexes[i].empty()) continue;

		//Get source id=0 of this component
		int index= connected_source_indexes[i][0];
		long int sindex= sourcesToBeMerged[index].source_index;
		long int windex= sourcesToBeMerged[index].worker_index;
		long int tindex= sourcesToBeMerged[index].task_index; 
		Source* source= (m_taskDataPerWorkers[windex][tindex]->sources_edge)[sindex];
		sourcesToBeRemoved.push_back(index);

		//Create a new source which merges the two
		Source* merged_source= new Source;
		*merged_source= *source;

		//Merge other sources in the group if any 
		int nMerged= 0;
		int source_type= source->Type;
		for(size_t j=1;j<connected_source_indexes[i].size();j++){
			int index_adj= connected_source_indexes[i][j];
			long int sindex_adj= sourcesToBeMerged[index_adj].source_index;
			long int windex_adj= sourcesToBeMerged[index_adj].worker_index;
			long int tindex_adj= sourcesToBeMerged[index_adj].task_index; 
			Source* source_adj= (m_taskDataPerWorkers[windex_adj][tindex_adj]->sources_edge)[sindex_adj];
				
			int status= merged_source->MergeSource(source_adj,copyPixels,checkIfAdjacent,computeStatPars,computeMorphPars);
			if(status<0){
				WARN_LOG("Failed to merge sources (i,j)=("<<index<<" {"<<sindex<<","<<windex<<","<<tindex<<"} , "<<index_adj<<" {"<<sindex_adj<<","<<windex_adj<<","<<tindex_adj<<"}), skip to next...");
				continue;
			}
			nMerged++;

			//Add this source to the list of edge sources to be removed
			sourcesToBeRemoved.push_back(index_adj);

		}//end loop of sources to be merged in this component

		//If at least one was merged recompute stats & pars of merged source
		if(nMerged>0) {
			INFO_LOG("Recomputing stats & moments of merged source in merge group "<<i<<" after #"<<nMerged<<" merged source...");
			if(merged_source->ComputeStats(computeRobustStats,forceRecomputing)<0){
				WARN_LOG("Failed to compute stats for merged source in merge group "<<i<<"...");
				continue;
			}
			if(merged_source->ComputeMorphologyParams()<0){
				WARN_LOG("Failed to compute morph pars for merged source in merge group "<<i<<"...");
				continue;
			}
		}//close if

		//Add merged source to collection
		m_SourcesMergedAtEdges.push_back(merged_source);

	}//end loop number of components

	INFO_LOG("#"<<m_SourcesMergedAtEdges.size()<<" merged sources at edge...");
	
	return 0;

}//close MergeSourcesAtEdge()


int SFinder::FindSourcesAtEdge()
{	
	//## Find if sources (both compact and extended) are at tile edge
	//## Those found at the edge are removed from the list and added to the edge list for further processing
	float xmin_s, xmax_s, ymin_s, ymax_s;
	
	//Loop over workers
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		//Loop over tasks per worker
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
	
			//Loop over sources found
			std::vector<Source*> sources_not_at_edges;
			int nSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources).size());
			if(nSources<=0) continue;
			int nEdgeSources= 0;

			for(int k=0;k<nSources;k++){
				//Get source coordinate range
				(m_taskDataPerWorkers[i][j]->sources)[k]->GetSourceRange(xmin_s,xmax_s,ymin_s,ymax_s);

				//Check if source is at the edge of its tile
				long int workerId= m_taskDataPerWorkers[i][j]->workerId;
				float xmin_tile= (m_taskDataPerWorkers[i][j])->ix_min;// + m_ImgXmin;
				float xmax_tile= (m_taskDataPerWorkers[i][j])->ix_max;// + m_ImgXmin;
				float ymin_tile= (m_taskDataPerWorkers[i][j])->iy_min;// + m_ImgYmin;
				float ymax_tile= (m_taskDataPerWorkers[i][j])->iy_max;// + m_ImgYmin; 
				bool isAtTileEdge= (m_taskDataPerWorkers[i][j]->sources)[k]->IsAtBoxEdge(xmin_tile,xmax_tile,ymin_tile,ymax_tile);
				/*
				long int xmin_tile= (m_taskDataPerWorkers[i][j])->ix_min;
				long int xmax_tile= (m_taskDataPerWorkers[i][j])->ix_max;
				long int ymin_tile= (m_taskDataPerWorkers[i][j])->iy_min;
				long int ymax_tile= (m_taskDataPerWorkers[i][j])->iy_max; 
				bool isAtTileEdgeX= (xmin_s==xmin_tile || xmax_s==xmax_tile);
				bool isAtTileEdgeY= (ymin_s==ymin_tile || ymax_s==ymax_tile);
				bool isAtTileEdge= (isAtTileEdgeX || isAtTileEdgeY);
				INFO_LOG("[PROC "<<m_procId<<"] - workerId="<<workerId<<", check if compact source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is at edge of its tile (x["<<xmin_tile<<","<<xmax_tile<<"] y["<<ymin_tile<<","<<ymax_tile<<"]), isAtTileEdgeX="<<isAtTileEdgeX<<", isAtTileEdgeY="<<isAtTileEdgeY<<", isAtTileEdge="<<isAtTileEdge);
				*/


				if(isAtTileEdge){
					INFO_LOG("[PROC "<<m_procId<<"] - workerId="<<workerId<<", source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is at edge of its tile (x["<<xmin_tile<<","<<xmax_tile<<"] y["<<ymin_tile<<","<<ymax_tile<<"])");
				}


				//Check if source is inside neighbour tile, e.g. is in overlapping area
				//NB: This is done only if source is not found at tile border previously
				bool isInOverlapArea= false;


				if(!isAtTileEdge){	
					//DEBUG_LOG("[PROC "<<m_procId<<"] - workerId="<<workerId<<", check if compact source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is inside neighbour tile...");
				
					for(size_t l=0;l<(m_taskDataPerWorkers[i][j]->neighborWorkerId).size();l++){	
						long int neighborTaskId= (m_taskDataPerWorkers[i][j]->neighborTaskId)[l];
						long int neighborWorkerId= (m_taskDataPerWorkers[i][j]->neighborWorkerId)[l];
					
						float xmin= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_min;//+ m_ImgXmin;
						float xmax= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_max;//+ m_ImgXmin;
						float ymin= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_min;//+ m_ImgYmin;
						float ymax= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_max;//+ m_ImgYmin;
						bool isOverlapping= (m_taskDataPerWorkers[i][j]->sources)[k]->HasBoxOverlap(xmin,xmax,ymin,ymax);
						if(isOverlapping){
							INFO_LOG("[PROC "<<m_procId<<"] - workerId="<<workerId<<", source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) overlaps with neighbor tile (x["<<xmin<<","<<xmax<<"] y["<<ymin<<","<<ymax<<"])");
							isInOverlapArea= true;
							break;
						}

						/*
						long int xmin= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_min;
						long int xmax= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_max;
						long int ymin= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_min;
						long int ymax= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_max;
		
						bool isOverlappingX= ( (xmin_s<=xmax && xmin_s>=xmin) || (xmax_s<=xmax && xmax_s>=xmin) );
						bool isOverlappingY= ( (ymin_s<=ymax && ymin_s>=ymin) || (ymax_s<=ymax && ymax_s>=ymin) );
							
						INFO_LOG("[PROC "<<m_procId<<"] - neighborWorkerId="<<neighborWorkerId<<", neighborTaskId="<<neighborTaskId<<", check if inside neighbor tile no. "<<j<<"(x["<<xmin<<","<<xmax<<"] y["<<ymin<<","<<ymax<<"]), isOverlappingX="<<isOverlappingX<<", isOverlappingY="<<isOverlappingY);
						if( isOverlappingX && isOverlappingY ){
							isInOverlapArea= true;
							break;
						}
						*/

					}//end loop neighbors
				}//close if

				//Tag source at edge is is located at the border of its tile or if it is located inside an overlapping area with another neighbor tile
				bool isAtEdge= (isAtTileEdge || isInOverlapArea);

				//Set edge flag in source
				if(isAtEdge) {
					(m_taskDataPerWorkers[i][j]->sources)[k]->SetEdgeFlag(true);
					(m_taskDataPerWorkers[i][j]->sources_edge).push_back( (m_taskDataPerWorkers[i][j]->sources)[k] );
					nEdgeSources++;
				}
				else {
					(m_taskDataPerWorkers[i][j]->sources)[k]->SetEdgeFlag(false);
					sources_not_at_edges.push_back( (m_taskDataPerWorkers[i][j]->sources)[k] );
				}
			}//end loop sources	

			//int nEdgeSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources_edge).size());
			INFO_LOG("[PROC "<<m_procId<<"] - #"<<nEdgeSources<<"/"<<nSources<<" sources are found at tile edges...");
	
			//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
			(m_taskDataPerWorkers[i][j]->sources).clear();
			(m_taskDataPerWorkers[i][j]->sources).insert( 
				(m_taskDataPerWorkers[i][j]->sources).end(),
				sources_not_at_edges.begin(),
				sources_not_at_edges.end()
			);
			sources_not_at_edges.clear();
			
		}//end loop tasks per worker
	}//end loop workers

	return 0;

}//close FindSourcesAtEdge()

}//close namespace
