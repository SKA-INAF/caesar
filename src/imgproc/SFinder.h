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
* @file SFinder.h
* @class SFinder
* @brief Source finder class
*
* Class to perform source finding 
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SFINDER_h
#define _SFINDER_h 1

#include <SysUtils.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

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
#include <limits>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif



namespace Caesar {

class Image;
class Source;
class ImgBkgData;
class TaskData;


class SFinder : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SFinder();


	public:
		
		/**
		* \brief Run source finder 
		*/
		int Run();
		
		
	public:
		/**
		* \brief Read only a tile from image
		*/
		int SetTileRead(double xmin,double xmax,double ymin,double ymax){
			if(xmin>=xmax || ymin>=ymax){
				m_ReadTile= false;
				return -1;
			}
			m_ReadTile= true;
			m_TileMinX= xmin;
			m_TileMaxX= xmax;
			m_TileMinY= ymin;
			m_TileMaxY= ymax;
			return 0;
		}	

	private:

		/**
		* \brief Clear data
		*/
		void Clear();

		/**
		* \brief Initialize class options from ConfigParser
		*/
		void InitOptions();
		/**
		* \brief Initialize class variables
		*/
		int Init();
		/**
		* \brief Save data to output file
		*/
		int Save();	


		/**
		* \brief Run source finder task
		*/
		int RunTask(TaskData* taskData,bool storeData=false);
		

		/**
		* \brief Set options from ConfigParser singleton
		*/
		int Configure();

		
		/**
		* \brief Compute Stats and Bkg info
		*/
		ImgBkgData* ComputeStatsAndBkg(Image* inputImg,bool useRange=false,double minThr=-std::numeric_limits<double>::infinity(),double maxThr=std::numeric_limits<double>::infinity());

		/**
		* \brief Read image
		*/
		Image* ReadImage(Caesar::FileInfo& info,std::string filename,std::string imgname="",long int ix_min=-1,long int ix_max=-1,long int iy_min=-1,long int iy_max=-1);

		/**
		* \brief Save detected sources to DS9 region file
		*/
		int SaveDS9RegionFile();

		/**
		* \brief Save detected source to catalog ascii file
		*/
		int SaveCatalogFile();

		/**
		* \brief Get smoothed image
		*/
		Image* ComputeSmoothedImage(Image* inputImg,int model);

		/**
		* \brief Find sources from input image
		*/
		int FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr,Image* searchedImg=0);
		/**
		* \brief Find compact sources
		*/
		Image* FindCompactSources(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData);

		/**
		* \brief Find compact sources
		*/
		Image* FindCompactSourcesRobust(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,int niter=10);

		/**
		* \brief Compute residual map
		*/
		Image* FindResidualMap(Image* inputImg,ImgBkgData* bkgData,std::vector<Source*> const & sources);

		/**
		* \brief Find extended sources
		*/
		Image* FindExtendedSources(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,bool storeData=false);
		
		/**
		* \brief Find extended sources with hierarchical clustering method
		*/
		Image* FindExtendedSources_HClust(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg=0,bool storeData=false);

		/**
		* \brief Find extended sources with Active Contour method
		*/
		Image* FindExtendedSources_AC(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg=0,bool storeData=false);

		/**
		* \brief Find extended sources with Wavelet Transform method
		*/
		Image* FindExtendedSources_WT(Image* inputImg,TaskData* taskData,Image* searchedImg=0);

		/**
		* \brief Find extended sources with Saliency Map thresholding method
		*/
		Image* FindExtendedSources_SalThr(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg=0,bool storeData=false);

		/**
		* \brief Compute blob mask image
		*/
		Image* ComputeBlobMaskImage(Image* inputImg);

		/**
		* \brief Compute edge image
		*/
		Image* ComputeEdgeImage(Image* inputImg,int model);

		/**
		* \brief Compute laplacian image
		*/
		Image* ComputeLaplacianImage(Image* inputImg);

		/**
		* \brief Select sources according to quality cuts given in configuration
		*/
		int SelectSources(std::vector<Source*>& sources);
		/**
		* \brief Tag a source as good or bad
		*/
		bool IsGoodSource(Source* aSource);
		/**
		* \brief Tag a source as point-like or not
		*/
		bool IsPointLikeSource(Source* aSource);

		/**
		* \brief Fit sources
		*/
		int FitSources(std::vector<Source*>& sources);

		#ifdef MPI_ENABLED
			/**
			* \brief Fit sources in parallel using MPI
			*/
			//int FitSourcesMPI(std::vector<Source*>& sources);	
		#endif

		/**
		* \brief Is fittable source
		*/
		bool IsFittableSource(Source* aSource);

		/**
		* \brief Print performance stats
		*/
		void PrintPerformanceStats();

		/**
		* \brief Prepare worker task data (for MPI run)
		*/
		int PrepareWorkerTasks();

		#ifdef MPI_ENABLED
			/**
			* \brief Collect task data from workers (for MPI run)
			*/
			int GatherTaskDataFromWorkers();
		#endif

		/**
		* \brief Merge sources found at each task in collections
		*/
		int MergeTaskData();

		/**
		* \brief Find sources at image edges (for MPI run)
		*/
		int FindSourcesAtEdge();

		/**
		* \brief Fit all task sources not found at tile edge (for MPI run)
		*/
		int FitTaskSources();

		/**
		* \brief Merge sources at edge (for MPI run)
		*/
		int MergeSourcesAtEdge();

		/**
		* \brief Merge task sources
		*/
		int MergeTaskSources(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData);

	public:
		
		
		//Input data
		std::string m_InputFileName;
		std::string m_InputImgName;
		std::string m_InputFileExtension;
		int m_InputFileType;
		Image* m_InputImg;
		int m_fitsHDUId;

		//Output data
		TApplication* m_Application;
		bool m_IsInteractiveRun;
		std::string m_OutputFileName;
		TFile* m_OutputFile;
		std::string m_catalogOutFileName;
		std::string m_catalogComponentsOutFileName;
		bool m_saveToFile;
		bool m_saveToCatalogFile;
		bool m_saveConfig;
		bool m_saveDS9Region;
		bool m_convertDS9RegionsToWCS;
		int m_ds9WCSType;
		std::string m_DS9CatalogFileName;
		int m_DS9RegionFormat;
		std::string m_DS9FitCatalogFileName;			
		TTree* m_SourceTree;
		bool m_saveSources;
		bool m_saveResidualMap;
		bool m_saveInputMap;
		bool m_saveSignificanceMap;
		bool m_saveBkgMap;
		bool m_saveNoiseMap;
		bool m_saveSaliencyMap;
		bool m_saveEdgenessMap;
		bool m_saveCurvatureMap;
		bool m_saveSegmentedMap;

		//Performance stats data
		TTree* m_PerfTree;
		double totTime;
		double initTime;
		double initTime_sum;
		double initTime_min;
		double initTime_max;
		double imageStatsTime;	
		double imageStatsTime_sum;
		double imageStatsTime_min;
		double imageStatsTime_max;
		double imageBkgTime;
		double imageBkgTime_sum;
		double imageBkgTime_min;
		double imageBkgTime_max;
		double readImageTime;
		double readImageTime_sum;
		double readImageTime_min;
		double readImageTime_max;
		double blobMaskTime;
		double blobMaskTime_min;
		double blobMaskTime_max;
		double blobMaskTime_sum;
		double blobFindingTime;
		double blobFindingTime_min;
		double blobFindingTime_max;
		double blobFindingTime_sum;
		double compactSourceTime;
		double compactSourceTime_min;
		double compactSourceTime_max;
		double compactSourceTime_sum;
		double sourceSelectionTime;
		double sourceSelectionTime_min;
		double sourceSelectionTime_max;
		double sourceSelectionTime_sum;
		double imgResidualTime;
		double imgResidualTime_min;
		double imgResidualTime_max;
		double imgResidualTime_sum;
		double extendedSourceTime;
		double extendedSourceTime_min;
		double extendedSourceTime_max;
		double extendedSourceTime_sum;
		double sourceFitTime;
		double sourceFitTime_min;
		double sourceFitTime_max;
		double sourceFitTime_sum;
		double saveTime;
		double virtMemPeak;	
		double virtMemPeak_min;
		double virtMemPeak_max;

		//Source
		Source* m_Source;
		std::vector<Source*> m_SourcesMergedAtEdges;
		std::vector<Source*> m_SourceCollection;

		//Read options
		bool m_ReadTile;
		double m_TileMinX;
		double m_TileMaxX;
		double m_TileMinY;
		double m_TileMaxY;
		float m_ImgXmin;
		float m_ImgYmin;

		//Image beam options
		double m_beamBmaj;
		double m_beamBmin;
		double m_beamBpa;
		double m_pixSizeX;
		double m_pixSizeY;
		double m_beamFWHM;
		double m_beamFWHMMin;
		double m_beamFWHMMax;
		double m_pixSize;
		double m_beamTheta;
		double m_fluxCorrectionFactor;

		//Read distributed options
		bool m_splitInTiles;
		long int m_TileSizeX;
		long int m_TileSizeY;
		bool m_UseTileOverlap;
		double m_TileStepSizeX;
		double m_TileStepSizeY;
		bool m_mergeSourcesAtEdge;
		bool m_mergeSources;
		bool m_mergeExtendedSources;
		bool m_mergeCompactSources;	

		//Stat computation
		bool m_useParallelMedianAlgo;

		//Bkg computation
		ImgBkgData* m_BkgData;
		Image* m_SignificanceMap;
		bool m_UseLocalBkg;	
		bool m_Use2ndPassInLocalBkg;
		bool m_SkipOutliersInLocalBkg;
		int m_LocalBkgMethod;
		int m_BkgEstimator;
		bool m_UseBeamInfoInBkg;
		double m_BoxSizeX;
		double m_BoxSizeY;
		double m_GridSizeX;
		double m_GridSizeY;


		//Residual map
		Image* m_ResidualImg;
		ImgBkgData* m_ResidualBkgData;
		//double m_DilateZBrightThr;
		double m_residualZHighThr;
		//double m_DilateZThr;
		double m_residualZThr;
		//bool m_DilateNestedSources;
		bool m_removeNestedSources;
		int m_dilateKernelSize;
		//int m_DilatedSourceType;		
		int m_removedSourceType;
		//int m_DilateSourceModel;
		int m_residualModel;
		//bool m_DilateRandomize;
		bool m_residualModelRandomize;
		int m_psSubtractionMethod;
		bool m_UseResidualInExtendedSearch;

		//Smoothing
		bool m_UsePreSmoothing;	
		int m_SmoothFilter;			
		int m_GausFilterKernSize;
		double m_GausFilterSigma;
		double m_GuidedFilterRadius;
		double m_GuidedFilterColorEps;

		//Compact source search
		bool m_SearchCompactSources;
		int m_NMinPix;
		double m_SeedThr;
		double m_MergeThr;
		bool m_MergeBelowSeed;
		bool m_SearchNegativeExcess;
		int m_compactSourceSearchNIters;
		double m_seedThrStep;

		//Nested source search
		Image* m_blobMask;
		int m_blobMaskMethod;
		bool m_SearchNestedSources;
		double m_SourceToBeamAreaThrToSearchNested;
		double m_NestedBlobThrFactor;
		double m_minNestedMotherDist;
		double m_maxMatchingPixFraction;
		double m_nestedBlobPeakZThr;
		double m_nestedBlobPeakZMergeThr;
		double m_nestedBlobMinScale;
		double m_nestedBlobMaxScale;
		double m_nestedBlobScaleStep;
		double m_nestedBlobKernFactor;
		
		//Source selection
		bool m_ApplySourceSelection;
		bool m_useMinBoundingBoxCut;
		double m_SourceMinBoundingBox;
		bool m_useCircRatioCut;
		double m_psCircRatioThr;
		bool m_useElongCut;
		double m_psElongThr;
		bool m_useEllipseAreaRatioCut;
		double m_psEllipseAreaRatioMinThr;
		double m_psEllipseAreaRatioMaxThr;
		bool m_useMaxNPixCut;
		double m_psMaxNPix;
		bool m_useNBeamsCut;
		double m_psNBeamsThr;

		//Source fitting
		bool m_fitSources;
		double m_nBeamsMaxToFit;
		int m_fitMaxNComponents;
		bool m_fitWithCentroidLimits;
		bool m_fixCentroidInPreFit;
		double m_fitCentroidLimit;
		bool m_fitWithFixedBkg;
		bool m_fitWithBkgLimits;
		double m_fitBkgLevel;	
		bool m_fitUseEstimatedBkgLevel;
		bool m_fitWithAmplLimits;
		bool m_fixAmplInPreFit;
		double m_fitAmplLimit;
		bool m_fixSigmaInPreFit;
		bool m_fitWithSigmaLimits;
		double m_fitSigmaLimit;
		bool m_fitWithFixedSigma;
		bool m_fitWithFixedTheta;
		bool m_fitWithThetaLimits;
		bool m_fixThetaInPreFit;
		double m_fitThetaLimit;
		bool m_useFluxZCutInFit;
		double m_fitZCutMin;
		int m_peakMinKernelSize;
		int m_peakMaxKernelSize;
		int m_peakKernelMultiplicityThr;
		int m_peakShiftTolerance;
		double m_peakZThrMin;
		double m_fitFcnTolerance;
		long int m_fitMaxIters;
		bool m_fitImproveConvergence;
		long int m_fitNRetries;
		bool m_fitDoFinalMinimizerStep;
		int m_fitFinalMinimizer;
		bool m_fitUseNestedAsComponents;
		double m_fitChi2RegPar;

		std::string m_fitMinimizer;		
		std::string m_fitMinimizerAlgo;
		int m_fitStrategy;
		int m_fitPrintLevel;
		double m_fitParBoundIncreaseStepSize;
		bool m_fitUseThreads;

		//Saliency computation
		Image* m_SaliencyImg;
		bool m_SaliencyUseOptimalThr;
		double m_SaliencyThrFactor;
		double m_SaliencyBkgThrFactor;
		double m_SaliencyImgThrFactor;
		int m_SaliencyResoMin;
		int m_SaliencyResoMax;
		int m_SaliencyResoStep;
		bool m_SaliencyUseRobustPars;
		bool m_SaliencyUseBkgMap;
		bool m_SaliencyUseNoiseMap;
		//bool m_SaliencyUseCurvInDiss;
		double m_SaliencyNNFactor;
		//double m_SaliencySpatialRegFactor;
		double m_SaliencyMultiResoCombThrFactor;
		double m_SaliencyDissExpFalloffPar;
		double m_SaliencySpatialDistRegPar;

		//Extended sources
		bool m_SearchExtendedSources;
		int m_ExtendedSearchMethod;

		//Wavelet Transform
		int m_wtScaleSearchMin;
		int m_wtScaleSearchMax;
		
		
		//Superpixel options
		int m_spSize;
		double m_spBeta;
		int m_spMinArea;
		bool m_spUseLogContrast;

		//Active-Contour main options	
		int m_acMethod;
		long int m_acNIters;
		int m_acInitLevelSetMethod;
		float m_acInitLevelSetSizePar;
		double m_acTolerance;
		
		//Chan-Vese options
		//int m_cvNIters;
		long int m_cvNItersInner;
		long int m_cvNItersReInit;
		double m_cvTimeStepPar;
		double m_cvWindowSizePar;
		double m_cvLambda1Par;
		double m_cvLambda2Par;
		double m_cvMuPar;
		double m_cvNuPar;
		double m_cvPPar;
		//bool m_cvInitContourToSaliencyMap;
		
		//LRAC options
		//int m_lracNIters;
		double m_lracLambdaPar;
		double m_lracRadiusPar;
		double m_lracEpsPar;
		//bool m_lracInitContourToSaliencyMap;
		
		//Hierachical clustering data
		Image* m_LaplImg;
		Image* m_EdgeImg;
		Image* m_SegmImg;
		int m_spMergingNSegmentsToStop;
		double m_spMergingRatio;
		double m_spMergingRegPar;
		double m_spMergingMaxDissRatio;
		double m_spMergingMaxDissRatio2ndNeighbours;
		double m_spMergingDissThreshold;
		int m_spMergingEdgeModel;		
		bool m_spMergingUse2ndNeighbours;
		bool m_spMergingIncludeSpatialPars;
		bool m_spMergingAddCurvDist;
		bool m_spMergingUseRobustPars;
		
		//MPI vars
		int m_nProc;
		int m_procId;
		int m_workerRanks;
		int m_nWorkers;
		bool m_mpiEnabled;
		bool m_mpiGroupsInitialized;
		#ifdef MPI_ENABLED
			MPI_Group m_WorldGroup;
			MPI_Group m_WorkerGroup;			
			MPI_Comm m_WorkerComm;
		#endif	

		//Task data
		std::vector< std::vector<TaskData*> > m_taskDataPerWorkers;

		//Task info tree
		TTree* m_TaskInfoTree;
		double m_xmin;
		double m_xmax;
		double m_ymin;
		double m_ymax;


	ClassDef(SFinder,1)

};//close SourceFinder


#ifdef __MAKECINT__
#pragma link C++ class SFinder+;
#endif

}//close namespace

#endif

