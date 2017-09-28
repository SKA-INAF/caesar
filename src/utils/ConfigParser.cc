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
* @file ConfigParser.cc
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#include <ConfigParser.h>
#include <Option.h>
#include <SysUtils.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

ClassImp(Caesar::ConfigParser)

namespace Caesar {

bool ConfigParser::m_HasRegisteredOptions;
std::string ConfigParser::m_ConfigFile= "";

ConfigParser::ConfigParser() {

	m_ConfigFile= "";
	m_HasRegisteredOptions= false;	

}//close costructor

ConfigParser::~ConfigParser(){

}//close destructor

int ConfigParser::Parse(std::string filename){
	
	//## Check and read file
	if(filename==""){
		cerr<<"ConfigParser::Parse(): ERROR: Empty filename string given!"<<endl;
		return -1;
	}

	Caesar::FileInfo info;
	if(!SysUtils::CheckFile(filename,info,false)){
		cerr<<"ConfigParser::Parse(): ERROR: Invalid config file specified (check if file actually exist)!"<<endl;
		return -1;
	}

	//## Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		cerr<<"ConfigParser::Parse(): ERROR: Failed to open config file "<<filename<<" for reading!"<<endl;
		return -1;
  }

	//## Parsing file
	m_ConfigFile= filename;
	cout<<"ConfigParser::Parse(): INFO: Reading and parsing file: "<<filename<<endl;

	std::string parsedline= "";

	while(std::getline(in,parsedline)) {
		//Check file
		if (!in.good()) break;

		//Check first character
		char first_char= *(parsedline.c_str());
		if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){//skip line
			continue;
		}
		
		istringstream line(parsedline);

		//Get option name
		std::string optionKey= "";    
		line >> optionKey;
    if(optionKey=="\n" || optionKey=="") continue;

		//Get option delimiter
		std::string optionDelimiter= "";   
		line >> optionDelimiter;
		if(optionDelimiter!="=") continue;

		//Get and set option value
		std::string optionValue= "";
		line >> optionValue;
		
		if(SetOptionFromString(optionKey,optionValue)<0){
			cerr<<"ConfigParser::Parse(): WARN: Failed to set option "<<optionKey<<" (unregistered?)...skip it!"<<endl;
			continue;
		}
		//cout<<"ConfigParser::Parse(): INFO: optionName: "<<optionKey<<", val="<<optionValue<<endl;

		if (!in.good()) break;
		
	}//close while

	//## Close file
	in.close();

	return 0;

}//close Parse()


int ConfigParser::RegisterPredefinedOptions(){

	cout<<"ConfigParser::RegisterPredefinedOptions(): INFO: Registering pre-defined options..."<<endl;

	//Register pre-defined options
	try {
		//=====================
		//==  Main options   ==
		//=====================
		REGISTER_OPTION(inputFile,std::string,"","","");
		REGISTER_OPTION(inputImage,std::string,"img","","");
		REGISTER_OPTION(fitsHDUId,int,1,1,1000);
		REGISTER_OPTION(outputFile,std::string,"Output.root","","");
		REGISTER_OPTION(ds9RegionFile,std::string,"DS9SourceRegion.reg","","");
		REGISTER_OPTION(saveDS9Region,bool,true,false,true);
		REGISTER_OPTION(DS9RegionFormat,int,2,0,3);
		REGISTER_OPTION(saveToFile,bool,true,false,true);
		REGISTER_OPTION(saveConfig,bool,true,false,true);
		REGISTER_OPTION(saveResidualMap,bool,true,false,true);
		REGISTER_OPTION(saveBkgMap,bool,true,false,true);
		REGISTER_OPTION(saveNoiseMap,bool,true,false,true);
		REGISTER_OPTION(saveSignificanceMap,bool,true,false,true);
		REGISTER_OPTION(saveInputMap,bool,true,false,true);
		REGISTER_OPTION(saveSaliencyMap,bool,true,false,true);
		REGISTER_OPTION(saveSources,bool,true,false,true);
		REGISTER_OPTION(saveEdgenessMap,bool,false,false,true);
		REGISTER_OPTION(saveCurvatureMap,bool,false,false,true);
		REGISTER_OPTION(saveSegmentedMap,bool,true,false,true);
		REGISTER_OPTION(saveToFITSFile,bool,false,false,true);
		REGISTER_OPTION(residualMapFITSFile,std::string,"ResidualMap.fits","","");
		REGISTER_OPTION(inputMapFITSFile,std::string,"InputMap.fits","","");
		REGISTER_OPTION(saliencyMapFITSFile,std::string,"SaliencyMap.fits","","");
		REGISTER_OPTION(bkgMapFITSFile,std::string,"BkgMap.fits","","");
		REGISTER_OPTION(noiseMapFITSFile,std::string,"NoiseMap.fits","","");
		REGISTER_OPTION(significanceMapFITSFile,std::string,"SignificanceMap.fits","","");
		REGISTER_OPTION(isInteractiveRun,bool,false,false,true);		

		REGISTER_OPTION(readTileImage,bool,false,false,true);
		REGISTER_OPTION(tileMinX,double,0,-1.e+6,1.e+6);
		REGISTER_OPTION(tileMaxX,double,0,-1.e+6,1.e+6);
		REGISTER_OPTION(tileMinY,double,0,-1.e+6,1.e+6);
		REGISTER_OPTION(tileMaxY,double,0,-1.e+6,1.e+6);


		//=======================================
		//==  Image beam options               ==
		//=======================================
		REGISTER_OPTION(beamFWHM,double,6.5,0,3600);//circular beam FWHM in arcsec
		REGISTER_OPTION(beamTheta,double,0,-0.0001,180);//circular beam rot in degrees
		REGISTER_OPTION(pixSize,double,1,0,3600);//pixel area in arcsec

		//=======================================
		//==  Distributed processing options   ==
		//=======================================
		REGISTER_OPTION(splitInTiles,bool,false,false,true);
		REGISTER_OPTION(nThreads,int,-1,-1,100);
		REGISTER_OPTION(tileSizeX,long int,1000,1,10000000);
		REGISTER_OPTION(tileSizeY,long int,1000,1,10000000);	
		REGISTER_OPTION(useTileOverlap,bool,false,false,true);
		REGISTER_OPTION(tileStepSizeX,double,1,0.001,1);
		REGISTER_OPTION(tileStepSizeY,double,1,0.001,1);
		REGISTER_OPTION(mergeSourcesAtEdge,bool,true,false,true);
		

		//======================
		//==  Logger options  ==
		//======================
		REGISTER_OPTION(loggerTarget,int,1,1,3);//1=CONSOLE, 2=FILE, 3=SYSLOG
		REGISTER_OPTION(loggerTag,std::string,"logger","","");
		REGISTER_OPTION(logLevel,std::string,"INFO","","");//INFO, WARN, DEBUG, WARN, ERROR, FATAL
		REGISTER_OPTION(logFile,std::string,"out.log","","");//only for FILE TARGET
		REGISTER_OPTION(appendToLogFile,bool,false,false,true);//only for FILE TARGET
		REGISTER_OPTION(maxLogFileSize,std::string,"10MB","","");//only for FILE TARGET
		REGISTER_OPTION(maxBackupLogFiles,int,1,1,1000000);//only for FILE TARGET
		REGISTER_OPTION(consoleTarget,std::string,"System.out","","");//only for CONSOLE TARGET (System.out, System.err)
		REGISTER_OPTION(syslogFacility,std::string,"local6","","");//only for SYSLOG TARGET (local0-7, user, syslog, ...)
		
		
		//===================
		//==  Bkg options  ==
		//===================
		REGISTER_OPTION(useLocalBkg,bool,true,false,true);
		REGISTER_OPTION(use2ndPassInLocalBkg,bool,true,false,true);
		REGISTER_OPTION(skipOutliersInLocalBkg,bool,false,false,true);
		REGISTER_OPTION(localBkgMethod,int,1,1,2);
		REGISTER_OPTION(bkgEstimator,int,2,1,4);
		REGISTER_OPTION(useBeamInfoInBkg,bool,true,false,true);
		REGISTER_OPTION(boxSizeX,double,20,0.01,1000);
		REGISTER_OPTION(boxSizeY,double,20,0.01,1000);
		REGISTER_OPTION(gridSizeX,double,0.2,0.,1.);
		REGISTER_OPTION(gridSizeY,double,0.2,0.,1.);
	
		//==============================
		//==  Filtering options       ==
		//==============================
		REGISTER_OPTION(usePreSmoothing,bool,true,false,true);
		REGISTER_OPTION(smoothFilter,int,2,1,2);
		REGISTER_OPTION(gausFilterKernSize,int,5,1,10001);
		REGISTER_OPTION(gausFilterSigma,double,1,0,1000);
		REGISTER_OPTION(guidedFilterRadius,double,12,0,1000);
		REGISTER_OPTION(guidedFilterColorEps,double,0.04,0,1000);
		

		//==============================
		//==  Source finding options  ==
		//==============================
		REGISTER_OPTION(searchCompactSources,bool,true,false,true);
		//REGISTER_OPTION(searchBrightSources,bool,false,false,true);
		//REGISTER_OPTION(searchFaintSources,bool,false,false,true);
		REGISTER_OPTION(searchNestedSources,bool,true,false,true);
		REGISTER_OPTION(minNPix,int,10,0,10000);
		REGISTER_OPTION(seedBrightThr,double,10,0,10000);	
		REGISTER_OPTION(seedThr,double,5,0,10000);
		REGISTER_OPTION(mergeThr,double,2.6,0,10000);
		REGISTER_OPTION(mergeBelowSeed,bool,false,false,true);
		REGISTER_OPTION(searchNegativeExcess,bool,false,false,true);
		REGISTER_OPTION(compactSourceSearchNIters,int,10,0,100);
		REGISTER_OPTION(seedThrStep,double,1,0,10);
		//REGISTER_OPTION(wtScaleFaint,int,1,1,10);

		//=====================================
		//==  Nested Source finding options  ==
		//=====================================
		REGISTER_OPTION(nestedBlobThrFactor,double,1,0,100);
		REGISTER_OPTION(minNestedMotherDist,double,2,0,100);
		REGISTER_OPTION(maxMatchingPixFraction,double,0.5,0,1);
		
		//=======================================
		//==  Extended Source finding options  ==
		//=======================================
		REGISTER_OPTION(searchExtendedSources,bool,false,false,true);
		REGISTER_OPTION(extendedSearchMethod,int,2,0,10);		
		REGISTER_OPTION(useResidualInExtendedSearch,bool,true,false,true);
		REGISTER_OPTION(wtScaleExtended,int,6,1,10);
		REGISTER_OPTION(activeContourMethod,int,2,0,3);		
		
		//================================
		//==  Source selection options  ==
		//================================
		REGISTER_OPTION(applySourceSelection,bool,true,false,true);
		REGISTER_OPTION(sourceMinBoundingBox,double,2,0,1000000);	
		REGISTER_OPTION(useMinBoundingBoxCut,bool,false,false,true);
		REGISTER_OPTION(useCircRatioCut,bool,false,false,true);
		REGISTER_OPTION(psCircRatioThr,double,0.4,0,1);
		REGISTER_OPTION(useElongCut,bool,false,false,true);
		REGISTER_OPTION(psElongThr,double,0.7,0,1);
		REGISTER_OPTION(useEllipseAreaRatioCut,bool,false,false,true);
		REGISTER_OPTION(psEllipseAreaRatioMinThr,double,0.6,0,10);
		REGISTER_OPTION(psEllipseAreaRatioMaxThr,double,1.4,0,10);
		REGISTER_OPTION(useMaxNPixCut,bool,false,false,true);
		REGISTER_OPTION(psMaxNPix,double,1000,0,1.e+7);
				
		//================================
		//==  Source residual options   ==
		//================================
		REGISTER_OPTION(dilateNestedSources,bool,true,false,true);		
		REGISTER_OPTION(dilateKernelSize,int,9,1,1001);		
		REGISTER_OPTION(dilatedSourceType,int,2,-1,3);
		REGISTER_OPTION(dilateSourceModel,int,1,1,2);
		REGISTER_OPTION(dilateRandomize,bool,false,false,true);

		//==================================
		//==  Source deblending options   ==
		//==================================
		REGISTER_OPTION(fitSources,bool,false,false,true);
		REGISTER_OPTION(fitMaxNComponents,int,3,0,100);			
		//REGISTER_OPTION(deblendCurvThr,double,0,-0.0001,1.00001);	
		//REGISTER_OPTION(deblendComponentMinNPix,double,5,0,100000);			
			
		
		//===================================
		//==  Chan-Vese algorithm options  ==
		//===================================
		REGISTER_OPTION(cvNIters,int,1000,0,100000);
		REGISTER_OPTION(cvTimeStepPar,double,0.007,0,1000);
		REGISTER_OPTION(cvWindowSizePar,double,1,0,1000);
		REGISTER_OPTION(cvLambda1Par,double,1,0,100);
		REGISTER_OPTION(cvLambda2Par,double,2,0,100);
		REGISTER_OPTION(cvMuPar,double,0.5,0,100);
		REGISTER_OPTION(cvNuPar,double,0,0,100);
		REGISTER_OPTION(cvPPar,double,1,0,100);

		//===================================================================
		//==  Linear Region-based Active Contour (LRAC) algorithm options  ==
		//===================================================================
		REGISTER_OPTION(lracNIters,int,1000,0,100000);
		REGISTER_OPTION(lracLambdaPar,double,0.1,0,1000);
		REGISTER_OPTION(lracRadiusPar,double,10,0,1000);
		REGISTER_OPTION(lracEpsPar,double,0.01,0,1000);//Convergence par

		//===================================
		//==  Saliency filtering options   ==
		//===================================
		REGISTER_OPTION(saliencyThrFactor,double,2.8,0,10);
		REGISTER_OPTION(saliencyBkgThrFactor,double,1,0.01,10);
		REGISTER_OPTION(saliencyImgThrFactor,double,1,0.01,10);
		REGISTER_OPTION(saliencyResoMin,int,20,1,1000);
		REGISTER_OPTION(saliencyResoMax,int,60,1,1000);
		REGISTER_OPTION(saliencyResoStep,int,10,1,100);
		REGISTER_OPTION(saliencyUseRobustPars,bool,true,false,true);
		REGISTER_OPTION(saliencyUseBkgMap,bool,true,false,true);
		REGISTER_OPTION(saliencyUseNoiseMap,bool,true,false,true);
		REGISTER_OPTION(saliencyUseCurvInDiss,bool,true,false,true);
		REGISTER_OPTION(saliencyNNFactor,double,0.2,0,1.0001);
		REGISTER_OPTION(saliencySpatialRegFactor,double,6,0,100);
		REGISTER_OPTION(saliencyMultiResoCombThrFactor,double,0.7,0,1);
		REGISTER_OPTION(saliencyDissExpFalloffPar,double,100,1.e-6,1.e+6);
		REGISTER_OPTION(saliencySpatialDistRegPar,double,1,0,1.e+6);

		
		//==============================================
		//==  Superpixel Generation options           ==
		//==============================================
		REGISTER_OPTION(spSize,int,20,5,10000);
		REGISTER_OPTION(spBeta,double,1,1.e-10,1.e+10);
		REGISTER_OPTION(spMinArea,int,10,1,10000);
		REGISTER_OPTION(spUseLogContrast,bool,false,false,true);

		//==============================================
		//==  Hierarchical clustering algo options    ==
		//==============================================
		REGISTER_OPTION(spMergingNSegmentsToStop,int,1,0,1000);
		REGISTER_OPTION(spMergingRatio,double,0.3,0,1.0001);
		REGISTER_OPTION(spMergingMaxDissRatio,double,1000,0,100000);
		REGISTER_OPTION(spMergingMaxDissRatio2ndNeighbours,double,1.05,0,100000);
		REGISTER_OPTION(spMergingDissThreshold,double,3,0,100000);
		REGISTER_OPTION(spMergingEdgeModel,int,2,1,2);
		REGISTER_OPTION(spMergingRegPar,double,0.5,-0.0001,1.0001);
		REGISTER_OPTION(spMergingIncludeSpatialPars,bool,false,false,true);
		REGISTER_OPTION(spMergingAddCurvDist,bool,true,false,true);
		REGISTER_OPTION(spMergingUseRobustPars,bool,false,false,true);
		REGISTER_OPTION(spMergingUse2ndNeighbours,bool,true,false,true);
		
		

		//Set has_registered flag (otherwise options are re-built)
		m_HasRegisteredOptions= true;

	}//close try block
	catch(...){
		cerr<<"ConfigParser::RegisterPredefinedOptions(): ERROR: Failed to load predefined options!"<<endl;
		return -1;
	}//close catch

	return 0;

}//close RegisterPredefinedOptions()

}//close namespace

