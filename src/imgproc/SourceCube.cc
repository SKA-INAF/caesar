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
* @file SourceCube.cc
* @class SourceCube
* @brief Source cube class
*
* Class representing an image source cube
* @author S. Riggi
* @date 26/03/2018
*/

#include <SourceCube.h>
#include <Source.h>
#include <CodeUtils.h>
#include <GraphicsUtils.h>

#include <Image.h>
#include <Contour.h>

#include <TObject.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TPad.h>
#include <TExec.h>
#include <TStyle.h>

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

ClassImp(Caesar::SourceCube)

namespace Caesar {


SourceCube::SourceCube() 
	: TNamed()
{
	//Initialize 
	Init();

}//close costructor

SourceCube::SourceCube(std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Init pars
	Init();

}//close constructor

SourceCube::~SourceCube()
{
	//Delete source collection
	DEBUG_LOG("Deleting source collection added in cube...");
	CodeUtils::DeletePtrCollection<Source>(m_sources);
	DEBUG_LOG("done!");

	DEBUG_LOG("Deleting source plot canvas...");
	CodeUtils::DeletePtr<TCanvas>(m_sourcePlot);
	DEBUG_LOG("done!");

	DEBUG_LOG("Deleting source SED plot canvas...");
	CodeUtils::DeletePtr<TCanvas>(m_sourceSEDPlot);
	DEBUG_LOG("done!");
	
}//close destructor


SourceCube::SourceCube(const SourceCube& sourceCube) 
{
  // Copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((SourceCube&)sourceCube).Copy(*this);
}

void SourceCube::Copy(TObject &obj) const 
{
	DEBUG_LOG("Copying parent TNamed...");
	TNamed::Copy((SourceCube&)obj);

	// Delete and copy sourceplot canvas
	/*
	if(((SourceCube&)obj).m_sourcePlot){
		delete ((SourceCube&)obj).m_sourcePlot;
		((SourceCube&)obj).m_sourcePlot= 0;
	}
	if(m_sourcePlot){
		((SourceCube&)obj).m_sourcePlot= new TCanvas;
		*((SourceCube&)obj).m_sourcePlot = *m_sourcePlot;
	}
	*/

	//Delete and copy source SED
	if(((SourceCube&)obj).m_sourceSED){
		delete ((SourceCube&)obj).m_sourceSED;
		((SourceCube&)obj).m_sourceSED= 0;
	}
	if(m_sourceSED){
		((SourceCube&)obj).m_sourceSED= new TGraphAsymmErrors;
		*((SourceCube&)obj).m_sourceSED = *m_sourceSED;
	}

	//Delete and copy source component SEDs
	for(size_t i=0;i<(((SourceCube&)obj).m_sourceComponentSED).size();i++){
		if( (((SourceCube&)obj).m_sourceComponentSED)[i] ){
			delete (((SourceCube&)obj).m_sourceComponentSED)[i];
			(((SourceCube&)obj).m_sourceComponentSED)[i]= 0;
		}
	}
	(((SourceCube&)obj).m_sourceComponentSED).clear();

	TGraphAsymmErrors* aSourceComponentSED= 0;
	for(size_t i=0;i<m_sourceComponentSED.size();i++){
		aSourceComponentSED= new TGraphAsymmErrors;
		*aSourceComponentSED= *(m_sourceComponentSED[i]);
		(((SourceCube&)obj).m_sourceComponentSED).push_back(aSourceComponentSED);
	}	
	

	//Copy source collection
	//Delete first any existing collection
	for(size_t i=0;i<(((SourceCube&)obj).m_sources).size();i++){
		if( (((SourceCube&)obj).m_sources)[i] ){
			delete (((SourceCube&)obj).m_sources)[i];
			(((SourceCube&)obj).m_sources)[i]= 0;
		}
	}
	(((SourceCube&)obj).m_sources).clear();

	Source* aSource= 0;
	for(size_t i=0;i<m_sources.size();i++){
		aSource= new Source;
		*aSource= *(m_sources[i]);
		(((SourceCube&)obj).m_sources).push_back(aSource);
	}

}//close Copy()

SourceCube& SourceCube::operator=(const SourceCube& sourceCube) 
{ 
	// Operator =
  if (this != &sourceCube) ((SourceCube&)sourceCube).Copy(*this);
  return *this;
}


void SourceCube::Init()
{
	//Init source cube data
	m_sources.clear();

	//Init source plot canvas
	m_sourcePlot= 0;
	m_sourceSEDPlot= 0;

	//Init source SEDs
	m_sourceSED= 0;
	m_sourceComponentSED.clear();

}//close Init()


int SourceCube::DoSourceImagePlot(bool useWCS,int coordSyst)
{
	//Delete previous canvas and primitives (if any)
	if(m_sourcePlot){
		int padCounter= 1;
		while(m_sourcePlot->GetPad(padCounter)) { 
			//Select current pad
			m_sourcePlot->cd(padCounter);

			//Retrieve image from pad primitives
			Image* simg= GraphicsUtils::FindImageFromPad();
			if(simg){
				INFO_LOG("Deleting source image retrieved as primitive from pad no. "<<padCounter);
				CodeUtils::DeletePtr<Image>(simg);
			}

			padCounter++; 	
		}//end loop pads

		//Delete canvas
		INFO_LOG("Deleting source plot canvas...");
		CodeUtils::DeletePtr<TCanvas>(m_sourcePlot);
	}//close if m_sourcePlot

	//Loop over source, get images and plot
	int nChannels= static_cast<int>(m_sources.size());
	int pixMargin= 30;
		
	struct SourceLabelInfo {
		SourceLabelInfo(){
			sname= "";
			freq= "";
		}
		std::string sname;
		std::string freq;
		TVector2 pos;
	};	
	
	//Draw canvas
	TString canvasName= Form("Plot_%s",this->GetName());
	m_sourcePlot= new TCanvas(canvasName,canvasName,1200,600);
	m_sourcePlot->cd();
	
	int nRows= 3; 
	int nCols= std::ceil((float)(nChannels)/(float)(nRows));
	m_sourcePlot->Divide(nRows,nCols);
	int lineColor= kBlack;
	int lineStyle= kSolid;

	for(int i=0;i<nChannels;i++){
		//Select pad where to draw
		m_sourcePlot->cd(i+1);

		//Set style
		gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadBottomMargin(0.1);
  	gStyle->SetPadLeftMargin(0.15);
		gStyle->SetPadRightMargin(0.15);
		
		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.15);
  	gPad->SetRightMargin(0.15);

		//Draw source image
		bool drawImg= true;
		bool drawContours= true;
		bool drawNested= true;
		bool drawFitComponents= true;
		m_sources[i]->Draw(pixMargin,eFluxMap,drawImg,drawContours,drawNested,drawFitComponents,lineColor,lineStyle,useWCS,coordSyst);

		/*
		//Get source image
		Image* simg= m_sources[i]->GetImage(eFluxMap,pixMargin);
		if(!simg){
			ERROR_LOG("Failed to get source image for cube channel "<<i+1<<"!");
			return nullptr;
		}
	
		//Get source info for label
		ImgMetaData* metadata= m_sources[i]->GetImageMetaData();

		SourceLabelInfo labelInfo;
		labelInfo.sname= std::string(m_sources[i]->GetName());
		std::string units= "";
		if(metadata){
			TString freqStr= Form("%f %s",metadata->Freq,(metadata->FreqUnit).c_str());
			labelInfo.freq= std::string(freqStr);
			units= metadata->BUnit;
		}


		//Select pad where to draw
		Plot->cd(i+1);

		//Set style
		gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadBottomMargin(0.1);
  	gStyle->SetPadLeftMargin(0.15);
		gStyle->SetPadRightMargin(0.15);
		
		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.15);
  	gPad->SetRightMargin(0.15);

		//Draw image
		TString histoName= Form("simg_ch%d",i+1);
		TH2D* htemp= simg->GetHisto2D(std::string(histoName));
		htemp->SetStats(0);
		htemp->GetXaxis()->SetTitle("X");
		htemp->GetYaxis()->SetTitle("Y");
		htemp->GetZaxis()->SetTitle(units.c_str());
		htemp->GetZaxis()->SetTitleSize(0.05);
		htemp->GetZaxis()->SetTitleOffset(0.9);
		htemp->Draw();
		gPad->Update();
		if(useWCS) htemp->Draw("COLAZ");
		else htemp->Draw("COLZ");

		//Append Caesar image to current pad
		simg->AppendPad();

		//Set WCS axis	
		if(useWCS){
			gPad->Update();
			TGaxis* xaxis_wcs= new TGaxis;
			TGaxis* yaxis_wcs= new TGaxis;
			bool useImageCoords= false;
			int status= GraphicsUtils::SetWCSAxis(simg,*xaxis_wcs,*yaxis_wcs,coordSyst,useImageCoords);
			if(status>=0){
				TExec* ex = new TExec("ex","GraphicsUtils::PadUpdater_PhysCoords()");
   			htemp->GetListOfFunctions()->Add(ex);
				xaxis_wcs->Draw("same");
				yaxis_wcs->Draw("same");
			}
			else{
				WARN_LOG("Failed to set gAxis!");
			}	
		}//close if useWCS
		*/
		
	}//end loop image to be drawn

	return 0;

}//close DrawSourceImages()

int SourceCube::DoSourceSEDs()
{
	//Delete existing source SEDs
	CodeUtils::DeletePtr<TGraphAsymmErrors>(m_sourceSED);
	CodeUtils::DeletePtrCollection<TGraphAsymmErrors>(m_sourceComponentSED);

	//Delete existing source SED plot
	CodeUtils::DeletePtr<TCanvas>(m_sourceSEDPlot);

	//Fill source SEDs
	int nChannels= static_cast<int>(m_sources.size());
	m_sourceSED= new TGraphAsymmErrors(nChannels);

	double fluxMin= 1.e+99;
	double fluxMax= -1.e-99;
	double nuMin= 1.e+99;
	double nuMax= -1.e-99;

	for(int i=0;i<nChannels;i++){
		//Get source metadata
		ImgMetaData* metadata= m_sources[i]->GetImageMetaData();
		if(!metadata){
			WARN_LOG("Failed to get metadata for source channel no. "<<i+1<<"!");
			return -1;
		}
	
		//Get source frequency
		double Nu= metadata->Freq;
		double dNu= metadata->dFreq;
		std::string FreqUnits= metadata->FreqUnit;
		if(FreqUnits=="Hz"){//convert to GHz
			Nu/= 1.e+9;
			dNu/= 1.e+9;
		}
		double lgNu= log10(Nu);
		double dlgNu= log10(TMath::E())*dNu/Nu;
		if(Nu<nuMin) nuMin= Nu;
		if(Nu>nuMax) nuMax= Nu;
		
		//Get source flux
		double flux= 0;
		double fluxErr= 0;
		if(m_sources[i]->GetFluxDensity(flux)<0){
			WARN_LOG("Failed to get flux density for source channel no. "<<i+1<<"!");
			return -1;
		}
		if(m_sources[i]->GetFluxDensityErr(fluxErr)<0){
			WARN_LOG("Failed to get flux density error for source channel no. "<<i+1<<"!");
			return -1;
		}
		double lgFlux= log10(flux);
		double lgFluxErr= log10(TMath::E())*fluxErr/flux;
		if(flux<fluxMin) fluxMin= flux;
		if(flux>fluxMax) fluxMax= flux;

		//Get estimated condon errors on components
		std::vector<double> componentFluxDensityErr_condon;
		if(m_sources[i]->GetCondonComponentFluxDensityErr(componentFluxDensityErr_condon)<0){
			WARN_LOG("Failed to estimated condon errors on components...");
		}

		INFO_LOG("Channel no. "<<i+1<<": nu="<<Nu<<", dnu="<<dNu<<", flux="<<flux<<", fluxErr="<<fluxErr<<", fluxErr(Condon)="<<componentFluxDensityErr_condon[0]);

		//Fill SED
		m_sourceSED->SetPoint(i,lgNu,lgFlux);
		m_sourceSED->SetPointError(i,dlgNu,dlgNu,lgFluxErr,lgFluxErr);
		//m_sourceSED->SetPoint(i,Nu,flux);
		//m_sourceSED->SetPointError(i,dNu,dNu,fluxErr,fluxErr);
		
	}//end loop channels
	
	//Get component SEDs
	TGraphAsymmErrors* componentSEDGraph= 0;

	for(size_t i=0;i<m_componentIndexes.size();i++){//loop over components matches
		
		componentSEDGraph= new TGraphAsymmErrors();

		for(size_t j=0;j<m_componentIndexes[i].size();j++){//loop over channels
			size_t sindex= m_componentIndexes[i][j].first;
			size_t cindex= m_componentIndexes[i][j].second;
			INFO_LOG("Match component no. "<<i+1<<": sindex="<<sindex<<", cindex="<<cindex);
		
			//Get component flux
			SourceFitPars fitPars= m_sources[sindex]->GetFitPars();	
			double flux= fitPars.GetComponentFluxDensity(cindex); 
			double fluxErr= fitPars.GetComponentFluxDensityErr(cindex); 
			double lgFlux= log10(flux);
			double lgFluxErr= log10(TMath::E())*fluxErr/flux;

			//Get source metadata
			ImgMetaData* metadata= m_sources[sindex]->GetImageMetaData();
			if(!metadata){
				WARN_LOG("Failed to get metadata for source channel no. "<<i+1<<"!");
				return -1;
			}

			//Get source frequency
			double Nu= metadata->Freq;
			double dNu= metadata->dFreq;
			std::string FreqUnits= metadata->FreqUnit;
			if(FreqUnits=="Hz"){//convert to GHz
				Nu/= 1.e+9;
				dNu/= 1.e+9;
			}
			double lgNu= log10(Nu);
			double dlgNu= log10(TMath::E())*dNu/Nu;

			//Fill SED
			componentSEDGraph->SetPoint(j,lgNu,lgFlux);
			componentSEDGraph->SetPointError(j,dlgNu,dlgNu,lgFluxErr,lgFluxErr);
			//componentSEDGraph->SetPoint(j,Nu,flux);
			//componentSEDGraph->SetPointError(j,dNu,dNu,fluxErr,fluxErr);
			
		}//end loop channel for this component

		//Add component SED to list
		m_sourceComponentSED.push_back(componentSEDGraph);

	}//end loop components


	//## Fill canvas
	TString canvasName= Form("SEDPlot_%s",this->GetName());
	m_sourceSEDPlot= new TCanvas(canvasName,canvasName,800,800);
	m_sourceSEDPlot->cd();

	double lgNu_min= log10(nuMin-0.6*nuMin);
	double lgNu_max= log10(nuMax+0.6*nuMax);
	double lgFlux_min= log10(fluxMin-0.6*fluxMin);
	double lgFlux_max= log10(fluxMax+0.6*fluxMax);

	TH2D* PlotBkg= new TH2D("PlotBkg","",100,lgNu_min,lgNu_max,100,lgFlux_min,lgFlux_max);
	PlotBkg->GetXaxis()->SetTitle("log10(#nu)");
	PlotBkg->GetYaxis()->SetTitle("log10(S)");	
	PlotBkg->SetStats(0);
	PlotBkg->Draw();

	//Draw source SED
	m_sourceSED->SetMarkerStyle(8);
	m_sourceSED->SetMarkerColor(kBlack);
	m_sourceSED->SetLineColor(kBlack);
	m_sourceSED->Draw("epsame");
	//m_sourceSED->Draw("AP");

	//Draw source component SEDs
	if(m_sourceComponentSED.size()>1){
		std::vector<int> componentColors {kRed,kGreen+1,kBlue,kMagenta,kYellow+1,kOrange,kGray};
		for(size_t i=0;i<m_sourceComponentSED.size();i++){	
			int color= kBlack;
			if(m_sourceComponentSED.size()<componentColors.size()){
				color= componentColors[i];
			}
			m_sourceComponentSED[i]->SetMarkerStyle(21);
			m_sourceComponentSED[i]->SetMarkerColor(color);
			m_sourceComponentSED[i]->SetLineColor(color);
			m_sourceComponentSED[i]->Draw("epsame");
		}
	}
	
	return 0;

}//close DoSourceSEDs()

}//close namespace 


