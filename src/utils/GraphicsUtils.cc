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
* @file GraphicsUtils.cc
* @class GraphicsUtils
* @brief Utility functions for graphics
*
* Utility functions for graphics
* @author S. Riggi
* @date 15/01/2016
*/


#include <GraphicsUtils.h>
#include <AstroUtils.h>
#include <Image.h>

#include <TObject.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TList.h>
#include <TF1.h>
#include <TPolyLine.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::GraphicsUtils)

namespace Caesar {

GraphicsUtils::GraphicsUtils(){

}

GraphicsUtils::~GraphicsUtils(){

}

void GraphicsUtils::SetPalette(int paletteStyle,int ncolors)
{

	switch(paletteStyle){
		case eRAINBOW :
			gStyle->SetPalette(55);
			break;
		case eBLACKWHITE :	
			SetBWPalette(ncolors);
			break;
		case eBLACKBODY :
			gStyle->SetPalette(53);
			break;
		case eHOT2COLD :
			SetHotColdPalette(ncolors);
			break;
		case eCOLD2HOT :
			SetColdHotPalette(ncolors);
			break;
		case eTHERMAL :
			SetThermalPalette(ncolors);
			break;
		default: 
			gStyle->SetPalette(55);
			break;
	}//close switch


}//close SetPalette()



int GraphicsUtils::SetWCSProjGrid(Image* img,std::vector<TPolyLine>& gridx,std::vector<TPolyLine>& gridy,int coordSystem){
	
	if(!img) {
		ERROR_LOG("Null ptr to image given!");
		return -1;
	}
	
	if(!gPad){
		ERROR_LOG("No pad available!");
		return -1;
	}

	//Get image ranges
	double xmin= gPad->GetUxmin();
	double xmax= gPad->GetUxmax();
	double ymin= gPad->GetUymin();
	double ymax= gPad->GetUymax();

	//Get image WCS
	if(!img->HasMetaData()){
		WARN_LOG("No meta-data present!");
		return -1;
	}
	ImgMetaData* metadata= img->GetMetaData();
	WorldCoor* wcs= metadata->GetWorldCoord(coordSystem);
	//WorldCoor* wcs= metadata->GetWorldCoord();
	if(!wcs){
		WARN_LOG("Cannot get WCS from image!");
		return -1;
	}		

	//Get range coord in WCS
	double xBR_wcs, xBL_wcs, xTR_wcs, xTL_wcs;
	double yBR_wcs, yBL_wcs, yTR_wcs, yTL_wcs;
	pix2wcs (wcs,xmin,ymin,&xBL_wcs,&yBL_wcs);//Bottom-Left corner
	pix2wcs (wcs,xmax,ymin,&xBR_wcs,&yBR_wcs);//Bottom-Right corner
	pix2wcs (wcs,xmin,ymax,&xTL_wcs,&yTL_wcs);//Top-Left corner
	pix2wcs (wcs,xmax,ymax,&xTR_wcs,&yTR_wcs);//Top-Right corner
	
	double xmin_wcs, xmax_wcs, ymin_wcs, ymax_wcs;	
	xmin_wcs= min(min(xBR_wcs,xBL_wcs),min(xTR_wcs,xTL_wcs));
	xmax_wcs= max(max(xBR_wcs,xBL_wcs),max(xTR_wcs,xTL_wcs));
	ymin_wcs= min(min(yBR_wcs,yBL_wcs),min(yTR_wcs,yTL_wcs));
	ymax_wcs= max(max(yBR_wcs,yBL_wcs),max(yTR_wcs,yTL_wcs));
	DEBUG_LOG("xmin="<<xmin<<" xmax="<<xmax<<" xmin_wcs="<<xmin_wcs<<" xmax_wcs="<<xmax_wcs<<" ymin="<<ymin<<" ymax="<<ymax<<" ymin_wcs="<<ymin_wcs<<" ymax_wcs="<<ymax_wcs);

	std::vector<double> x_list;
	std::vector<double> y_list;
	//double thisx= xmin_wcs;
	//double stepx= 0.1*fabs(xmax_wcs-xmin_wcs);
	//double thisy= ymin_wcs;
	//double stepy= 0.1*fabs(ymax_wcs-ymin_wcs);
	double thisx= xmin;
	double stepx= 0.1*fabs(xmax-xmin);
	double thisy= ymin;
	double stepy= 0.1*fabs(ymax-ymin);
	//while(thisx<xmax_wcs){
	while(thisx<xmax){
		x_list.push_back(thisx);
		cout<<"x="<<thisx<<endl;
		thisx+= stepx;
	}
	//while(thisy<ymax_wcs){
	while(thisy<ymax){
		y_list.push_back(thisy);
		cout<<"y="<<thisy<<endl;
		thisy+= stepy;
	}

		
	
	//Convert wcs to desired type
	char* flag = (char*)("");
	std::string wcsType= metadata->GetWCSType();
	if(coordSystem==eGALACTIC)
		flag = (char*)("GALACTIC");	
	else if(coordSystem==eJ2000)
		flag = (char*)("FK5");
	else if(coordSystem==eB1950)
		flag = (char*)("FK4");
	else if(coordSystem==-1 && wcsType!="")
		flag = (char*)(wcsType.c_str());
			
	if(strcmp(flag,"")!=0) {
		wcsininit (wcs,flag);			
	}
	
	
	//## Generate grid
	double stepx_wcs= 0.05*fabs(xmax_wcs-xmin_wcs);
	double stepy_wcs= 0.05*fabs(ymax_wcs-ymin_wcs);
	DEBUG_LOG("stepx_wcs="<<stepx_wcs<<" stepy_wcs="<<stepy_wcs);

	for(unsigned int i=0;i<x_list.size();i++){
		double xstart= x_list[i];
		double ystart= ymin;
		double xstart_wcs, ystart_wcs;
		//wcsoutinit (wcs,flag);
		pix2wcs (wcs,xstart,ystart,&xstart_wcs,&ystart_wcs); 
		DEBUG_LOG("Grid no. "<<i<<" (xstart,ystart)=("<<xstart<<","<<ystart<<") ==> ("<<xstart_wcs<<","<<ystart_wcs<<")");
	
		//Move along wcs proj
		int nPts= 0;
		double x= xstart;
		double y= ystart;
		double x_wcs= xstart_wcs;
		double y_wcs= ystart_wcs;
		int offset;
		TPolyLine thisPolyLine;

		while(img->HasBin(x,y)){	
			DEBUG_LOG("wcs("<<x_wcs<<","<<y_wcs<<") ==> ("<<x<<","<<y<<") offset="<<offset);
			thisPolyLine.SetPoint(nPts,x,y);	
			nPts++;
			
			x_wcs+= stepx_wcs;
			//wcsininit (wcs,flag);
			wcs2pix(wcs, x_wcs, y_wcs,&x,&y,&offset);
		}//end loop points
		gridx.push_back(thisPolyLine);

	}//end loop x
	

	for(unsigned int i=0;i<y_list.size();i++){
		double xstart= xmin;
		double ystart= y_list[i];
		double xstart_wcs, ystart_wcs;
		//wcsoutinit (wcs,flag);
		pix2wcs (wcs,xstart,ystart,&xstart_wcs,&ystart_wcs); 
		DEBUG_LOG("Grid no. "<<i<<" (xstart,ystart)=("<<xstart<<","<<ystart<<") ==> ("<<xstart_wcs<<","<<ystart_wcs<<")");
	
		//Move along wcs proj
		int nPts= 0;
		double x= xstart;
		double y= ystart;
		double x_wcs= xstart_wcs;
		double y_wcs= ystart_wcs;
		int offset;
		TPolyLine thisPolyLine;

		while(img->HasBin(x,y)){	
			DEBUG_LOG("wcs("<<x_wcs<<","<<y_wcs<<") ==> ("<<x<<","<<y<<") offset="<<offset);
			thisPolyLine.SetPoint(nPts,x,y);	
			nPts++;
			
			y_wcs+= stepy_wcs;
			//wcsininit (wcs,flag);
			wcs2pix(wcs, x_wcs, y_wcs,&x,&y,&offset);
		}//end loop points
		gridy.push_back(thisPolyLine);

	}//end loop x
	
	
	return 0;

}//close SetWCSProjGrid()


int GraphicsUtils::SetWCSAxis(Image* img,TGaxis& xaxis_wcs,TGaxis& yaxis_wcs,int coordSystem){

	if(!img) {
		ERROR_LOG("Null ptr to image given!");
		return -1;
	}
	
	//Get image ranges
	double xmin= gPad->GetUxmin();
	double xmax= gPad->GetUxmax();
	double ymin= gPad->GetUymin();
	double ymax= gPad->GetUymax();
	DEBUG_LOG("xmin/xmax="<<xmin<<"/"<<xmax<<" ymin/ymax="<<ymin<<"/"<<ymax);

	//Get image WCS
	if(!img->HasMetaData()){
		ERROR_LOG("No meta-data present!");
		return -1;
	}
	ImgMetaData* metadata= img->GetMetaData();
	WorldCoor* wcs= metadata->GetWorldCoord(coordSystem);
	if(!wcs){
		ERROR_LOG("Cannot get WCS from image!");
		return -1;
	}	
	std::string wcsType= metadata->GetWCSType();

	//Get range coord in WCS
	DEBUG_LOG("Set pixel2wcs coords...");
	double xmin_wcs, xmax_wcs, ymin_wcs, ymax_wcs;	
	AstroUtils::PixelToWCSCoords(img,wcs,xmin,ymin,xmin_wcs,ymin_wcs); 
	AstroUtils::PixelToWCSCoords(img,wcs,xmax,ymin,xmax_wcs,ymin_wcs);
	AstroUtils::PixelToWCSCoords(img,wcs,xmin,ymax,xmin_wcs,ymax_wcs);		
	double minx_wcs= min(xmin_wcs, xmax_wcs); 
	double maxx_wcs= max(xmin_wcs, xmax_wcs); 
	double miny_wcs= min(ymin_wcs, ymax_wcs); 
	double maxy_wcs= max(ymin_wcs, ymax_wcs);

	TF1* fcn= (TF1*)gROOT->GetFunction("invertXFcn");
	if(!fcn){
		fcn= new TF1("invertXFcn","-x");
	}
	fcn->SetRange(minx_wcs,maxx_wcs);
	
	//Clear up
	delete wcs;
	wcs= 0;

	//## Update pad and make axis
	gPad->Update();
	
	//--> xaxis	
	//double SMALL_VAL= 1.e-4;
	xaxis_wcs.SetName("xaxis_wcs");
	xaxis_wcs.SetTitleSize(0.055);
	xaxis_wcs.SetTitleOffset(0.8);
	
	xaxis_wcs.SetTextFont(52);
	xaxis_wcs.SetLabelSize(0.045);
	//xaxis_wcs.SetLabelSize(0.03);
	xaxis_wcs.SetLabelFont(42);
	//xaxis_wcs.SetNdivisions(505);
	//xaxis_wcs.SetNdivisions(50510);
	xaxis_wcs.SetNdivisions(510);
	xaxis_wcs.SetOption("+");
	xaxis_wcs.SetX1(gPad->GetUxmin());
	xaxis_wcs.SetX2(gPad->GetUxmax());
	xaxis_wcs.SetY1(gPad->GetUymin());
	xaxis_wcs.SetY2(gPad->GetUymin());	
	//xaxis_wcs.SetY1(gPad->GetUymin()-40);
	//xaxis_wcs.SetY2(gPad->GetUymin()-40);	
	xaxis_wcs.SetWmin(minx_wcs);
	xaxis_wcs.SetWmax(maxx_wcs);
	
	if(wcsType=="FK5" || wcsType=="FK4"){
		//xaxis_wcs.SetTitle("#alpha (deg)");
		xaxis_wcs.SetTitle("Right Ascension (deg)");

		//Invert axis
		xaxis_wcs.SetFunction("invertXFcn");
	}	
	else if(wcsType=="GALACTIC"){
		xaxis_wcs.SetTitle("l (deg)");	
	}
		
	//--> yaxis
	yaxis_wcs.SetName("yaxis_wcs");
	yaxis_wcs.SetTitleSize(0.055);
	yaxis_wcs.SetTitleOffset(1.45);
	yaxis_wcs.SetTextFont(52);
	//yaxis_wcs.SetLabelSize(0.03);
	yaxis_wcs.SetLabelSize(0.045);
	yaxis_wcs.SetLabelFont(42);
	yaxis_wcs.SetX1(gPad->GetUxmin());
	yaxis_wcs.SetX2(gPad->GetUxmin());
	//yaxis_wcs.SetX1(gPad->GetUxmin()-40);
	//yaxis_wcs.SetX2(gPad->GetUxmin()-40);
	
	yaxis_wcs.SetY1(gPad->GetUymin());
	yaxis_wcs.SetY2(gPad->GetUymax());
	yaxis_wcs.SetWmin(miny_wcs);
	yaxis_wcs.SetWmax(maxy_wcs);
	//yaxis_wcs.SetNdivisions(50510);
	//yaxis_wcs.SetNdivisions(50505);
	yaxis_wcs.SetNdivisions(510);
	yaxis_wcs.SetOption("-");
	if(wcsType=="FK5" || wcsType=="FK4"){
		//yaxis_wcs.SetTitle("#delta (deg)");
		yaxis_wcs.SetTitle("Declination (deg)");
	}	
	else if(wcsType=="GALACTIC"){
		yaxis_wcs.SetTitle("b (deg)");
	}	
	

	return 0;

}//close GetWCSAxis()

int GraphicsUtils::PadUpdater(){

	//## Check pad	
	if(!gPad){
		ERROR_LOG("No pad available!");
		return -1;
	}

	//## Retrieve current pad event and check it is an "axis" change event
	//int event = gPad->GetEvent();
  //int px = gPad->GetEventX();
  //int py = gPad->GetEventY();
	//cout<<"GraphicsUtils::UpdateGAxis(): INFO: event="<<event<<endl;
	//...
	//...

	//## Update gaxis if any
	if(UpdateGAxis()<0){
		WARN_LOG("Failed to update gAxis for current pad!");
	}

	return 0;

}//close PadUpdater()

Image* GraphicsUtils::FindImageFromPad(){

	//## Find image 
	TList* primitiveList= gPad->GetListOfPrimitives();
	if(!primitiveList){
		cerr<<"GraphicsUtils::FindImageFromPad(): WARN: Cannot retrieve the list of primitives!"<<endl;
		return 0;
	}

	Caesar::Image* img= 0;
	TString imgName= "";
	for(int i=0;i<primitiveList->GetSize();i++){
  	TObject* obj = (TObject*)primitiveList->At(i);
  	if(obj->ClassName() == std::string("Caesar::Image") ){
    	img = (Caesar::Image*)obj;
    	imgName = (img->GetName()).c_str();
			break;
  	}  	
	}//end loop primitives

	return img;

}//close FindImageFromPad()


int GraphicsUtils::UpdateGAxis(){

	/*
	//## Check pad	
	if(!gPad){
		cerr<<"GraphicsUtils::UpdateGAxis(): WARN: No pad available!"<<endl;
		return -1;
	}

	//## Retrieve current pad event and check it is an "axis" change event
	int event = gPad->GetEvent();
  int px = gPad->GetEventX();
  int py = gPad->GetEventY();
	//cout<<"GraphicsUtils::UpdateGAxis(): INFO: event="<<event<<endl;
	//...
	//...

	//## Find image 
	Caesar::Image* img= 0;
	TString imgName= "";
	TList* primitiveList= gPad->GetListOfPrimitives();
	if(!primitiveList){
		cerr<<"GraphicsUtils::UpdateGAxis(): WARN: Cannot retrieve the list of primitives!"<<endl;
		return -1;
	}

	for(int i=0;i<primitiveList->GetSize();i++){
  	TObject* obj = (TObject*)primitiveList->At(i);
  	if(obj->ClassName() == std::string("Caesar::Image") ){
    	img = (Caesar::Image*)obj;
    	imgName = img->GetName();
  	}
	}//end loop primitives
	*/

	//## Find gaxis
	TGaxis* xaxis_wcs= (TGaxis*)gPad->FindObject("xaxis_wcs");
	TGaxis* yaxis_wcs= (TGaxis*)gPad->FindObject("yaxis_wcs");
	if(!xaxis_wcs || !yaxis_wcs) {
		ERROR_LOG("Cannot get current gaxis!");
		return 0;
	}

	//## Find image from pad
	Image* img= FindImageFromPad();
	if(!img){
		ERROR_LOG("Cannot retrieve image from current pad!");
		return -1;
	}
	
	//## Get image WCS
	if(!img->HasMetaData()){
		ERROR_LOG("No meta-data present!");
		return -1;
	}
	ImgMetaData* metadata= img->GetMetaData();
	WorldCoor* wcs= metadata->GetWorldCoord();
	if(!wcs){
		ERROR_LOG("Cannot get WCS from image!");
		return -1;
	}	


	/*
	//## Set gaxis
	int status= SetWCSAxis(img,*xaxis,*yaxis);
	if(status<0){
		cerr<<"GraphicsUtils::UpdateGAxis(): WARN: Failed to update current gaxis!"<<endl;
		return -1;
	}
	*/

	//Get image ranges
	double xmin= gPad->GetUxmin();
	double xmax= gPad->GetUxmax();
	double ymin= gPad->GetUymin();
	double ymax= gPad->GetUymax();
	DEBUG_LOG("xmin/xmax="<<xmin<<"/"<<xmax<<" ymin/ymax="<<ymin<<"/"<<ymax);

	//Get range coord in WCS
	DEBUG_LOG("Find pixel2wcs coords crrespnding to new range...");
	double xmin_wcs, xmax_wcs, ymin_wcs, ymax_wcs;	
	AstroUtils::PixelToWCSCoords(img,wcs,xmin,ymin,xmin_wcs,ymin_wcs); 
	AstroUtils::PixelToWCSCoords(img,wcs,xmax,ymin,xmax_wcs,ymin_wcs);
	AstroUtils::PixelToWCSCoords(img,wcs,xmin,ymax,xmin_wcs,ymax_wcs);		
	double minx_wcs= min(xmin_wcs, xmax_wcs); 
	double maxx_wcs= max(xmin_wcs, xmax_wcs); 
	double miny_wcs= min(ymin_wcs, ymax_wcs); 
	double maxy_wcs= max(ymin_wcs, ymax_wcs);
	
	//Update axis range
	xaxis_wcs->SetX1(gPad->GetUxmin());
	xaxis_wcs->SetX2(gPad->GetUxmax());
	xaxis_wcs->SetY1(gPad->GetUymin());
	xaxis_wcs->SetY2(gPad->GetUymin());	
	//xaxis_wcs->SetY1(gPad->GetUymin()-40);
	//xaxis_wcs->SetY2(gPad->GetUymin()-40);	
	
	xaxis_wcs->SetWmin(minx_wcs);
	xaxis_wcs->SetWmax(maxx_wcs);

	yaxis_wcs->SetX1(gPad->GetUxmin());
	yaxis_wcs->SetX2(gPad->GetUxmin());
	//yaxis_wcs->SetX1(gPad->GetUxmin()-40);
	//yaxis_wcs->SetX2(gPad->GetUxmin()-40);
	
	yaxis_wcs->SetY1(gPad->GetUymin());
	yaxis_wcs->SetY2(gPad->GetUymax());
	yaxis_wcs->SetWmin(miny_wcs);
	yaxis_wcs->SetWmax(maxy_wcs);

	gPad->Update();
	
	return 0;

}//close UpdateGAxis()


int GraphicsUtils::SetThermalPalette(int ncolors){
	
	//Define palette colors	
	const int nc= 433;	
	double r[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00392157,0.00392157,0.00784314,0.00784314,0.0117647,0.0156863,0.0156863,0.0196078,0.0235294,0.027451,0.0313725,0.0352941,0.0392157,0.0431373,0.0470588,0.0509804,0.0509804,0.054902,0.0627451,0.0705882,0.0745098,0.0823529,0.0901961,0.0980392,0.105882,0.109804,0.117647,0.12549,0.133333,0.141176,0.14902,0.156863,0.164706,0.172549,0.180392,0.188235,0.196078,0.203922,0.211765,0.219608,0.223529,0.231373,0.235294,0.243137,0.247059,0.254902,0.258824,0.266667,0.270588,0.278431,0.286275,0.290196,0.298039,0.305882,0.309804,0.317647,0.321569,0.329412,0.337255,0.345098,0.352941,0.360784,0.364706,0.372549,0.380392,0.388235,0.392157,0.4,0.407843,0.415686,0.423529,0.427451,0.435294,0.439216,0.443137,0.45098,0.458824,0.466667,0.470588,0.478431,0.486275,0.494118,0.498039,0.505882,0.513725,0.517647,0.52549,0.529412,0.537255,0.541176,0.545098,0.552941,0.560784,0.568627,0.576471,0.584314,0.588235,0.596078,0.6,0.607843,0.611765,0.615686,0.623529,0.627451,0.635294,0.639216,0.643137,0.65098,0.654902,0.658824,0.662745,0.666667,0.670588,0.678431,0.682353,0.686275,0.690196,0.690196,0.694118,0.698039,0.701961,0.705882,0.709804,0.713725,0.717647,0.721569,0.72549,0.729412,0.729412,0.733333,0.737255,0.741176,0.745098,0.74902,0.74902,0.752941,0.752941,0.756863,0.756863,0.760784,0.764706,0.764706,0.768627,0.772549,0.776471,0.776471,0.780392,0.784314,0.788235,0.792157,0.792157,0.796078,0.796078,0.8,0.803922,0.807843,0.807843,0.811765,0.811765,0.815686,0.819608,0.819608,0.823529,0.823529,0.827451,0.827451,0.831373,0.831373,0.835294,0.835294,0.839216,0.843137,0.847059,0.847059,0.85098,0.854902,0.854902,0.858824,0.858824,0.862745,0.866667,0.866667,0.870588,0.870588,0.87451,0.87451,0.87451,0.878431,0.878431,0.878431,0.882353,0.886275,0.886275,0.890196,0.890196,0.894118,0.894118,0.894118,0.898039,0.898039,0.898039,0.901961,0.905882,0.905882,0.905882,0.909804,0.909804,0.909804,0.913725,0.913725,0.917647,0.917647,0.921569,0.921569,0.921569,0.921569,0.92549,0.92549,0.92549,0.92549,0.929412,0.929412,0.929412,0.933333,0.933333,0.933333,0.933333,0.937255,0.937255,0.937255,0.937255,0.941176,0.941176,0.941176,0.945098,0.945098,0.945098,0.945098,0.945098,0.945098,0.945098,0.945098,0.94902,0.94902,0.94902,0.952941,0.952941,0.952941,0.952941,0.956863,0.956863,0.956863,0.956863,0.956863,0.956863,0.956863,0.960784,0.960784,0.960784,0.960784,0.964706,0.964706,0.964706,0.968627,0.968627,0.968627,0.968627,0.972549,0.972549,0.972549,0.972549,0.972549,0.972549,0.972549,0.976471,0.976471,0.976471,0.976471,0.976471,0.976471,0.976471,0.976471,0.980392,0.980392,0.980392,0.984314,0.984314,0.984314,0.984314,0.988235,0.988235,0.988235,0.988235,0.992157,0.992157,0.992157,0.992157,0.992157,0.992157,0.992157,0.992157,0.992157,0.992157,0.992157,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,0.996078,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	double g[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00392157,0.00392157,0.00392157,0.00392157,0.00392157,0.00392157,0.00392157,0.00784314,0.00784314,0.00784314,0.0117647,0.0117647,0.0156863,0.0156863,0.0156863,0.0196078,0.0196078,0.0196078,0.0235294,0.0235294,0.0235294,0.027451,0.027451,0.0313725,0.0352941,0.0392157,0.0392157,0.0431373,0.0470588,0.0470588,0.0509804,0.054902,0.0588235,0.0627451,0.0666667,0.0705882,0.0745098,0.0745098,0.0784314,0.0823529,0.0862745,0.0901961,0.0941176,0.0941176,0.0980392,0.101961,0.105882,0.109804,0.109804,0.113725,0.117647,0.12549,0.129412,0.133333,0.137255,0.141176,0.145098,0.14902,0.152941,0.156863,0.164706,0.168627,0.172549,0.180392,0.184314,0.184314,0.188235,0.192157,0.196078,0.2,0.203922,0.207843,0.211765,0.215686,0.219608,0.223529,0.227451,0.231373,0.235294,0.239216,0.243137,0.247059,0.254902,0.258824,0.262745,0.266667,0.270588,0.27451,0.278431,0.282353,0.286275,0.290196,0.298039,0.298039,0.301961,0.301961,0.305882,0.309804,0.313725,0.317647,0.321569,0.32549,0.329412,0.337255,0.341176,0.345098,0.34902,0.352941,0.356863,0.360784,0.360784,0.364706,0.368627,0.372549,0.376471,0.380392,0.384314,0.388235,0.392157,0.396078,0.4,0.4,0.403922,0.407843,0.411765,0.415686,0.419608,0.419608,0.423529,0.427451,0.431373,0.435294,0.439216,0.443137,0.447059,0.45098,0.454902,0.458824,0.462745,0.466667,0.470588,0.478431,0.482353,0.486275,0.494118,0.498039,0.501961,0.505882,0.509804,0.513725,0.517647,0.521569,0.52549,0.529412,0.533333,0.533333,0.537255,0.541176,0.545098,0.54902,0.552941,0.552941,0.556863,0.560784,0.564706,0.568627,0.572549,0.576471,0.580392,0.584314,0.588235,0.596078,0.6,0.603922,0.611765,0.615686,0.623529,0.627451,0.631373,0.635294,0.639216,0.643137,0.65098,0.654902,0.658824,0.666667,0.670588,0.67451,0.678431,0.682353,0.686275,0.690196,0.694118,0.698039,0.701961,0.705882,0.709804,0.713725,0.721569,0.72549,0.72549,0.729412,0.733333,0.737255,0.741176,0.745098,0.752941,0.756863,0.760784,0.764706,0.768627,0.772549,0.776471,0.780392,0.784314,0.788235,0.792157,0.792157,0.796078,0.8,0.803922,0.807843,0.811765,0.811765,0.815686,0.819608,0.827451,0.831373,0.835294,0.839216,0.843137,0.847059,0.85098,0.854902,0.854902,0.858824,0.862745,0.862745,0.866667,0.870588,0.870588,0.87451,0.878431,0.882353,0.886275,0.886275,0.890196,0.894118,0.894118,0.898039,0.901961,0.901961,0.905882,0.909804,0.913725,0.917647,0.921569,0.921569,0.92549,0.929412,0.933333,0.933333,0.933333,0.937255,0.937255,0.941176,0.941176,0.945098,0.945098,0.945098,0.94902,0.94902,0.94902,0.952941,0.956863,0.956863,0.956863,0.960784,0.960784,0.960784,0.964706,0.964706,0.968627,0.968627,0.972549,0.972549,0.972549,0.972549,0.976471,0.976471,0.976471,0.980392,0.980392,0.984314,0.988235,0.988235,0.992157,0.992157,0.992157,0.996078,0.996078,0.996078,0.996078,1};
	double b[]={0.0392157,0.0784314,0.117647,0.145098,0.164706,0.180392,0.196078,0.211765,0.227451,0.243137,0.258824,0.27451,0.290196,0.309804,0.321569,0.333333,0.341176,0.34902,0.360784,0.368627,0.380392,0.388235,0.396078,0.403922,0.411765,0.419608,0.431373,0.439216,0.45098,0.454902,0.458824,0.462745,0.466667,0.470588,0.47451,0.482353,0.486275,0.490196,0.494118,0.501961,0.505882,0.513725,0.517647,0.521569,0.52549,0.529412,0.537255,0.537255,0.541176,0.545098,0.54902,0.552941,0.556863,0.556863,0.560784,0.564706,0.568627,0.572549,0.576471,0.576471,0.580392,0.584314,0.584314,0.588235,0.588235,0.588235,0.588235,0.592157,0.592157,0.592157,0.592157,0.596078,0.596078,0.596078,0.6,0.6,0.6,0.603922,0.603922,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.611765,0.611765,0.611765,0.611765,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.615686,0.611765,0.611765,0.611765,0.611765,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.607843,0.603922,0.603922,0.603922,0.6,0.6,0.6,0.6,0.596078,0.596078,0.596078,0.596078,0.592157,0.592157,0.588235,0.588235,0.584314,0.584314,0.584314,0.584314,0.584314,0.584314,0.580392,0.576471,0.576471,0.576471,0.572549,0.572549,0.572549,0.568627,0.568627,0.564706,0.564706,0.560784,0.556863,0.556863,0.552941,0.54902,0.545098,0.541176,0.537255,0.533333,0.529412,0.52549,0.521569,0.521569,0.517647,0.509804,0.505882,0.501961,0.494118,0.486275,0.482353,0.47451,0.470588,0.462745,0.458824,0.454902,0.447059,0.443137,0.435294,0.431373,0.419608,0.411765,0.403922,0.396078,0.392157,0.384314,0.376471,0.368627,0.360784,0.352941,0.341176,0.329412,0.317647,0.305882,0.290196,0.278431,0.266667,0.254902,0.239216,0.227451,0.215686,0.2,0.188235,0.176471,0.164706,0.14902,0.137255,0.12549,0.113725,0.109804,0.105882,0.0980392,0.0941176,0.0862745,0.0823529,0.0784314,0.0745098,0.0705882,0.0627451,0.0588235,0.054902,0.0509804,0.0470588,0.0470588,0.0431373,0.0392157,0.0392157,0.0352941,0.0352941,0.0313725,0.0313725,0.0313725,0.027451,0.027451,0.0235294,0.0235294,0.0196078,0.0196078,0.0196078,0.0156863,0.0156863,0.0156863,0.0156863,0.0117647,0.0117647,0.0117647,0.0117647,0.0117647,0.0117647,0.0117647,0.00784314,0.00784314,0.00784314,0.00784314,0.00392157,0.00392157,0.00392157,0.00392157,0.00392157,0.00392157,0.00392157,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00392157,0.00392157,0.00392157,0.00392157,0.00784314,0.00784314,0.0117647,0.0156863,0.0156863,0.0196078,0.0235294,0.0313725,0.0352941,0.0392157,0.0392157,0.0431373,0.0470588,0.0509804,0.054902,0.054902,0.0627451,0.0705882,0.0784314,0.0862745,0.0980392,0.105882,0.117647,0.12549,0.133333,0.141176,0.14902,0.156863,0.168627,0.180392,0.192157,0.207843,0.219608,0.235294,0.247059,0.262745,0.27451,0.286275,0.301961,0.313725,0.329412,0.341176,0.356863,0.372549,0.388235,0.403922,0.415686,0.431373,0.447059,0.466667,0.482353,0.501961,0.521569,0.541176,0.556863,0.572549,0.588235,0.603922,0.619608,0.635294,0.65098,0.666667,0.686275,0.701961,0.713725,0.729412,0.741176,0.756863,0.768627,0.780392,0.792157,0.803922,0.819608,0.831373,0.847059,0.858824,0.87451,0.886275,0.898039,0.909804,0.921569,0.933333,0.945098,0.956863,0.964706};
	double stop[]={0,0.00231481,0.00462963,0.00694444,0.00925926,0.0115741,0.0138889,0.0162037,0.0185185,0.0208333,0.0231481,0.025463,0.0277778,0.0300926,0.0324074,0.0347222,0.037037,0.0393519,0.0416667,0.0439815,0.0462963,0.0486111,0.0509259,0.0532407,0.0555556,0.0578704,0.0601852,0.0625,0.0648148,0.0671296,0.0694444,0.0717593,0.0740741,0.0763889,0.0787037,0.0810185,0.0833333,0.0856481,0.087963,0.0902778,0.0925926,0.0949074,0.0972222,0.099537,0.101852,0.104167,0.106481,0.108796,0.111111,0.113426,0.115741,0.118056,0.12037,0.122685,0.125,0.127315,0.12963,0.131944,0.134259,0.136574,0.138889,0.141204,0.143519,0.145833,0.148148,0.150463,0.152778,0.155093,0.157407,0.159722,0.162037,0.164352,0.166667,0.168981,0.171296,0.173611,0.175926,0.178241,0.180556,0.18287,0.185185,0.1875,0.189815,0.19213,0.194444,0.196759,0.199074,0.201389,0.203704,0.206019,0.208333,0.210648,0.212963,0.215278,0.217593,0.219907,0.222222,0.224537,0.226852,0.229167,0.231481,0.233796,0.236111,0.238426,0.240741,0.243056,0.24537,0.247685,0.25,0.252315,0.25463,0.256944,0.259259,0.261574,0.263889,0.266204,0.268519,0.270833,0.273148,0.275463,0.277778,0.280093,0.282407,0.284722,0.287037,0.289352,0.291667,0.293981,0.296296,0.298611,0.300926,0.303241,0.305556,0.30787,0.310185,0.3125,0.314815,0.31713,0.319444,0.321759,0.324074,0.326389,0.328704,0.331019,0.333333,0.335648,0.337963,0.340278,0.342593,0.344907,0.347222,0.349537,0.351852,0.354167,0.356481,0.358796,0.361111,0.363426,0.365741,0.368056,0.37037,0.372685,0.375,0.377315,0.37963,0.381944,0.384259,0.386574,0.388889,0.391204,0.393519,0.395833,0.398148,0.400463,0.402778,0.405093,0.407407,0.409722,0.412037,0.414352,0.416667,0.418981,0.421296,0.423611,0.425926,0.428241,0.430556,0.43287,0.435185,0.4375,0.439815,0.44213,0.444444,0.446759,0.449074,0.451389,0.453704,0.456019,0.458333,0.460648,0.462963,0.465278,0.467593,0.469907,0.472222,0.474537,0.476852,0.479167,0.481481,0.483796,0.486111,0.488426,0.490741,0.493056,0.49537,0.497685,0.5,0.502315,0.50463,0.506944,0.509259,0.511574,0.513889,0.516204,0.518519,0.520833,0.523148,0.525463,0.527778,0.530093,0.532407,0.534722,0.537037,0.539352,0.541667,0.543981,0.546296,0.548611,0.550926,0.553241,0.555556,0.55787,0.560185,0.5625,0.564815,0.56713,0.569444,0.571759,0.574074,0.576389,0.578704,0.581019,0.583333,0.585648,0.587963,0.590278,0.592593,0.594907,0.597222,0.599537,0.601852,0.604167,0.606481,0.608796,0.611111,0.613426,0.615741,0.618056,0.62037,0.622685,0.625,0.627315,0.62963,0.631944,0.634259,0.636574,0.638889,0.641204,0.643519,0.645833,0.648148,0.650463,0.652778,0.655093,0.657407,0.659722,0.662037,0.664352,0.666667,0.668981,0.671296,0.673611,0.675926,0.678241,0.680556,0.68287,0.685185,0.6875,0.689815,0.69213,0.694444,0.696759,0.699074,0.701389,0.703704,0.706019,0.708333,0.710648,0.712963,0.715278,0.717593,0.719907,0.722222,0.724537,0.726852,0.729167,0.731481,0.733796,0.736111,0.738426,0.740741,0.743056,0.74537,0.747685,0.75,0.752315,0.75463,0.756944,0.759259,0.761574,0.763889,0.766204,0.768519,0.770833,0.773148,0.775463,0.777778,0.780093,0.782407,0.784722,0.787037,0.789352,0.791667,0.793981,0.796296,0.798611,0.800926,0.803241,0.805556,0.80787,0.810185,0.8125,0.814815,0.81713,0.819444,0.821759,0.824074,0.826389,0.828704,0.831019,0.833333,0.835648,0.837963,0.840278,0.842593,0.844907,0.847222,0.849537,0.851852,0.854167,0.856481,0.858796,0.861111,0.863426,0.865741,0.868056,0.87037,0.872685,0.875,0.877315,0.87963,0.881944,0.884259,0.886574,0.888889,0.891204,0.893519,0.895833,0.898148,0.900463,0.902778,0.905093,0.907407,0.909722,0.912037,0.914352,0.916667,0.918981,0.921296,0.923611,0.925926,0.928241,0.930556,0.93287,0.935185,0.9375,0.939815,0.94213,0.944444,0.946759,0.949074,0.951389,0.953704,0.956019,0.958333,0.960648,0.962963,0.965278,0.967593,0.969907,0.972222,0.974537,0.976852,0.979167,0.981481,0.983796,0.986111,0.988426,0.990741,0.993056,0.99537,0.997685,1};
		
	int FI = TColor::CreateGradientColorTable(nc,stop,r,g,b,ncolors);
	int ThermalPalette[ncolors];
	for (int i=0;i<ncolors;i++) {
		int colorId= FI+i;
		ThermalPalette[i] = colorId;
	}

	gStyle->SetNumberContours(ncolors);
	gStyle->SetPalette(ncolors,ThermalPalette);

	return 0;

}//close SetThermalPalette()

int GraphicsUtils::SetHotColdPalette(int ncolors){

	//Define palette colors	
	const int nc = 12;
	double r[]    = {0.164, 0.150, 0.250, 0.450, 0.670, 0.880, 1.000, 1.000,1.000,0.970,0.850,0.650};
  double g[]    = {0.043, 0.306, 0.630, 0.853, 0.973, 1.000,1.000,0.880,0.679,0.430,0.150,0.000};
  double b[]    = {0.850,1.000,1.000,1.000,1.000,1.000,0.750,0.600,0.450,0.370,0.196,0.130};
  
  double r_HotToCold[nc];
  double g_HotToCold[nc];
  double b_HotToCold[nc];
  double stop_HotToCold[] = {0.,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
	//double stop_HotToCold[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.6,0.7,0.8,1};

  for (int i=0;i<nc;i++) {
   r_HotToCold[i]= r[nc-i-1];
   g_HotToCold[i]= g[nc-i-1];
   b_HotToCold[i]= b[nc-i-1];
  }

	int FI = TColor::CreateGradientColorTable(nc, stop_HotToCold, r_HotToCold, g_HotToCold, b_HotToCold, ncolors);
	int HotToColdPalette[ncolors];

  for (int i=0;i<ncolors;i++) {
  	HotToColdPalette[i] = FI + i;
  }
	
	gStyle->SetNumberContours(ncolors);
	gStyle->SetPalette(ncolors,HotToColdPalette);


	return 0;

}//close SetHotColdPalette()


int GraphicsUtils::SetColdHotPalette(int ncolors){

	//Define palette colors	
	const int nc = 12;
	double r[]    = {0.164, 0.150, 0.250, 0.450, 0.670, 0.880, 1.000, 1.000,1.000,0.970,0.850,0.650};
  double g[]    = {0.043, 0.306, 0.630, 0.853, 0.973, 1.000,1.000,0.880,0.679,0.430,0.150,0.000};
  double b[]    = {0.850,1.000,1.000,1.000,1.000,1.000,0.750,0.600,0.450,0.370,0.196,0.130};
  double stop[] = {0.,0.2,0.3,0.4,0.5,0.7,0.75,0.8,0.85,0.9,0.95,1};

	int ColdToHotPalette[ncolors];
  int FI = TColor::CreateGradientColorTable(nc,stop,r,g,b,ncolors);
  for (int i=0;i<ncolors;i++) {
		ColdToHotPalette[i] = FI+i;
	}

	gStyle->SetNumberContours(ncolors);
	gStyle->SetPalette(ncolors,ColdToHotPalette);

	return 0;

}//close SetColdHotPalette()


int GraphicsUtils::SetBWPalette(int ncolors){

	//Define palette colors	
	const int nc = 2;
  double r[]   = { 0.00, 1.00};
  double g[] = { 0.00, 1.00};
  double b[]  = { 0.00, 1.00};
  double stop[] = { 0.00, 1.00};
                                                                                
  int GrayPalette[ncolors];
  int FI= TColor::CreateGradientColorTable(nc,stop,r,g,b,ncolors);	
	for (int i=0;i<ncolors;i++) {
  	GrayPalette[i] = FI + i;
  }

	gStyle->SetNumberContours(ncolors);
	gStyle->SetPalette(ncolors,GrayPalette);

	return 0;

}//close SetBWPalette()


}//close namespace
