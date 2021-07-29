#include <Image.h>
#include <Source.h>

#include <TTree.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Caesar;


int ExtractSourceCutout(std::string imgfilename,std::string catalogfile,int cutoutSize=100,int offset=0,bool randomizeOffset=false, std::string sourcefilelist="",std::string label="source",bool writeHeader=false,std::string fileprefix="",bool maskAllSources=false)
{
	//Require cutout size to be "even"
	if(cutoutSize%2!=0){
		cerr<<"ERROR: Require an even cutout size!"<<endl;
		exit(1);
	}
	if(offset>=cutoutSize){
		cerr<<"ERROR: Offset larger than cutout size!"<<endl;
		exit(1);
	}

	//Get current directory
	std::string currentDir= gSystem->pwd();
	
	delete gRandom;	
	gRandom= new TRandom3(0);

	//=================================
	//==      READ IMAGE 
	//=================================
	//- Read image filename
	Image* img= new Image;
	if(img->ReadFITS(imgfilename)<0){
		cerr<<"ERROR: Failed to open image file "<<imgfilename<<"!"<<endl;
		return -1;
	}

	//=================================
	//==      READ SOURCES 
	//=================================
	//- Read source catalog file
	TFile* inputFile= new TFile(catalogfile.c_str(),"READ");
	if(!inputFile){
		cerr<<"ERROR: Cannot get sources from file "<<catalogfile<<"!"<<endl;
		exit(1);
	}

	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");	
	if(!SourceInfo){
		cerr<<"ERROR: Cannot get sources from file "<<catalogfile<<"!"<<endl;
		exit(1);
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	cout<<"# "<<SourceInfo->GetEntries()<<" sources to be drawn..."<<endl;
	
	//- Store source list 
	std::vector<Source*> sources;
	
	for(int i=0;i<SourceInfo->GetEntries();i++){
		SourceInfo->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	//=================================
	//==      READ SOURCE NAME LIST 
	//=================================
	TTree* SourceNameInfo= 0;
	bool selectSources= false;
	std::vector<std::string> selectedSourceNames;
	std::vector<Source*> sources_sel;

	if(sourcefilelist!=""){
		selectSources= true;

		SourceNameInfo= new TTree;
		SourceNameInfo->ReadFile(sourcefilelist.c_str(),"sname/C");

		char sname[100];
		SourceNameInfo->SetBranchAddress("sname",sname);

		for(int i=0;i<SourceNameInfo->GetEntries();i++){
			SourceNameInfo->GetEntry(i);
			std::string sourceName= std::string(sname);
			selectedSourceNames.push_back(sourceName);
		}

		//Select sources
		for(size_t i=0;i<sources.size();i++)
		{
			Source* source= sources[i];
			std::string sname= source->GetName();
			std::vector<std::string>::iterator it= std::find(selectedSourceNames.begin(), selectedSourceNames.end(),sname);
			bool found= (it!=selectedSourceNames.end());	
			if(found) sources_sel.push_back(source);
		}

		cout<<"# "<<sources_sel.size()<<" sources selected for cutout ..."<<endl;

	}//close if
	else{
		sources_sel= sources;
	}

	//=================================
	//==      COMPUTE SOURCE MASK
	//=================================
	cout<<"INFO: Computing selected source mask..."<<endl;
	bool copyMetaData= true;
	bool resetStats= true;
	Image* smask= 0;

	if(maskAllSources){
		smask= img->GetCloned("smask",copyMetaData,resetStats);	
		smask->Reset();

		for(size_t i=0;i<sources_sel.size();i++){
			Source* source= sources_sel[i];
			for(int l=0;l<source->GetNPixels();l++){
				Pixel* pixel= source->GetPixel(l);
				long int id= pixel->id;
				double x= pixel->x;
				double y= pixel->y;
				long int gBinId= img->FindBin(x,y);
				if(gBinId<0){
					cerr<<"WARN: Cannot find gbin of pixel ("<<x<<","<<y<<"), skip..."<<endl;
					continue;
				}		
				smask->SetPixelValue(gBinId,1);

			}//end loop pixels		
		}//end loop sel sources
	}//close if

	/*
	bool isBinary= true;	
	bool invertMask= false;
	Image* smask= img->GetSourceMask(sources_sel,isBinary,invertMask);
	if(!smask){
		cerr<<"ERROR: Failed to compute source mask!"<<endl;
		exit(1);
	} 
	*/

	//=================================
	//==      MAKE CUTOUTS 
	//=================================
	cout<<"INFO: Making cutouts..."<<endl;
	FILE* fout= fopen("sources_bb.dat","a");//append mode
	if(writeHeader) fprintf(fout,"# filepath,x1,y1,x2,y2,className");

	FILE* fout_mask= fopen("sources_mask.dat","a");//append mode
	if(writeHeader) fprintf(fout_mask,"# filepath,filepath_mask,className");

	int cutoutHalfWidth= (int)(cutoutSize/2.);

	for(size_t i=0;i<sources_sel.size();i++)
	{
		Source* source= sources_sel[i];
		std::string sname= source->GetName();
		
		/*
		//Select source
		if(selectSources){
			std::vector<std::string>::iterator it= std::find(selectedSourceNames.begin(), selectedSourceNames.end(),sname);
			bool found= (it!=selectedSourceNames.end());	
			if(!found) continue;
		}
		*/

		
		double x0_s= source->X0;
		double y0_s= source->Y0;		
		long int gbin= img->FindBin(x0_s,y0_s);
		long int x0= img->GetBinX(gbin);
		long int y0= img->GetBinY(gbin);
		
		//Apply image offset
		//x0+= minX;
		//y0+= minY;

		//Apply cutout offset
		int offset_x= offset;
		int offset_y= offset;
		if(randomizeOffset){	
			double rand_x= gRandom->Uniform(-offset,offset);
			double rand_y= gRandom->Uniform(-offset,offset);
			offset_x= (int)rand_x;
			offset_y= (int)rand_y;
			x0+= offset_x;
			y0+= offset_y;
		}
		else{
			x0+= offset_x;
			y0+= offset_y;
		}

		//Extract image cutout
		std::string cutoutImgName= "";
		std::string cutoutImgName_mask= "";
		if(fileprefix=="") {
			cutoutImgName= Form("cutout_%s_dx%d_dy%d.fits",sname.c_str(),offset_x,offset_y);	
			cutoutImgName_mask= Form("mask_%s_dx%d_dy%d.fits",sname.c_str(),offset_x,offset_y);
			if(offset==0 && !randomizeOffset){
				cutoutImgName= Form("cutout_%s.fits",sname.c_str());
				cutoutImgName_mask= Form("mask_%s.fits",sname.c_str());
			}
		}
		else {
			cutoutImgName= Form("cutout_%s_%s_dx%d_dy%d.fits",fileprefix.c_str(),sname.c_str(),offset_x,offset_y);
			cutoutImgName_mask= Form("mask_%s_%s_dx%d_dy%d.fits",fileprefix.c_str(),sname.c_str(),offset_x,offset_y);
			if(offset==0 && !randomizeOffset){
				cutoutImgName= Form("cutout_%s_%s.fits",fileprefix.c_str(),sname.c_str());
				cutoutImgName_mask= Form("mask_%s_%s.fits",fileprefix.c_str(),sname.c_str());
			}
		}

		long int xmin= x0 - cutoutHalfWidth;
		long int xmax= x0 + cutoutHalfWidth-1;
		//long int xmax= x0 + cutoutHalfWidth;
		long int ymin= y0 - cutoutHalfWidth;
		long int ymax= y0 + cutoutHalfWidth-1;
		//long int ymax= y0 + cutoutHalfWidth;
		
		cout<<"(x0,y0)=("<<x0<<","<<y0<<"), xmin/max="<<xmin<<"/"<<xmax<<", ymin/ymax="<<ymin<<"/"<<ymax<<endl;
	
		Image* cutoutImg= img->GetTile(xmin,xmax,ymin,ymax);


		//Find mask
		Image* cutoutImg_mask= 0;
		if(maskAllSources){
			cutoutImg_mask= smask->GetTile(xmin,xmax,ymin,ymax);
		}
		else{		
			cutoutImg_mask= img->GetTile(xmin,xmax,ymin,ymax);
			cutoutImg_mask->Reset();

			for(int l=0;l<source->GetNPixels();l++){
				Pixel* pixel= source->GetPixel(l);
				long int id= pixel->id;
				double x= pixel->x - xmin;
				double y= pixel->y - ymin;
				//long int gBinId= cutoutImg_mask->FindBin(x,y);
				long int gBinId= cutoutImg_mask->GetBin(x,y);
				if(gBinId<0){
					cerr<<"WARN: Cannot find gbin of pixel ("<<x<<","<<y<<"), skip..."<<endl;
					continue;
				}		
				cutoutImg_mask->SetPixelValue(gBinId,1);

			}//end loop pixels
		}//close if

		if(!cutoutImg || !cutoutImg_mask){
			cerr<<"ERROR: Failed to extract cutout for source "<<sname<<", skip to next!"<<endl;
			continue;
		}

		//Get source bounding box
		float sourceXmin= 0;
		float sourceXmax= 0;
		float sourceYmin= 0;
		float sourceYmax= 0;
		source->GetSourceRange(sourceXmin,sourceXmax,sourceYmin,sourceYmax);
	
		long int gbin_P1= img->FindBin(sourceXmin,sourceYmin);
		long int gbin_P2= img->FindBin(sourceXmax,sourceYmax);
		long int x1= img->GetBinX(gbin_P1);
		long int y1= img->GetBinY(gbin_P1);
		long int x2= img->GetBinX(gbin_P2);
		long int y2= img->GetBinY(gbin_P2);

		cout<<"INFO: Extracting cutout for source "<<sname<<" (file="<<cutoutImgName<<"), Pos("<<x0_s<<","<<y0_s<<"), PixPos("<<x0<<","<<y0<<"), Source BBox("<<sourceXmin<<","<<sourceXmax<<","<<sourceYmin<<","<<sourceYmax<<"), Source Pix BBox("<<x1<<","<<x2<<","<<y1<<","<<y2<<")"<<", Cutout xrange=["<<xmin<<","<<xmax<<"], Cutout yrange=["<<ymin<<","<<ymax<<"]"<<endl;

		/*
		double x1= sourceXmin - xmin;
		double x2= sourceXmax - xmin;
		double y1= sourceYmin - ymin;
	 	double y2= sourceYmax - ymin;	

		int sourceBB_x1= (int)(x1);
		int sourceBB_x2= (int)(x2);
		int sourceBB_y1= (int)(y1);
	 	int sourceBB_y2= (int)(y2);
		*/
		
		//ML algorithms expecting (x1,y1)=TOP LEFT  and (x2,y2) BOTTOM RIGHT
		//NB: Here the y-axis origin is inverted so we take (x1,y1)=BOTTOM LEFT and (x2,y2)=TOP RIGHT
		int sourceBB_tl_x= x1 - xmin;
		int sourceBB_tl_y= y2 - ymin;

		int sourceBB_br_x= x2 - xmin;		
	 	int sourceBB_br_y= y1 - ymin;

		int sourceBB_bl_x= x1 - xmin;		
	 	int sourceBB_bl_y= y1 - ymin;

		int sourceBB_tr_x= x2 - xmin;
		int sourceBB_tr_y= y2 - ymin;

		//Write cutout to file
		cutoutImg->WriteFITS(cutoutImgName);
		cutoutImg_mask->WriteFITS(cutoutImgName_mask);
				

		//Write bounding box info
		std::string cutoutImgName_fullpath= currentDir + std::string("/") + cutoutImgName;
		std::string cutoutImgName_mask_fullpath= currentDir + std::string("/") + cutoutImgName_mask;
		//fprintf(fout,"%s,%d,%d,%d,%d,%s\n",cutoutImgName_fullpath.c_str(),sourceBB_tl_x,sourceBB_tl_y,sourceBB_br_x,sourceBB_br_y,label.c_str());

		fprintf(fout,"%s,%d,%d,%d,%d,%s\n",cutoutImgName_fullpath.c_str(),sourceBB_bl_x,sourceBB_bl_y,sourceBB_tr_x,sourceBB_tr_y,label.c_str());
		fprintf(fout_mask,"%s,%s,%s\n",cutoutImgName_fullpath.c_str(),cutoutImgName_mask_fullpath.c_str(),label.c_str());

		std::stringstream ss;
		ss<<"polygon ";
		ss<<sourceBB_bl_x+1<<" "<<sourceBB_bl_y+1<<" ";//bottom-left
		ss<<sourceBB_br_x+1<<" "<<sourceBB_br_y+1<<" ";//bottom-right
		ss<<sourceBB_tr_x+1<<" "<<sourceBB_tr_y+1<<" ";//top-right
		ss<<sourceBB_tl_x+1<<" "<<sourceBB_tl_y+1<<" ";//top-left
		ss<<"# text={"<<sname<<"}";

		
		//Write bounding box DS9 region
		std::string regionFileName= "";
		if(fileprefix=="") {
			regionFileName= Form("ds9_%s_dx%d_dy%d.reg",sname.c_str(),offset_x,offset_y);
			if(offset==0 && !randomizeOffset){
				regionFileName= Form("ds9_%s.reg",sname.c_str());
			}
		}
		else {
			regionFileName= Form("ds9_%s_%s_dx%d_dy%d.reg",fileprefix.c_str(),sname.c_str(),offset_x,offset_y);
			if(offset==0 && !randomizeOffset){
				regionFileName= Form("ds9_%s_%s.reg",fileprefix.c_str(),sname.c_str());
			}
		}
		
		FILE* fout_region= fopen(regionFileName.c_str(),"w");
		fprintf(fout_region,"global color=yellow font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
		fprintf(fout_region,"image\n");
		fprintf(fout_region,"%s\n",ss.str().c_str());
		fclose(fout_region);

	}//end loop sources

	fclose(fout);
	fclose(fout_mask);
	

	return 0;

}//close macro


