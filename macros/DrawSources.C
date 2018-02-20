#include <Image.h>
#include <Source.h>

#include <iostream>

using namespace std;
using namespace Caesar;


void DrawSources(std::string filename=""){

	//## Read input file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()) {
		cerr<<"ERROR: Cannot open file "<<filename<<"!"<<endl;
		exit(1);
	}

	//## Get input image
	Image* img= (Image*)inputFile->Get("img");
	if(!img){
		cerr<<"ERROR: Cannot get image from file "<<filename<<"!"<<endl;
		exit(1);
	}
	
	//## Get source tree
	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");	
	if(!SourceInfo){
		cerr<<"ERROR: Cannot get sources from file "<<filename<<"!"<<endl;
		exit(1);
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	cout<<"# "<<SourceInfo->GetEntries()<<" sources to be drawn..."<<endl;
	
	//## Store source list 
	std::vector<Source*> sources;
	
	for(int i=0;i<SourceInfo->GetEntries();i++){
		SourceInfo->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	//## Draw plot
	bool useCurrentCanvas= false;
	bool drawFull= false;
	
	TCanvas* Plot= new TCanvas("Plot","Plot");
	Plot->cd();

	TH2D* img_histo= img->GetHisto2D("histo");
	img_histo->SetStats(0);
	img_histo->Draw("COLZ");

	for(size_t i=0;i<sources.size();i++){
		int source_type= sources[i]->Type;
		bool isAtEdge= sources[i]->IsAtEdge();
		double x= sources[i]->X0;
		double y= sources[i]->Y0;

		int color= kBlack;
		int lineStyle= kSolid;

		if(source_type==Source::eExtended) color= kGreen;
		else if(source_type==Source::ePointLike) color= kRed;
		else if(source_type==Source::eCompact) color= kBlue;	
		else color= kGray;
		if(isAtEdge) {
			color+= 2;
			lineStyle= kDashed;
		}
		
		cout<<"Source no. "<<i+1<<" (x,y)=("<<x<<","<<y<<"), type="<<source_type<<" (color="<<color<<")"<<endl;
	
		sources[i]->Draw(false,false,false,color,lineStyle);		
	}//end loop sources

}//close macro

