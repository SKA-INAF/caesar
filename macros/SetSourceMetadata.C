#include <Source.h>
#include <SourceExporter.h>
#include <iostream>

using namespace std;
using namespace Caesar;

void SetSourceMetadata(std::string filename,std::string imgfile,std::string outfile="output.root")
{
	//## Read image file
	Image* img= new Image;	
	if(img->ReadFITS(imgfile)<0){
		cerr<<"ERROR: Cannot read FITS image from file "<<imgfile<<"!"<<endl;
		exit(1);
	}
	ImgMetaData* imgmetadata= img->GetMetaData();
	if(!imgmetadata){
		cerr<<"ERROR: Cannot get metadata from image "<<imgfile<<"!"<<endl;
		exit(1);
	}
	double Freq= imgmetadata->Freq;
	std::string FreqUnit= imgmetadata->FreqUnit;
	double dFreq= imgmetadata->dFreq;
	double FreqRef= imgmetadata->FreqRef;

	//## Read input file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()) {
		cerr<<"ERROR: Cannot open file "<<filename<<"!"<<endl;
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
	Source* source= 0;
	std::vector<Source*> sources;
	
	TFile* outputFile= new TFile(outfile.c_str(),"RECREATE");
	TTree* outputTree= new TTree("SourceInfo","SourceInfo");
	outputTree->Branch("Source",&source);

	for(int i=0;i<SourceInfo->GetEntries();i++){
		SourceInfo->GetEntry(i);


		//Update metadata spectral info
		/*
		ImgMetaData* metadata= aSource->GetImageMetaData();
		cout<<"BEFORE: Source "<<aSource->GetName()<<": Freq="<<metadata->Freq<<", unit="<<metadata->FreqUnit<<endl;
		metadata->Freq= Freq;
		metadata->FreqUnit= FreqUnit;
		metadata->dFreq= dFreq;
		metadata->FreqRef= FreqRef;
		cout<<"AFTER: Source "<<aSource->GetName()<<": Freq="<<aSource->GetImageMetaData()->Freq<<", unit="<<aSource->GetImageMetaData()->FreqUnit<<endl;
		*/

		//copy source
		source= new Source;
		*source= *aSource;
	
		//source->SetImageMetaData(metadata);
		
		//Update metadata spectral info
		ImgMetaData* metadata= source->GetImageMetaData();
		cout<<"BEFORE: Source "<<source->GetName()<<": Freq="<<metadata->Freq<<", unit="<<metadata->FreqUnit<<endl;
		metadata->Freq= Freq;
		metadata->FreqUnit= FreqUnit;
		metadata->dFreq= dFreq;
		metadata->FreqRef= FreqRef;
		cout<<"AFTER: Source "<<source->GetName()<<": Freq="<<source->GetImageMetaData()->Freq<<", unit="<<source->GetImageMetaData()->FreqUnit<<endl;
		

		//Fill tree
		outputTree->Fill();

		//Add to collection
		sources.push_back(source);

	}//end loop sources

	//Close file
	if(outputFile){
		outputFile->cd();
		outputTree->Write();
		outputFile->Close();
	}

}//close macro
