
#include <Source.h>
#include <SourceExporter.h>
#include <iostream>

using namespace std;
using namespace Caesar;


void WriteSourceDataTree(std::string filename="",std::string outfilename="sourceDataTree.root",std::string outfilename_comp="sourceCompDataTree.root")
{
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
	std::vector<Source*> sources;
	
	for(int i=0;i<SourceInfo->GetEntries();i++){
		SourceInfo->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);

	}//end loop sources

	//## Write source collection to simple TTree
	bool dumpNestedSourceInfo= true;
	int wcsType= eJ2000;
	WCS* wcs= 0;
	bool writeAdditionalSourceInfo= true;

	cout<<"INFO: Writing source data TTree..."<<endl;
	SourceExporter::WriteToROOT(outfilename,sources,dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo);

	cout<<"INFO: Writing source component data TTree..."<<endl;
	SourceExporter::WriteComponentsToROOT(outfilename_comp,sources,dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo);

}//close macro


