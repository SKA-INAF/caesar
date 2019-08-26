#include <Source.h>
#include <SourceExporter.h>
#include <iostream>

using namespace std;
using namespace Caesar;


void WriteSourceJsonCatalog(std::string filename="",std::string outfilename="sourceCatalog.json",std::string outfilename_comp="sourceCompCatalog.json")
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
	cout<<"INFO: Writing source json catalog ..."<<endl;
	SourceExporter::WriteToJson(outfilename,sources);

	//cout<<"INFO: Writing source component json catalog ..."<<endl;
	//SourceExporter::WriteComponentsToJson(outfilename_comp,sources);

}//close macro

