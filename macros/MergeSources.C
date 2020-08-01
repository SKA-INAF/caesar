int ReadMergeDataFile(std::vector<std::vector<std::string>>& mergedSourceNames,std::string fileName);

int MergeSources(std::string catalogFileName,std::string mergeSourceFileName,std::string outputFileName="output.root")
{
	//Create output TTree
	TFile* outputFile= new TFile(outputFileName.c_str(),"RECREATE");
	Source* aSource= 0;
	TTree* outputTree= new TTree("SourceInfo","SourceInfo");
	outputTree->Branch("Source",&aSource);


	//===================================
	//==    READ MERGE DATA FILE
	//===================================
	std::vector<std::vector<std::string>> mergedSourceNames;
	if(ReadMergeDataFile(mergedSourceNames,mergeSourceFileName)<0)
	{
		ERROR_LOG("Failed to read merge data file "<<mergeSourceFileName<<"!");
		return -1;		
	}

	//=======================================
	//==    READ SOURCE CATALOG DATA FILE
	//=======================================
	//- Read input file
	TFile* inputFile= new TFile(catalogFileName.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()) {
		ERROR_LOG("Cannot open file "<<catalogFileName<<"!");
		return -1;
	}
	
	//- Get source tree
	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");	
	if(!SourceInfo){
		ERROR_LOG("Cannot get sources from file "<<catalogFileName<<"!");
		return -1;
	}

	
	SourceInfo->SetBranchAddress("Source",&aSource);

	
	
	//- Store source list 
	std::map<std::string,Source*> sourceMap;
	
	for(int i=0;i<SourceInfo->GetEntries();i++)
	{
		SourceInfo->GetEntry(i);
		
		//Add mother source
		Source* source= new Source;
		*source= *aSource;
		
		std::string sname= source->GetName();
		sourceMap[sname]= source;

		//Add nested sources
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			Source* nestedSource= new Source;
			*nestedSource= *(nestedSources[j]);

			std::string sname_nested= nestedSource->GetName();			
			sourceMap[sname_nested]= nestedSource;

		}//end loop nested sources

	}//end loop sources

	INFO_LOG("# "<<sourceMap.size()<<" sources ...");
	
	//=======================================
	//==    MERGE SOURCE
	//=======================================
	for(size_t j=0;j<mergedSourceNames.size();j++){

		bool sourcesFound= true;
		std::vector<Source*> mergedSources;

		for(size_t k=0;k<mergedSourceNames[j].size();k++){
			std::string sname= mergedSourceNames[j][k];
				
			auto it = sourceMap.find(sname);
			if(it==sourceMap.end()){
				WARN_LOG("Cannot find source "<<sname<<"!");
				sourcesFound= false;
				break;
			}
			Source* mergedSource= sourceMap[sname];
	
			mergedSources.push_back(mergedSource);

		}//end loop sources to be merged

		if(!sourcesFound || mergedSources.empty()){
			WARN_LOG("One/more sources to be merged not found in map, skip this...");
			continue;
		}

		aSource= new Source;
		*aSource= *(mergedSources[0]);

		bool copyPixels= false;
		bool checkIfAdjacent= false;
		bool computeStatPars= false;
		bool computeMorphPars= false;
		bool sumMatchingPixels= false;
		std::stringstream ss;
		ss<<mergedSources[0]->GetName();
		for(size_t k=1;k<mergedSources.size();k++){
			aSource->MergeSource(mergedSources[k],copyPixels,checkIfAdjacent,computeStatPars,computeMorphPars,sumMatchingPixels);
			ss<<"_"<<mergedSources[k]->GetName();
		}

		//Set name of merged source
		aSource->SetName(ss.str());
		
		//Compute stats of merged source
		bool computeRobustStats= true;
		bool forceRecomputing= true;
			
		aSource->ComputeStats(computeRobustStats,forceRecomputing);

		//Fill TTree
		outputTree->Fill();

	}//end loop sources

	outputFile->cd();
	outputTree->Write();
	outputFile->Close();

	return 0;

}//close macro



int ReadMergeDataFile(std::vector<std::vector<std::string>>& mergedSourceNames,std::string fileName)
{
	//Clear data
	mergedSourceNames.clear();

	//Check filename
	std::ifstream fileStream(fileName);	
	if (fileStream.fail() || !fileStream.is_open()){
		ERROR_LOG("Failed to open file "<<fileName<<" for reading...");
		return -1;
	}

	//Read file line-by-line and parse
	std::string line;
	int lineCounter= 0;

	while (std::getline(fileStream, line)) 
	{  	
		
		std::istringstream iss(line);	
		lineCounter++;
		
		//Parse fields
		std::string field;
		std::vector<std::string> data;
		std::stringstream ss;
		ss<<"{";
		while (std::getline(iss, field, ' ')) 
		{
			if(field=="\n" || field=="") continue;
			
  		ss<<field<<",";
			data.push_back(field);			
		}	
		ss<<"}";
		cout<<ss.str()<<". data.size()="<<data.size()<<endl;
		mergedSourceNames.push_back(data);
		
	}//end file read
				
	
	return 0;

}//close macro
