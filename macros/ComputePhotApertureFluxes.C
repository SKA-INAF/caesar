
struct MacroOptions
{
	int bkgEstimator;
	bool dilateSMask;
	int dilateKernSize;
	int boxThicknessMin;
	int boxThicknessMax;
	int boxThicknessStep;

	MacroOptions()
	{
		SetDefaults();
	}

	void SetDefaults()
	{
		bkgEstimator= eMedianBkg;
		dilateSMask= false;
		dilateKernSize= 7;
		boxThicknessMin= 2;
		boxThicknessMax= 10;
		boxThicknessStep= 1;
	}

};

std::string sourceName="S172";

/*
tileMinX = 600                       
tileMaxX = 9600                                                                                         
tileMinY = 700                                                                                         
tileMaxY = 9000
*/

int ComputePhotApertureFluxes(std::string filename,std::string imgfilename,int tileMinX=-1,int tileMaxX=-1,int tileMinY=-1,int tileMaxY=-1,std::string sname_sel="",MacroOptions options=MacroOptions())
{
	//## Read input image file
	Image* img= new Image;
	if(img->ReadFITS(imgfilename,1,tileMinX,tileMaxX,tileMinY,tileMaxY)<0){
		cerr<<"ERROR: Failed to read image "<<imgfilename<<endl;
		return -1;
	}


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

		//if(aSource->GetName()!=sourceName) continue;//DEBUG!!!
		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);

	}//end loop sources


	//## Compute source mask	
	bool isBinary= true;
	bool invert= true;
	Image* smask= img->GetSourceMask(sources,isBinary,invert);
	if(!smask){
		cerr<<"ERROR: Failed to compute source mask!"<<endl;
		return -1;
	}

	
	//## Dilate mask
	Image* smask_final= 0;
	if(options.dilateSMask){
		int niters= 1;
		bool skipZeroPixels= false;
		smask_final= smask->GetMorphDilatedImage(options.dilateKernSize,niters,skipZeroPixels);
		if(!smask_final){
			cerr<<"ERROR: Failed to dilate source mask!"<<endl;
			return -1;
		}
	}
	else{
		smask_final= smask;
	}

	//## Loop over sources and compute bkg info around them
	for(size_t i=0;i<sources.size();i++)	
	{
		Source* source= sources[i];
		std::string sname= source->GetName();
		if(sname_sel!="" && sname!=sname_sel) continue;//skip source we are not interested

		double nPix= source->GetNPixels();	
		double bkgRMSSum= source->GetBkgRMSSum();
		double bkgSum= source->GetBkgSum();	
		double S= source->GetS();
		double beamArea= source->GetBeamFluxIntegral();
		double flux_nobkgSubtr= S/beamArea;
		double flux= (S-bkgSum)/beamArea;
		bool hasFitInfo= source->HasFitInfo();
		double flux_fit= 0;
		double fluxDensity_fit= 0;
		if(hasFitInfo){
			source->GetFluxDensity(fluxDensity_fit);
			flux_fit= fluxDensity_fit/beamArea;
		}
	
		//Compute flux with phot aperture method
		std::vector<double> fluxes_apert;
		for(int boxThickness=options.boxThicknessMin;boxThickness<=options.boxThicknessMax;boxThickness+=options.boxThicknessStep)
		{
			
			BkgSampleData bkgSampleData;
			img->GetBkgInfoAroundSource(bkgSampleData,source,boxThickness,options.bkgEstimator,smask_final);

			double bkgLevel= bkgSampleData.bkgLevel;
			double bkgRMS= bkgSampleData.bkgRMS;
			double nPix_bkg= bkgSampleData.npix;

			double S_apert= (S-nPix*bkgLevel);
			double flux_apert= S_apert/beamArea;
			fluxes_apert.push_back(flux_apert);

			cout<<"Source "<<sname<<": flux(mJy)="<<flux*1.e+3<<" (no bkg sub="<<flux_nobkgSubtr*1.e+3<<"), flux_fit(mJy)="<<flux_fit*1.e+3<<", flux_apert(mJy)="<<flux_apert*1.e+3<<", nPix_apert="<<nPix_bkg<<", bkgLevel_apert(mJy)="<<bkgLevel*1.e+3<<", nPix="<<nPix<<", bkgLevel="<<bkgSum/nPix*1.e+3<<endl;
		}

		double fluxMean_apert= 0;
		double fluxSigma_apert= 0;
		StatsUtils::ComputeMeanAndRMS(fluxMean_apert,fluxSigma_apert,fluxes_apert);
		double fluxMeanErr_apert= fluxSigma_apert/sqrt(fluxes_apert.size());
		cout<<"Source "<<sname<<": <flux_apert>(mJy)="<<fluxMean_apert*1.e+3<<" +- "<<fluxMeanErr_apert*1.e+3<<endl;
	
	}//end loop sources
	
	
	
	/*
	//DEBUG
	Image* tileImg= img->GetTile(1430-tileMinX,1520-tileMinX,5060-tileMinY,5140-tileMinY);
	tileImg->WriteFITS("prova.fits");

	Image* tileImg_mask= smask->GetTile(1430-tileMinX,1520-tileMinX,5060-tileMinY,5140-tileMinY);
	tileImg_mask->WriteFITS("prova2.fits");

	smask->WriteFITS("prova3.fits");

	Image* smask_noinvert= img->GetSourceMask(sources,false,false);
	smask_noinvert->WriteFITS("prova4.fits");
	*/

	return 0;

}//close macro
