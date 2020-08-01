

int DrawSourcePosDiff(std::string filename)
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

	Source* source= 0;
	SourceInfo->SetBranchAddress("Source",&source);

	cout<<"# "<<SourceInfo->GetEntries()<<" sources to be drawn..."<<endl;
	
	//## Store source pos list 
	//std::vector<Source*> sources;
	WCS* wcs= 0;
	std::vector<TVector2> sourcePosList;
	std::vector<TVector2> sourceCompPosList;

	for(int i=0;i<SourceInfo->GetEntries();i++)
	{
		SourceInfo->GetEntry(i);

		//Get WCS
		if(!wcs){
			wcs= source->GetWCS();
			if(!wcs){
				cerr<<"ERROR: Cannot get WCS!"<<endl;
				return -1;
			}
		}

		//Get source position in WCS
		double x0= source->X0;
		double y0= source->Y0;
		double x0_wcs= 0;
		double y0_wcs= 0;
		AstroUtils::PixelToWCSCoords(x0_wcs,y0_wcs,wcs,x0,y0);

		sourcePosList.push_back(TVector2(x0_wcs,y0_wcs));

		//Get fit pars
		bool hasFitPars= source->HasFitInfo();
		if(!hasFitPars) continue;
		SourceFitPars fitPars= source->GetFitPars();
		int nComponents= source->GetNFitComponents();

		for(int j=0;j<nComponents;j++)
		{	
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(j);
			if(!isSelected) continue;

			//Get component WCS position
			x0= 0;
			y0= 0;	
			fitPars.GetComponentPosition(x0,y0,j);
			x0_wcs= 0;
			y0_wcs= 0;
			AstroUtils::PixelToWCSCoords(x0_wcs,y0_wcs,wcs,x0,y0); 

			sourcePosList.push_back(TVector2(x0_wcs,y0_wcs));

		}//end loop components 		
	}//end loop sources


	//Compute source mutual difference
	std::vector<double> posDiffList;

	for(size_t i=0;i<sourcePosList.size();i++)
	{
		if(i%1000==0) cout<<"--> Source no. "<<i+1<<"/"<<sourcePosList.size()<<" processed..."<<endl;
		TVector2 pos_i= sourcePosList[i];
		double x_i= pos_i.X();
		double y_i= pos_i.Y();

		for(size_t j=i+1;j<sourcePosList.size();j++)
		{
			TVector2 pos_j= sourcePosList[j];
			double x_j= pos_j.X();
			double y_j= pos_j.Y();

			double dpos= sqrt( pow(x_i-x_j,2) + pow(y_i-y_j,2) )*3600.;
			posDiffList.push_back(dpos);
			
		}//end loop sources
	}//end loop sources

	//Compute stats
	Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(posDiffList);
	double medianPos= stats.median;
	double minPos= stats.minVal;
	double maxPos= stats.maxVal;
	double iqrPos= stats.Q3-stats.Q1;

	cout<<"source position diff stats: median="<<medianPos<<", min/max(arcsec)="<<minPos<<"/"<<maxPos<<", iqr="<<iqrPos<<endl;

	//Fill histo
	TH1D* posDiffHisto= new TH1D("posDiffHisto","",100,minPos-1,maxPos+1);
	posDiffHisto->Sumw2();
	
	for(size_t i=0;i<posDiffList.size();i++)
	{
		posDiffHisto->Fill(posDiffList[i]);

	}//end loop data

	TCanvas* Plot= new TCanvas("Plot","Relative source difference in arcsec");
	Plot->cd();

	posDiffHisto->GetXaxis()->SetTitle("r (arcsec)");
	posDiffHisto->DrawNormalized("hist");

	

	return 0;

}//close DrawSourcePosDiff()
