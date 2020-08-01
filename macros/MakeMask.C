

int MakeMask(std::string filename,double seedThr=5,double mergeThr=2.6,int minNPix=5,double deltaThr=0,bool draw=false)
{
	//Read FITS
	cout<<"Read FITS..."<<endl;
	Image* img= new Image;
	if(img->ReadFITS(filename)<0){
		cerr<<"ERROR: Failed to read FITS"<<endl;
		return -1;
	}

	//Clone image 
	Image* img_tmp= img->GetCloned("img_tmp");

	//Compute stats
	cout<<"Compute stats..."<<endl;
	img_tmp->ComputeStats(true);
	img_tmp->PrintStats();

	//Compute background
	cout<<"Compute bkg ..."<<endl;
	ImgBkgData* bkgData= img_tmp->ComputeBkg(eMedianBkg,false);
	//ImgBkgData* bkgData= img->ComputeBkg(eMedianBkg,true,20,20,5,5,true,true);

	
	//Find significance map
	cout<<"Compute significance map ..."<<endl;
	Image* zmap= img_tmp->GetSignificanceMap(bkgData,false);

	//Find source (iter 0)
	cout<<"Find sources (iter0)"<<endl;
	std::vector<Source*> sources_iter0;
	img_tmp->FindCompactSource(sources_iter0,zmap,bkgData,seedThr,mergeThr,minNPix);
	
	//Mask sources
	cout<<"Mask sources ..."<<endl;
	double maskValue= 0.;
	img_tmp->MaskSources(sources_iter0,maskValue);	

	//Recompute stats excluding masked pixels
	cout<<"Recompute stats excluding masked pixels..."<<endl;
	img_tmp->ComputeStats(true,true);
	img_tmp->PrintStats();

	//Recompute bkg
	cout<<"Recompute bkg..."<<endl;
	ImgBkgData* bkgData_noSources= img_tmp->ComputeBkg(eMedianBkg,false);
	
	//Recompute zmap from initial map using recomputed bkg
	cout<<"Recompute significance..."<<endl;
	Image* zmap_noSources= img->GetSignificanceMap(bkgData_noSources,false);

	//Compute sources with improved bkg
	cout<<"Find sources ..."<<endl;
	std::vector<Source*> sources;
	double seedThr_iter= seedThr-deltaThr;
	img->FindCompactSource(sources,zmap_noSources,bkgData_noSources,seedThr_iter,mergeThr,minNPix);
	
	//Set source type
	for(size_t i=0;i<sources.size();i++)
	{
		sources[i]->SetType(eCompact);
	}

	//Get mask
	cout<<"Get mask ..."<<endl;
	Image* maskImg= img->GetSourceMask(sources,true);
	//Image* maskImg= img->GetSourceMask(sources,false);

		
	//Plot sources
	if(draw)
	{
		cout<<"Plot sources..."<<endl;
		img->Plot(sources);

		//Plot mask
		TCanvas* MaskPlot= new TCanvas("MaskPlot","MaskPlot");
		MaskPlot->cd();

		TH2D* maskHisto= maskImg->GetHisto2D("maskHisto");
		maskHisto->SetStats(0);
		maskHisto->Draw("COLZ");

		//Plot significance map
		TCanvas* ZMapPlot= new TCanvas("ZMapPlot","ZMapPlot");
		ZMapPlot->cd();

		TH2D* zmapHisto= zmap_noSources->GetHisto2D("zmapHisto");
		zmapHisto->SetStats(0);
		zmapHisto->Draw("COLZ");
	}
	

	//Write mask to file
	std::string baseFileName= CodeUtils::ExtractFileNameFromPath(filename,true);
	std::string maskOutFileName= baseFileName + "_mask.fits";
	cout<<"Writing mask file "<<maskOutFileName<<" ..."<<endl;
	
	maskImg->WriteFITS(maskOutFileName);

	//Write regions
	std::string regionOutFile= baseFileName + "_sources.reg";
	cout<<"Writing DS9 region file "<<regionOutFile<<" ..."<<endl;
	SourceExporter::WriteToDS9(regionOutFile,sources);

	
	return 0;

}//close MakeMask()
