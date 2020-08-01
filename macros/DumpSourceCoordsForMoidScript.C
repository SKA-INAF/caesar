

int DumpSourceCoordsForMoidScript(std::string filename="",std::string outfilename="sources_moid.dat",std::string outfilename_comp="sourceComponents_moid.dat")
{
	//- Open output file
	FILE* fout= fopen(outfilename.c_str(),"w");
	FILE* fout_comp= fopen(outfilename_comp.c_str(),"w");

	//- Create files for single and multi-components
	std::string baseFileName= CodeUtils::ExtractFileNameFromPath(outfilename,true);
	std::string outfilename_single= baseFileName + std::string("_single.dat");
	std::string outfilename_multiple= baseFileName + std::string("_multiple.dat");
	FILE* fout_single= fopen(outfilename_single.c_str(),"w");
	FILE* fout_multiple= fopen(outfilename_multiple.c_str(),"w");

	std::string baseFileName_comp= CodeUtils::ExtractFileNameFromPath(outfilename_comp,true);
	std::string outfilename_comp_single= baseFileName_comp + std::string("_single.dat");
	std::string outfilename_comp_multiple= baseFileName_comp + std::string("_multiple.dat");
	FILE* fout_comp_single= fopen(outfilename_comp_single.c_str(),"w");
	FILE* fout_comp_multiple= fopen(outfilename_comp_multiple.c_str(),"w");

	//- Read input file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()) {
		cerr<<"ERROR: Cannot open file "<<filename<<"!"<<endl;
		return -1;
	}
	
	//- Get source tree
	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");	
	if(!SourceInfo){
		cerr<<"ERROR: Cannot get sources from file "<<filename<<"!"<<endl;
		return -1;
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	cout<<"INFO: # "<<SourceInfo->GetEntries()<<" sources to be read..."<<endl;
	
	//- Read source data
	WCS* wcs_gal= 0;
	double rmax= -1.e-99;
	std::string sname_rmax= "";
	double x0_rmax= 0;
	double y0_rmax= 0;	
	
	double rmax_comp= -1.e-99;
	std::string sname_rmax_comp= "";
	double x0_rmax_comp= 0;
	double y0_rmax_comp= 0;	

	double rmax_tol= 0.2;//add tolerance to max radius
	

	for(int i=0;i<SourceInfo->GetEntries();i++)
	{
		SourceInfo->GetEntry(i);

		//Check if has selected fit components
		if(aSource->GetNSelFitComponents()<=0) continue;

		//Get galactic WCS
		if(!wcs_gal){
			wcs_gal= aSource->GetWCS(eGALACTIC);
			if(!wcs_gal){
				cerr<<"ERROR: Cannot get WCS!"<<endl;
				return -1;
			}
		}

		//Get name
		std::string name= aSource->GetName();

		//Get coordinates
		double x0= aSource->X0;
		double y0= aSource->Y0;
		float xmin= 0;
		float xmax= 0;
		float ymin= 0;
		float ymax= 0;
		aSource->GetSourceRange(xmin,xmax,ymin,ymax);

		//Get WCS coordinates
		double l_wcs= 0;
		double b_wcs= 0;
		AstroUtils::PixelToWCSCoords(l_wcs,b_wcs,wcs_gal,x0,y0); 
	
		double l1_wcs= 0;
		double b1_wcs= 0;
		AstroUtils::PixelToWCSCoords(l1_wcs,b1_wcs,wcs_gal,xmin,ymin);
		double l2_wcs= 0;
		double b2_wcs= 0;
		AstroUtils::PixelToWCSCoords(l2_wcs,b2_wcs,wcs_gal,xmax,ymax); 

		double dist= sqrt( pow(l2_wcs-l_wcs,2) + pow(b2_wcs-b_wcs,2) );
		if(dist>rmax){
			rmax= dist;
			sname_rmax= name;
			x0_rmax= x0;
			y0_rmax= y0;
		}

		//Get components
		bool hasFitPars= aSource->HasFitInfo();
		if(!hasFitPars) continue;
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= aSource->GetNFitComponents();
		for(int j=0;j<nComponents;j++)
		{	
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(j);
			if(!isSelected) continue;

			//Get component WCS position
			x0= 0;
			y0= 0;
			fitPars.GetComponentPosition(x0,y0,j);
			l_wcs= 0;
			b_wcs= 0;
			AstroUtils::PixelToWCSCoords(l_wcs,b_wcs,wcs_gal,x0,y0);

			//Get WCS ellipse pars
			double x0_wcs= 0;
			double y0_wcs= 0;
			double bmaj_wcs= 0;
			double bmin_wcs= 0;
			double pa_wcs= 0;
			fitPars.GetComponentFitWCSEllipsePars(j,x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs);

			if(bmaj_wcs>rmax_comp){
				rmax_comp= bmaj_wcs;
				sname_rmax_comp= name;
				x0_rmax_comp= x0;
				y0_rmax_comp= y0;
			}

		}//end loop components
		
	}//end loop

	//- Add tolerance
	//rmax*= (1+rmax_tol);
	rmax_comp*= (1-rmax_tol);

	cout<<"INFO: Island rmax name="<<sname_rmax<<", (x,y)=("<<x0_rmax<<","<<y0_rmax<<"), rmax(deg)="<<rmax<<", rmax(arcsec)="<<rmax*3600<<endl;
	cout<<"INFO: Component rmax name="<<sname_rmax_comp<<" (x,y)=("<<x0_rmax_comp<<","<<y0_rmax_comp<<"), rmax(arcsec)="<<rmax_comp<<endl;

	
	for(int i=0;i<SourceInfo->GetEntries();i++)
	{
		SourceInfo->GetEntry(i);

		//Check if has selected fit components
		if(aSource->GetNSelFitComponents()<=0) continue;

		//Get galactic WCS
		if(!wcs_gal){
			wcs_gal= aSource->GetWCS(eGALACTIC);
			if(!wcs_gal){
				cerr<<"ERROR: Cannot get WCS!"<<endl;
				return -1;
			}
		}

		//Get name
		std::string name= aSource->GetName();

		//Get WCS coordinates
		double x0= aSource->X0;
		double y0= aSource->Y0;
		double l_wcs= 0;
		double b_wcs= 0;
		AstroUtils::PixelToWCSCoords(l_wcs,b_wcs,wcs_gal,x0,y0); 

		//Get components
		bool hasFitPars= aSource->HasFitInfo();
		if(!hasFitPars) continue;
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= aSource->GetNFitComponents();
		int nSelComponents= aSource->GetNSelFitComponents();

		//Print to file
		fprintf(fout,"%f %f %f %s %s\n",l_wcs,b_wcs,rmax*3600,"gal",name.c_str());
		if(nSelComponents==1) fprintf(fout_single,"%f %f %f %s %s\n",l_wcs,b_wcs,rmax*3600,"gal",name.c_str());
		else fprintf(fout_multiple,"%f %f %f %s %s\n",l_wcs,b_wcs,rmax*3600,"gal",name.c_str());

		
		for(int j=0;j<nComponents;j++)
		{	
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(j);
			if(!isSelected) continue;

			//Get name
			std::string compId= std::to_string(j+1);
			std::string cname= name + std::string("_fitcomp") + compId;

			//Get component WCS position
			x0= 0;
			y0= 0;
			fitPars.GetComponentPosition(x0,y0,j);
			l_wcs= 0;
			b_wcs= 0;
			AstroUtils::PixelToWCSCoords(l_wcs,b_wcs,wcs_gal,x0,y0);

			//Print to file
			fprintf(fout_comp,"%f %f %f %s %s\n",l_wcs,b_wcs,rmax_comp,"gal",cname.c_str());

			if(nSelComponents==1) fprintf(fout_comp_single,"%f %f %f %s %s\n",l_wcs,b_wcs,rmax_comp,"gal",cname.c_str());
			else fprintf(fout_comp_multiple,"%f %f %f %s %s\n",l_wcs,b_wcs,rmax_comp,"gal",cname.c_str());


		}//end loop components

	}//end loop sources


	//- Close file
	fclose(fout);
	fclose(fout_comp);
	fclose(fout_comp_single);
	fclose(fout_comp_multiple);	

	return 0;

}//close macro
