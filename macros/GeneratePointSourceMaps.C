

int GeneratePointSourceMaps(std::string filename,int nGenImgs=1,double sourceDensity=100,double Smin=1.e-4,double Smax=1,double Sslope=1.6)
{
	//- Read input map
	Image* img= new Image;
	if(img->ReadFITS(filename)<0){
		cerr<<"ERROR: Failed to read input map "<<filename<<"!"<<endl;
		return -1;
	}

	//- Set output filenames
	//std::string filename_base= CodeUtils::ExtractFileNameFromPath(filename,true);
	std::string outputfile_prefix= "simmap";
	std::string outputfile_info_prefix= "sources";

	
	for(int i=0;i<nGenImgs;i++)
	{
		//Clone map
		Image* simimg= img->GetCloned("",true,false);
		
		//Add point sources
		cout<<"INFO: Adding point sources (RUN "<<i+1<<") ..."<<endl;
		std::string outputfile_info= Form("%s-RUN%d.dat",outputfile_info_prefix.c_str(),i+1); 
		simimg->AddSimPointSourceDensity(sourceDensity,Smin,Smax,Sslope,0,0,outputfile_info);
		
		//Write FITS file
		std::string outputfile= Form("%s-RUN%d.fits",outputfile_prefix.c_str(),i+1);		
		cout<<"INFO: Writing outfile "<<outputfile<<" ..."<<endl;
		simimg->WriteFITS(outputfile);

		//Delete map	
		cout<<"DEBUG: Deleting sim map..."<<endl;	
		delete simimg;
		simimg= 0;	
	

	}//end loop generated maps

	
	cout<<"END RUN"<<endl;


	return 0;

}//close macro
