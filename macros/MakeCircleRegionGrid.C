

int MakeCircleRegionGrid(double ra_start,double ra_stop,double dra_arcmin,double dec_start,double dec_stop,double ddec_arcmin,double radius_arcmin,std::string outputFile="ds9.reg")
{
	//Open output file
	FILE* fout= fopen(outputFile.c_str(),"w");
	
	//Write header
	fprintf(fout,"global color=green font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"fk5\n");

	//Loop over grid
	double dra= dra_arcmin/60.;//in deg
	double ddec= ddec_arcmin/60.;//in deg
	double radius= radius_arcmin*60.;//in arcsec
	double radius_deg= radius_arcmin/60.;//in deg
	int nsteps_ra= (ra_stop-ra_start)/dra;
	int nsteps_dec= (dec_stop-dec_start)/ddec;
	int counter= 1;

	for(int i=0;i<=nsteps_dec;i++)	
	{
		double dec= dec_start + i*ddec;
		double ra_0= ra_start;
		if(i%2==1) ra_0+= dra/2.;
			
		for(int j=0;j<=nsteps_ra;j++)
		{
			double ra= ra_0 + j*dra;
			
			std::stringstream ss;
			ss<<"circle("<<ra<<","<<dec<<","<<radius<<"\""<<") # text={R"<<counter<<"}";
			fprintf(fout,"%s\n",ss.str().c_str());

			counter++;
		}//end ra loop
	}//end dec loop

	//Close file
	fclose(fout);

	return 0;

}//close macros
