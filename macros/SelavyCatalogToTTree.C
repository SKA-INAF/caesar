#include <iostream>

using namespace std;

//==========================================
//==     CLASS DEFINITIONS
//==========================================
struct IslandInfo 
{
	std::string islandId;
	std::string islandName;
	int nComponents;
	std::string x_hms;
	std::string y_hms;
	double x;//deg
	double y;//deg
	double freq;//MHz
	double bmaj_wcs;//arcsec
	double bmin_wcs;//arcsec
	double pa_wcs;//deg
	double flux;//mJy
	double flux_err;//mJy
	double peakFlux;//mJy/beam
	double avgBkg;
	double noiseRMS;
	double maxResidual;
	double minResidual;
	double meanResidual;
	double rmsResidual;
	double stdevResidual;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double npix;
	double solidAngle;//arcmin^2
	double beamArea;
	double x_avg;
	double y_avg;
	double x_cen;
	double y_cen;
	double x_peak;
	double y_peak;
	int flag_i1;
	int flag_i2;
	int flag_i3;
	int flag_i4;

	IslandInfo()
	{	
		islandId= "";
		islandName= "";
		nComponents= 0;
		x_hms= "";
		y_hms= "";
		x= 0;
		y= 0;
		freq= 0;
		bmaj_wcs= 0;
		bmin_wcs= 0;
		pa_wcs= 0;
		flux= 0;
		flux_err= 0;
		peakFlux= 0;
		avgBkg= 0;	
		noiseRMS= 0;
		maxResidual= 0;
		minResidual= 0;
		meanResidual= 0;
		rmsResidual= 0;
		stdevResidual= 0;
		xmin= 0;
		xmax= 0;
		ymin= 0;
		ymax= 0;
		npix= 0;
		solidAngle= 0;
		beamArea= 0;
		x_avg= 0;
		y_avg= 0;
		x_cen= 0;
		y_cen= 0;
		x_peak= 0;
		y_peak= 0;
		flag_i1= 0;
		flag_i2= 0;
		flag_i3= 0;
		flag_i4= 0;
	}

};//close IslandInfo


struct SourceInfo 
{
	std::string islandId;
	std::string componentId;
	std::string componentName;
	std::string x_hms;
	std::string y_hms;
	double x;//deg
	double y;//deg
	double x_err;//arcsec
	double y_err;//arcsec
	double freq;//MHz
	double peakFlux;//mJy/beam
	double peakFlux_err;//mJy/beam
	double flux;//mJy
	double flux_err;//mJy
	double bmaj_wcs;//arcsec
	double bmin_wcs;//arcsec
	double pa_wcs;//deg
	double bmaj_wcs_err;//arcsec
	double bmin_wcs_err;//arcsec
	double pa_wcs_err;//deg
	double bmaj_deconv_wcs;//arcsec
	double bmin_deconv_wcs;//arcsec
	double pa_deconv_wcs;//deg
	double bmaj_deconv_wcs_err;//arcsec
	double bmin_deconv_wcs_err;//arcsec
	double pa_deconv_wcs_err;//deg
	double chi2;
	double rms_fit;//mJy/beam
	double spectralIndex;
	double spectralIndex_err;
	double spectralCurv;
	double spectralCurv_err;
	double rms;//mJy/beam
	int hasSiblings;
	int fitFailed;
	int spectralIndexFromTT;
	int c4Flag;
	double snr;	
	
	

	SourceInfo()
	{
		islandId= "";
		componentId= "";
		componentName= "";
		x_hms= "";
		y_hms= "";
		x= 0;
		y= 0;
		x_err= 0;
		y_err= 0;
		freq= 0;
		peakFlux= 0;
		peakFlux_err= 0;
		flux= 0;
		flux_err= 0;
		bmaj_wcs= 0;
		bmin_wcs= 0;
		pa_wcs= 0;
		bmaj_wcs_err= 0;
		bmin_wcs_err= 0;
		pa_wcs_err= 0;
		bmaj_deconv_wcs= 0;
		bmin_deconv_wcs= 0;
		pa_deconv_wcs= 0;
		bmaj_deconv_wcs_err= 0;
		bmin_deconv_wcs_err= 0;
		pa_deconv_wcs_err= 0;
		chi2= 0;
		rms_fit= 0;
		spectralIndex= 0;
		spectralIndex_err= 0;
		spectralCurv= 0;
		spectralCurv_err= 0;
		hasSiblings= 0;
		rms= 0;
		fitFailed= 0;
		spectralIndexFromTT= 0;
		c4Flag= 0;

		snr= 0;		
	}


};//close struct SourceInfo

//==========================================
//==     METHODS
//==========================================
int ReadIslandData(std::string filename, std::vector<std::string> excludedPatterns= {});
int ReadData(std::string filename, std::vector<std::string> excludedPatterns= {});
static int ParseIslandData(const std::vector<std::string>& fields);
static int ParseData(const std::vector<std::string>& fields);
bool HasPattern(std::string,std::string);
void Init();
void Save();

//==========================================
//==     VARS
//==========================================
const int gNDataCols_selavy= 37;
const int gNDataCols_selavy_island= 38;
std::vector<std::string> gSkipLinePatterns_selavy {};
TFile* gOutputFile= 0;
std::string gOutputFileName= "";
TTree* gOutputTree= 0;
SourceInfo gSourceInfo;
IslandInfo gIslandInfo;
bool gReadIslandData= false;

int SelavyCatalogToTTree(std::string fileName,std::string fileName_out="output.root",bool readIslandData=false)
{
	//================================
	//==      CHECK ARGS
	//================================
	if(fileName=="" || fileName_out==""){
		cerr<<"ERROR: Empty input/output filename given!"<<endl;
		return -1;
	}
	gOutputFileName= fileName_out;
	gReadIslandData= readIslandData;

	//================================
	//==      INIT
	//================================
	cout<<"INFO: Init macro data ..."<<endl;
	Init();

	//================================
	//==      READ CATALOG DATA
	//================================
	if(gReadIslandData) ReadIslandData(fileName,gSkipLinePatterns_selavy);	
	else ReadData(fileName,gSkipLinePatterns_selavy);	

	//================================
	//==      SAVE
	//================================
	Save();

	return 0;

}//close macro


void Save()
{
	//Save TTree to file
	if(gOutputFile && gOutputFile->IsOpen()){
		gOutputFile->cd();
		if(gOutputTree) gOutputTree->Write();
		gOutputFile->Close();
	}

}//close Save()

void Init()
{
	//Create output file
	gOutputFile= new TFile(gOutputFileName.c_str(),"RECREATE");
	gOutputFile->cd();

	//Create output TTree
	gOutputTree= new TTree("data","data");

	if(gReadIslandData){
		gOutputTree->Branch("islandId",&gIslandInfo.islandId);
		gOutputTree->Branch("islandName",&gIslandInfo.islandName);
		gOutputTree->Branch("nComponents",&gIslandInfo.nComponents);
		gOutputTree->Branch("x_hms",&gIslandInfo.x_hms);
		gOutputTree->Branch("y_hms",&gIslandInfo.y_hms);
		gOutputTree->Branch("x",&gIslandInfo.x);
		gOutputTree->Branch("y",&gIslandInfo.y);
		gOutputTree->Branch("freq",&gIslandInfo.freq);
		gOutputTree->Branch("bmaj_wcs",&gIslandInfo.bmaj_wcs);
		gOutputTree->Branch("bmin_wcs",&gIslandInfo.bmin_wcs);
		gOutputTree->Branch("pa_wcs",&gIslandInfo.pa_wcs);
		gOutputTree->Branch("flux",&gIslandInfo.flux);
		gOutputTree->Branch("flux_err",&gIslandInfo.flux_err);	
		gOutputTree->Branch("peakFlux",&gIslandInfo.peakFlux);
		gOutputTree->Branch("avgBkg",&gIslandInfo.avgBkg);
		gOutputTree->Branch("noiseRMS",&gIslandInfo.noiseRMS);
		gOutputTree->Branch("maxResidual",&gIslandInfo.maxResidual);
		gOutputTree->Branch("minResidual",&gIslandInfo.minResidual);
		gOutputTree->Branch("meanResidual",&gIslandInfo.meanResidual);
		gOutputTree->Branch("rmsResidual",&gIslandInfo.rmsResidual);
		gOutputTree->Branch("stdevResidual",&gIslandInfo.stdevResidual);
		gOutputTree->Branch("xmin",&gIslandInfo.xmin);
		gOutputTree->Branch("xmax",&gIslandInfo.xmax);
		gOutputTree->Branch("ymin",&gIslandInfo.ymin);
		gOutputTree->Branch("ymax",&gIslandInfo.ymax);
		gOutputTree->Branch("npix",&gIslandInfo.npix);
		gOutputTree->Branch("solidAngle",&gIslandInfo.solidAngle);
		gOutputTree->Branch("beamArea",&gIslandInfo.beamArea);
		gOutputTree->Branch("x_avg",&gIslandInfo.x_avg);
		gOutputTree->Branch("y_avg",&gIslandInfo.y_avg);
		gOutputTree->Branch("x_cen",&gIslandInfo.x_cen);
		gOutputTree->Branch("y_cen",&gIslandInfo.y_cen);
		gOutputTree->Branch("x_peak",&gIslandInfo.x_peak);
		gOutputTree->Branch("y_peak",&gIslandInfo.y_peak);
		gOutputTree->Branch("flag_i1",&gIslandInfo.flag_i1);
		gOutputTree->Branch("flag_i2",&gIslandInfo.flag_i2);
		gOutputTree->Branch("flag_i3",&gIslandInfo.flag_i3);
		gOutputTree->Branch("flag_i4",&gIslandInfo.flag_i4);

	}//close if read island data
	else{
		gOutputTree->Branch("islandId",&gSourceInfo.islandId);
		gOutputTree->Branch("componentId",&gSourceInfo.componentId);
		gOutputTree->Branch("componentName",&gSourceInfo.componentName);
		gOutputTree->Branch("x_hms",&gSourceInfo.x_hms);
		gOutputTree->Branch("y_hms",&gSourceInfo.y_hms);
		gOutputTree->Branch("x",&gSourceInfo.x);
		gOutputTree->Branch("y",&gSourceInfo.y);
		gOutputTree->Branch("x_err",&gSourceInfo.x_err);
		gOutputTree->Branch("y_err",&gSourceInfo.y_err);
		gOutputTree->Branch("freq",&gSourceInfo.freq);
		gOutputTree->Branch("peakFlux",&gSourceInfo.peakFlux);
		gOutputTree->Branch("peakFlux_err",&gSourceInfo.peakFlux_err);
		gOutputTree->Branch("flux",&gSourceInfo.flux);
		gOutputTree->Branch("flux_err",&gSourceInfo.flux_err);	
		gOutputTree->Branch("snr",&gSourceInfo.snr);
		gOutputTree->Branch("bmaj_wcs",&gSourceInfo.bmaj_wcs);
		gOutputTree->Branch("bmin_wcs",&gSourceInfo.bmin_wcs);
		gOutputTree->Branch("pa_wcs",&gSourceInfo.pa_wcs);
		gOutputTree->Branch("bmaj_wcs_err",&gSourceInfo.bmaj_wcs_err);
		gOutputTree->Branch("bmin_wcs_err",&gSourceInfo.bmin_wcs_err);
		gOutputTree->Branch("pa_wcs_err",&gSourceInfo.pa_wcs_err);
		gOutputTree->Branch("bmaj_deconv_wcs",&gSourceInfo.bmaj_deconv_wcs);
		gOutputTree->Branch("bmin_deconv_wcs",&gSourceInfo.bmin_deconv_wcs);
		gOutputTree->Branch("pa_deconv_wcs",&gSourceInfo.pa_deconv_wcs);
		gOutputTree->Branch("bmaj_deconv_wcs_err",&gSourceInfo.bmaj_deconv_wcs_err);
		gOutputTree->Branch("bmin_deconv_wcs_err",&gSourceInfo.bmin_deconv_wcs_err);
		gOutputTree->Branch("pa_deconv_wcs_err",&gSourceInfo.pa_deconv_wcs_err);
		gOutputTree->Branch("chi2_fit",&gSourceInfo.chi2);
		gOutputTree->Branch("rms_fit",&gSourceInfo.rms_fit);
		gOutputTree->Branch("spectralIndex",&gSourceInfo.spectralIndex);
		gOutputTree->Branch("spectralIndex_err",&gSourceInfo.spectralIndex_err);
		gOutputTree->Branch("spectralCurv",&gSourceInfo.spectralCurv);
		gOutputTree->Branch("spectralCurv_err",&gSourceInfo.spectralCurv_err);
		gOutputTree->Branch("rms",&gSourceInfo.rms);
		gOutputTree->Branch("multicomponentFit",&gSourceInfo.hasSiblings);
		gOutputTree->Branch("fitFailed",&gSourceInfo.fitFailed);
		gOutputTree->Branch("spectralIndexFromTT",&gSourceInfo.spectralIndexFromTT);
		gOutputTree->Branch("c4Flag",&gSourceInfo.c4Flag);
	}

}//close Init()

int ReadData(std::string filename,std::vector<std::string> excludedPatterns)
{
	//Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		cerr<<"ERROR: Failed to open config file "<<filename<<" for reading!"<<endl;
		return -1;
  }

	//Parsing file
	cout<<"INFO: Reading and parsing file: "<<filename<<endl;

	std::string parsedline= "";	
	int line_counter= 0;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;

		//Parse line fields
		istringstream ss(parsedline);
		std::string field= "";
		bool skipLine= false;
		std::vector<std::string> fields;

		while(ss >> field)
    {
			//Skip first line
			char first_char= *(field.c_str());
			if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
				skipLine= true;
				break;
			}

			//Search for pattern and skip if present
			for(size_t i=0;i<excludedPatterns.size();i++){
				bool hasPattern= HasPattern(field,excludedPatterns[i]);
				if(hasPattern) {
					skipLine= true;
					break;
				}
			}
			if(skipLine) break;
	
			//cout<<field<<"  ";
			fields.push_back(field);

		}//end loop line
		
		if(skipLine) continue;
		//else cout<<endl;

		
		//Process line
		if(ParseData(fields)<0){
			cerr<<"ERROR: Failed to parse line "<<line_counter<<" of file "<<filename<<"!"<<endl;
		}	

		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	return 0;

}//close ReadData()





int ReadIslandData(std::string filename,std::vector<std::string> excludedPatterns)
{
	//Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		cerr<<"ERROR: Failed to open config file "<<filename<<" for reading!"<<endl;
		return -1;
  }

	//Parsing file
	cout<<"INFO: Reading and parsing file: "<<filename<<endl;

	std::string parsedline= "";	
	int line_counter= 0;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;

		//Parse line fields
		istringstream ss(parsedline);
		std::string field= "";
		bool skipLine= false;
		std::vector<std::string> fields;

		while(ss >> field)
    {
			//Skip first line
			char first_char= *(field.c_str());
			if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
				skipLine= true;
				break;
			}

			//Search for pattern and skip if present
			for(size_t i=0;i<excludedPatterns.size();i++){
				bool hasPattern= HasPattern(field,excludedPatterns[i]);
				if(hasPattern) {
					skipLine= true;
					break;
				}
			}
			if(skipLine) break;
	
			//cout<<field<<"  ";
			fields.push_back(field);

		}//end loop line
		
		if(skipLine) continue;
		//else cout<<endl;

		
		//Process line
		if(ParseIslandData(fields)<0){
			cerr<<"ERROR: Failed to parse line "<<line_counter<<" of file "<<filename<<"!"<<endl;
		}	

		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	return 0;

}//close ReadIslandData()

int ParseData(const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols_selavy){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols_selavy<<" expected!"<<endl;
		return -1;
	}

	//#   island_id    component_id component_name ra_hms_cont dec_dms_cont ra_deg_cont dec_deg_cont    ra_err    dec_err  freq  flux_peak flux_peak_err  flux_int flux_int_err maj_axis min_axis pos_ang maj_axis_err min_axis_err pos_ang_err maj_axis_deconv min_axis_deconv pos_ang_deconv maj_axis_deconv_err min_axis_deconv_err pos_ang_deconv_err chi_squared_fit rms_fit_gauss spectral_index spectral_curvature spectral_index_err spectral_curvature_err  rms_image has_siblings fit_is_estimate spectral_index_from_TT flag_c4
	
	//Set var values
	std::string blobName= std::string(fields[0]);
	std::string compId= std::string(fields[1]);
	std::string compName= std::string(fields[2]);
	std::string x_hms= std::string(fields[3]);
	std::string y_hms= std::string(fields[4]);
	double x_wcs= atof(fields[5].c_str());
	double y_wcs= atof(fields[6].c_str());
	double x_wcs_err= atof(fields[7].c_str());
	double y_wcs_err= atof(fields[8].c_str());
	double freq= atof(fields[9].c_str());
	double A= atof(fields[10].c_str())/1000.;//in Jy/beam
	double A_err= atof(fields[11].c_str())/1000.;//in Jy/beam
	double flux= atof(fields[12].c_str())/1000;//in Jy
	double flux_err= atof(fields[13].c_str())/1000;//in Jy

	double bmaj_wcs= atof(fields[14].c_str());
	double bmin_wcs= atof(fields[15].c_str());
	double pa_wcs= atof(fields[16].c_str());
	double bmaj_wcs_err= atof(fields[17].c_str());
	double bmin_wcs_err= atof(fields[18].c_str());
	double pa_wcs_err= atof(fields[19].c_str());

	double bmaj_deconv_wcs= atof(fields[20].c_str());
	double bmin_deconv_wcs= atof(fields[21].c_str());
	double pa_deconv_wcs= atof(fields[22].c_str());
	double bmaj_deconv_wcs_err= atof(fields[23].c_str());
	double bmin_deconv_wcs_err= atof(fields[24].c_str());
	double pa_deconv_wcs_err= atof(fields[25].c_str());

	double chi2= atof(fields[26].c_str());
	double rms_fit= atof(fields[27].c_str());
	double spectralIndex= atof(fields[28].c_str());
	double spectralCurv= atof(fields[29].c_str());
	double spectralIndex_err= atof(fields[30].c_str());
	double spectralCurv_err= atof(fields[31].c_str());
	
	double rms= atof(fields[32].c_str())/1000;//in Jy/beam
	int hasSiblings= atoi(fields[33].c_str());	
	int fitFailed= atoi(fields[34].c_str());
	int spectralIndexFromTT= atoi(fields[35].c_str());
	int c4Flag= atoi(fields[36].c_str());
	double snr= A/rms;
	

	//Set data struct
	gSourceInfo.islandId= blobName;
	gSourceInfo.componentId= compId;
	gSourceInfo.componentName= compName;
	gSourceInfo.x_hms= x_hms;
	gSourceInfo.y_hms= y_hms;
	gSourceInfo.x= x_wcs;
	gSourceInfo.y= y_wcs;
	gSourceInfo.x_err= x_wcs_err;
	gSourceInfo.y_err= y_wcs_err;
	gSourceInfo.freq= freq;
	gSourceInfo.peakFlux= A;
	gSourceInfo.peakFlux_err= A_err;
	gSourceInfo.flux= flux;
	gSourceInfo.flux_err= flux_err;
	gSourceInfo.bmaj_wcs= bmaj_wcs;
	gSourceInfo.bmin_wcs= bmin_wcs;
	gSourceInfo.pa_wcs= pa_wcs;
	gSourceInfo.bmaj_wcs_err= bmaj_wcs_err;
	gSourceInfo.bmin_wcs_err= bmin_wcs_err;
	gSourceInfo.pa_wcs_err= pa_wcs_err;

	gSourceInfo.bmaj_deconv_wcs= bmaj_deconv_wcs;
	gSourceInfo.bmin_deconv_wcs= bmin_deconv_wcs;
	gSourceInfo.pa_deconv_wcs= pa_deconv_wcs;
	gSourceInfo.bmaj_deconv_wcs_err= bmaj_deconv_wcs_err;
	gSourceInfo.bmin_deconv_wcs_err= bmin_deconv_wcs_err;
	gSourceInfo.pa_deconv_wcs_err= pa_deconv_wcs_err;
	gSourceInfo.chi2= chi2;
	gSourceInfo.rms_fit= rms_fit;
	
	gSourceInfo.spectralIndex= spectralIndex;
	gSourceInfo.spectralIndex_err= spectralIndex_err;
	gSourceInfo.spectralCurv= spectralCurv;
	gSourceInfo.spectralCurv_err= spectralCurv_err;
	
	gSourceInfo.rms= rms;
	gSourceInfo.hasSiblings= hasSiblings;
	gSourceInfo.fitFailed= fitFailed;
	gSourceInfo.spectralIndexFromTT= spectralIndexFromTT;
	gSourceInfo.c4Flag= c4Flag;

	gSourceInfo.snr= snr;
	

	//Fill TTree
	if(gOutputTree) gOutputTree->Fill();

	return 0;

}//close ParseData()




int ParseIslandData(const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols_selavy_island){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols_selavy<<" expected!"<<endl;
		return -1;
	}

	//#   island_id    island_name n_components ra_hms_cont dec_dms_cont ra_deg_cont dec_deg_cont    freq maj_axis min_axis pos_ang flux_int flux_int_err  flux_peak mean_background background_noise max_residual min_residual mean_residual rms_residual stdev_residual x_min x_max y_min y_max    n_pix solid_angle beam_area   x_ave   y_ave   x_cen   y_cen x_peak y_peak flag_i1 flag_i2 flag_i3 flag_i4

	//Set var values
	std::string islandId= std::string(fields[0]);
	std::string islandName= std::string(fields[1]);
	int nComponents= atoi(fields[2]);
	std::string x_hms= std::string(fields[3]);
	std::string y_hms= std::string(fields[4]);
	double x_wcs= atof(fields[5].c_str());
	double y_wcs= atof(fields[6].c_str());
	double freq= atof(fields[7].c_str());
	double bmaj_wcs= atof(fields[8].c_str());
	double bmin_wcs= atof(fields[9].c_str());
	double pa_wcs= atof(fields[10].c_str());	
	double flux= atof(fields[11].c_str())/1000;//in Jy
	double flux_err= atof(fields[12].c_str())/1000;//in Jy
	double A= atof(fields[13].c_str())/1000.;//in Jy/beam
	double avgBkg= atof(fields[14].c_str());
	double noiseRMS= atof(fields[15].c_str());

	double max_residual= atof(fields[16].c_str());
	double min_residual= atof(fields[17].c_str());
	double mean_residual= atof(fields[18].c_str());
	double rms_residual= atof(fields[19].c_str());
	double stdev_residual= atof(fields[20].c_str());
	double xmin= atof(fields[21].c_str());
	double xmax= atof(fields[22].c_str());
	double ymin= atof(fields[23].c_str());
	double ymax= atof(fields[24].c_str());
	double npix= atof(fields[25].c_str());
	double solidAngle= atof(fields[26].c_str());
	double beamArea= atof(fields[27].c_str());
	
	double x_ave= atof(fields[28].c_str());
	double y_ave= atof(fields[29].c_str());
	double x_cen= atof(fields[30].c_str());
	double y_cen= atof(fields[31].c_str());
	double x_peak= atof(fields[32].c_str());
	double y_peak= atof(fields[33].c_str());

	int flag_i1= atoi(fields[34].c_str());
	int flag_i2= atoi(fields[35].c_str());
	int flag_i3= atoi(fields[36].c_str());
	int flag_i4= atoi(fields[37].c_str());

	//Set data struct
	gIslandInfo.islandId= islandId;
	gIslandInfo.islandName= islandName;
	gIslandInfo.nComponents= nComponents;
	gIslandInfo.x_hms= x_hms;
	gIslandInfo.y_hms= y_hms;
	gIslandInfo.x= x;
	gIslandInfo.y= y;
	gIslandInfo.freq= freq;
	gIslandInfo.bmaj_wcs= bmaj_wcs
	gIslandInfo.bmin_wcs= bmin_wcs
	gIslandInfo.pa_wcs= pa_wcs;
	gIslandInfo.flux= flux;
	gIslandInfo.flux_err= flux_err;
	gIslandInfo.peakFlux= A;
	gIslandInfo.avgBkg= avgBkg;	
	gIslandInfo.noiseRMS= noiseRMS;
	gIslandInfo.maxResidual= max_residual;
	gIslandInfo.minResidual= min_residual;
	gIslandInfo.meanResidual= mean_residual;
	gIslandInfo.rmsResidual= rms_residual;
	gIslandInfo.stdevResidual= stdev_residual;
	gIslandInfo.xmin= xmin;
	gIslandInfo.xmax= xmax;
	gIslandInfo.ymin= ymin;
	gIslandInfo.ymax= ymax;
	gIslandInfo.npix= npix;
	gIslandInfo.solidAngle= solidAngle;
	gIslandInfo.beamArea= beamArea;
	gIslandInfo.x_avg= x_avg;
	gIslandInfo.y_avg= y_avg;
	gIslandInfo.x_cen= x_cen;
	gIslandInfo.y_cen= y_cen;
	gIslandInfo.x_peak= x_peak;
	gIslandInfo.y_peak= y_peak;
	gIslandInfo.flag_i1= flag_i1;
	gIslandInfo.flag_i2= flag_i2;
	gIslandInfo.flag_i3= flag_i3;
	gIslandInfo.flag_i4= flag_i4;

	//Fill TTree
	if(gOutputTree) gOutputTree->Fill();

	return 0;

}//close ParseIslandData()

bool HasPattern(std::string str,std::string pattern)
{
	std::size_t found = str.find(pattern);
  if (found!=std::string::npos) return true;

	return false;

}//close HasPattern()
