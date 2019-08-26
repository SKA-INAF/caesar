
using namespace std;
using namespace Caesar;

//Functions
int ApplySelection(std::string sourceFileName,std::string regionFileName);
int ReadList(std::vector<std::string>& fileNames, std::string filelist);

//Variables
std::vector<DS9Region*> m_regions;
const int gNDataCols= 6;

int SelectSimSourcesInMosaic(std::string sourceFileList,std::string regionFileList,std::string mosaicRegionFile)
{
	//##################################
	//##  PARSE REGIONS
	//##################################
	//Read and parse DS9 regions
	if (DS9RegionParser::Parse(m_regions,mosaicRegionFile)<0){
		cerr<<"ERROR: Failed to read and parse DS9 region file "<<mosaicRegionFile<<"!"<<endl;
		return -1;
	}

	//Check if empty
	if(m_regions.empty()){
		cerr<<"No regions read from file "<<mosaicRegionFile<<"!"<<endl;
		return -1;
	}
	else{
		cout<<"INFO: #"<<m_regions.size()<<" regions read from file "<<mosaicRegionFile<<"..."<<endl;
	}
	
	//##################################
	//##  PARSE REGION FILELIST
	//##################################
	std::vector<std::string> regionFileNames;
	if(ReadList(regionFileNames,regionFileList)<0){
		cerr<<"ERROR: Failed to read list "<<regionFileList<<endl;
		return -1;
	}	
	cout<<"INFO: #"<<regionFileNames.size()<<" region files in list..."<<endl;

	//##################################
	//##  PARSE SOURCE FILELIST
	//##################################
	std::vector<std::string> sourceFileNames;
	if(ReadList(sourceFileNames,sourceFileList)<0){
		cerr<<"ERROR: Failed to read list "<<sourceFileList<<endl;
		return -1;
	}	
	cout<<"INFO: #"<<sourceFileNames.size()<<" source files in list..."<<endl;

	//...
	//...

	return 0;

}//close macro


int ReadList(std::vector<std::string>& fileNames, std::string filelist)
{
	//Clear collection
	fileNames.clear();

	//Open file for reading
	ifstream in;  
  in.open(filelist.c_str());
  if(!in.good()) {
		cerr<<"ERROR: Failed to open config file "<<filelist<<" for reading!"<<endl;
		return -1;
  }

	//Parsing file
	cout<<"INFO: Reading and parsing file: "<<filelist<<endl;

	std::string parsedline= "";	
	int line_counter= 0;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;

		//Parse line fields
		//istringstream ss(parsedline);
		fileNames.push_back(parsedline);

		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();


	return 0;
	
}//close ReadList()

int ParseSourceData(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols<<" expected!"<<endl;
		return -1;
	}
	
	//Set var values
	info.name= std::string(fields[0]);
	info.x= atof(fields[3].c_str());//wcs x
	info.y= atof(fields[4].c_str());//wcs y
	info.flux= atof(fields[5].c_str());

	return 0;

}//close ParseData()


int ApplySelection(std::string sourceFileName,std::string regionFileName)
{
	

	return 0;

}//close ApplySelection()
