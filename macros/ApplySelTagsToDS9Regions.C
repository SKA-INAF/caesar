//==========================================
//==     METHODS
//==========================================
int ReadSelectionFile(std::string filename);
int ReadRegionFile(std::string filename,std::vector<std::string> skipLinePatterns,std::vector<std::string> keepLinePatterns);
bool HasPattern(std::string,std::string);
int ParseRegion(std::string& regionStr,const std::vector<std::string>& fields);

//==========================================
//==     VARS
//==========================================
struct SelectionTag {
	int tag;
	int edgeFlag;

	SelectionTag(int _tag,int _edgeFlag)
		: tag(_tag), edgeFlag(_edgeFlag)	
	{}

};
std::map<std::string,SelectionTag> gSelectionTagMap;
std::vector<std::string> gSkipLinePatterns {};
std::vector<std::string> gKeepLinePatterns {"image","global"};
std::vector<std::string> gParseLinePatterns {"ellipse","polygon"};
FILE* gOutFile= 0;
std::string gOutFileName= "ds9_augm.reg";

int ApplySelTagsToDS9Regions(std::string ds9file,std::string selectionfile)
{
	//Open output file
	gOutFile= fopen(gOutFileName.c_str(),"w");

	//Read selection file
	if(ReadSelectionFile(selectionfile)<0){
		cerr<<"ERROR: Failed to read selection file "<<selectionfile<<"!"<<endl;
		return -1;
	}

	//Read ds9 region
	if(ReadRegionFile(ds9file,gSkipLinePatterns,gKeepLinePatterns)<0){
		cerr<<"ERROR: Failed to read ds9 region file "<<ds9file<<"!"<<endl;
		return -1;
	}	

	//Close out file
	fclose(gOutFile);

	return 0;

}//close macro


int ReadSelectionFile(std::string filename)
{
	//Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		cerr<<"ERROR: Failed to open selection file "<<filename<<" for reading!"<<endl;
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
		bool keepLine= false;
		std::vector<std::string> fields;
		int field_counter= 0;

		while(ss >> field)
    {
			field_counter++;

			//Skip first line
			char first_char= *(field.c_str());
			if(field_counter==1){
				if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
					cout<<"skipped line detected (first_char="<<first_char<<", field="<<field<<", line="<<parsedline<<")..."<<endl;
					skipLine= true;
					break;
				}
			}

			//Collect fields
			fields.push_back(field);
			

		}//end loop line

		//Skip line?		
		if(skipLine) {
			continue;
		}
		
		//Skip line if no fields parsed
		if(fields.empty()) continue;
		if(fields.size()!=5){
			cerr<<"ERROR: Invalid number of fields parsed ("<<fields.size()<<") when 5 expected!"<<endl;
			continue;
		}

		//Process line
		long int sourceId= std::stol(fields[0]);
		int componentId= std::stoi(fields[1]);
		int nestedId= std::stoi(fields[2]);	
 		int edgeFlag= std::stoi(fields[3]);
		int tag= std::stoi(fields[4]);

		std::stringstream sourceNameStream;
		sourceNameStream<<"S"<<sourceId;
		if(nestedId!=0) sourceNameStream<<"_N"<<nestedId;
		sourceNameStream<<"_fitcomp"<<componentId;
	
		std::string sourceName= sourceNameStream.str();
		cout<<"INFO: sourceName="<<sourceName<<endl;

		//Fill map
		gSelectionTagMap.insert( std::make_pair(sourceName,SelectionTag(tag,edgeFlag)) );

		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	return 0;

}//close ReadSelectionFile()


int ReadRegionFile(std::string filename,std::vector<std::string> skipLinePatterns,std::vector<std::string> keepLinePatterns)
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
	std::stringstream keeplinesstream;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;

		//Parse line fields
		istringstream ss(parsedline);
		std::string field= "";
		bool skipLine= false;
		bool keepLine= false;
		std::vector<std::string> fields;
		int field_counter= 0;

		while(ss >> field)
    {
			field_counter++;

			//Skip first line
			char first_char= *(field.c_str());
			if(field_counter==1){
				if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
					cout<<"skipped line detected (first_char="<<first_char<<", field="<<field<<", line="<<parsedline<<")..."<<endl;
					skipLine= true;
					break;
				}
			}

			//Search for pattern and skip if present
			for(size_t i=0;i<skipLinePatterns.size();i++){
				bool hasPattern= HasPattern(field,skipLinePatterns[i]);
				if(hasPattern) {
					skipLine= true;
					break;
				}
			}
			if(skipLine) break;

			//Search for pattern to keep line
			for(size_t i=0;i<keepLinePatterns.size();i++){
				bool hasPattern= HasPattern(field,keepLinePatterns[i]);
				if(hasPattern) {
					keepLine= true;
					break;
				}
			}
			if(keepLine) break;

			//Collect fields
			fields.push_back(field);
			

		}//end loop line

		//Skip line?		
		if(skipLine) {
			cout<<"INFO: skipline: "<<parsedline<<endl;
			continue;
		}
		
		//Keep line	
		if(keepLine){	
			cout<<"INFO: keep line: "<<parsedline<<endl;
			keeplinesstream<<parsedline<<endl;
			fprintf(gOutFile,"%s\n",parsedline.c_str());
			continue;
		}		
	
		//Skip line if no fields parsed
		if(fields.empty()) continue;

		//Process line		
		std::string regionStr= "";
 		int status= ParseRegion(regionStr,fields);
		if(status==0){
			fprintf(gOutFile,"%s\n",regionStr.c_str());
		}
		else{
			cerr<<"ERROR: Failed to parse line "<<line_counter<<" of file "<<filename<<"!"<<endl;
		}	

		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	return 0;

}//close ReadRegionFile()

int ParseRegion(std::string& regionStr,const std::vector<std::string>& fields)
{
	//Init return str
	regionStr= "";

	//Check number of fields present
	cout<<"#"<<fields.size()<<" fields present in line..."<<endl;
	for(size_t i=0;i<fields.size();i++){
		cout<<"fields no. "<<i+1<<": "<<fields[i]<<endl;
	}

	/*
	//Parse region parameters until '#'
	std::string regionParStr= "";
	std::stringstream ss;
	int posIter= 0;
	for(size_t i=0;i<fields.size();i++){
		if(fields[i]=="#"){
			posIter= i+1;
			break;
		}
		ss<<fields[i]<<" ";
	}
	regionParStr= ss.str();
	cout<<"regionParStr: "<<regionParStr<<endl;
	*/

	//Find source name
	std::string sourceName= "";
	std::size_t pos_start;
	std::size_t pos_end;

	for(size_t i=0;i<fields.size();i++){
		pos_start = fields[i].find("text={S");	
		if(pos_start==std::string::npos) continue;
		
		pos_end= fields[i].find_last_of("}");
		if(pos_end==std::string::npos) continue;
		sourceName= fields[i].substr(pos_start+6,pos_end-pos_start-6);
		break;

	}//end loop field

	//Check if source name was properly parsed
	if(sourceName=="") return -1;

	//Extract source id
	pos_start = sourceName.find_first_of("S");	
	pos_end = sourceName.find_first_of("_");	
	std::string sourceId_str= sourceName.substr(pos_start+1,pos_end-pos_start-1);
	if(sourceId_str=="") return -1;
	long int sourceId= std::stol(sourceId_str);
	
	//Extract nested source id
	std::string nestedId_str= "";
	pos_start = sourceName.find("_N");
	if(pos_start!=std::string::npos){
		pos_end = sourceName.find("_fitcomp");
		nestedId_str= sourceName.substr(pos_start+2,pos_end-pos_start-2);
	}
	
	int nestedId= 0;
	if(nestedId_str!="") nestedId= std::stoi(nestedId_str);

	//Extract component id	
	pos_start = sourceName.find("_fitcomp");
	std::string componentId_str= sourceName.substr(pos_start+8,sourceName.length()-pos_start-8);
	if(componentId_str=="") return -1;
	int componentId= stoi(componentId_str);
	
	cout<<"sourceName="<<sourceName<<", sourceId="<<sourceId<<", nestedId="<<nestedId<<", componentId="<<componentId<<endl;

	//Retrieve tag and flags from map	
	std::string edgeTag_str= "";
	std::string sourceTag_str= "tag={not-classified}";
	auto it= gSelectionTagMap.find(sourceName);
	if(!gSelectionTagMap.empty() && it!=gSelectionTagMap.end()){//key not found
		SelectionTag st= it->second;
		int edgeFlag= st.edgeFlag;
		int tag= st.tag;
		
		if(edgeFlag==1) edgeTag_str= "tag={at-edge}";
		if(tag==1) sourceTag_str= "tag={real}";
		else if(tag==2) sourceTag_str= "tag={candidate}";
		else if(tag==3) sourceTag_str= "tag={fake}";

	}//close if

	//Append tags to fields
	std::vector<std::string> fields_augm= fields;
	if(edgeTag_str!="") fields_augm.push_back(edgeTag_str);
	if(sourceTag_str!="") fields_augm.push_back(sourceTag_str);
	
	//Populate output string
	for(size_t i=0;i<fields_augm.size()-1;i++){
		regionStr+= fields_augm[i] + std::string(" ");
	}
	regionStr+= fields_augm[fields_augm.size()-1];
	cout<<"INFO: regionStr="<<regionStr<<endl;

	return 0;

}//close ParseRegion()





bool HasPattern(std::string str,std::string pattern)
{
	std::size_t found = str.find(pattern);
  if (found!=std::string::npos) return true;

	return false;

}//close HasPattern()
