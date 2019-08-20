#include <Source.h>
#include <Logger.h>

#include <TObject.h>
#include <TEllipse.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>
using namespace std;
using namespace Caesar;

std::string TrimStr(const std::string& str,const std::string& whitespace = " \t")
{
	const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
		return ""; // no content

  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);

}//close TrimStr()

std::string ReduceStr(const std::string& str,
                   const std::string& fill = " ",
                   const std::string& whitespace = " \t")
{
	// trim first
  auto result = TrimStr(str, whitespace);

  // replace sub ranges
 	auto beginSpace = result.find_first_of(whitespace);
  while (beginSpace != std::string::npos){
  	const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
    const auto range = endSpace - beginSpace;
    result.replace(beginSpace, range, fill);

    const auto newStart = beginSpace + fill.length();
    beginSpace = result.find_first_of(whitespace, newStart);
  }

  return result;

}//close ReduceStr()


int ParseMGPSCatalogFile(std::string fileName,std::string outfileName="parsed.dat")
{
	//Check filename
	std::ifstream fileStream(fileName);	
	if (fileStream.fail() || !fileStream.is_open()){
		ERROR_LOG("Failed to open file "<<fileName<<" for reading...");
		return -1;
	}

	//Init output file
	FILE* fout= fopen(outfileName.c_str(),"w");

	//Read file line-by-line and parse
	std::vector<std::string> fileNames;
	std::string line;
	int nsources= 0;
	int lineCounter= 0;

	while (std::getline(fileStream, line)) {
  	
		//Check first character
		char first_char= *(line.c_str());
		if(first_char!='|'){//skip line
			continue;
		}
		std::istringstream iss(line);	
		lineCounter++;
		
		//Parse header (first line after |)
		if(lineCounter==1){
			INFO_LOG("Parsing header...");
			std::string field;
			int fieldId= 0;
			while (std::getline(iss, field, '|')) {
				fieldId++;
				if(field=="\n" || field=="" || field.empty()) continue;
  			cout<<"field="<<field<<endl;
			}	
			cout<<endl;

			continue;
		}//close if

		else{//parse source field
			cout<<"Source line: "<<line<<endl;

			//Parse fields
			std::string field;
			int fieldCounter= 0;
			while (std::getline(iss, field, '|')) {
				if(field=="\n" || field=="") continue;
				if(fieldCounter!=0) fprintf(fout,"  ");
				std::string field_parsed= ReduceStr(field);
  			cout<<"field="<<field_parsed<<endl;
				fprintf(fout,"%s",field_parsed.c_str());
				fieldCounter++;
			}	
			fprintf(fout,"\n");
			cout<<endl;
	
			nsources++;
		}//close else     

	}//end file read
				
	INFO_LOG("#"<<nsources<<" sources read...");

	//Close file
	fclose(fout);



	return 0;

}//close macro

