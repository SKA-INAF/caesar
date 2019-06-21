// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************
/**
* @file CutParser.cc
* @class CutParser
* @brief Parse the cut file containing cut parameters
* 
* @author S. Riggi
* @date 20/06/2019
*/

#include <CutParser.h>
#include <Cut.h>
#include <Consts.h>
#include <SysUtils.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

ClassImp(Caesar::CutParser)

namespace Caesar {

bool CutParser::m_hasRegisteredCuts;
std::string CutParser::m_cutFile= "";

CutParser::CutParser() 
{
	m_cutFile= "";
	m_hasRegisteredCuts= false;	

}//close costructor

CutParser::~CutParser()
{

}//close destructor

int CutParser::Parse(std::string filename)
{	
	//## Check and read file
	if(filename==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty filename string given!");
		#endif
		return -1;
	}

	Caesar::FileInfo info;
	if(!SysUtils::CheckFile(filename,info,false)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid config file specified (check if file actually exist)!");
		#endif
		return -1;
	}

	//## Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open config file "<<filename<<" for reading!");
		#endif
		return -1;
  }

	//## Parsing file
	m_cutFile= filename;
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading and parsing config file "<<filename<<" ...");
	#endif

	std::string parsedline= "";

	while(std::getline(in,parsedline)) 
	{
		//Check file
		if (!in.good()) break;

		//Check first character
		char first_char= *(parsedline.c_str());
		if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){//skip line
			continue;
		}
		
		istringstream line(parsedline);

		//Get cut name
		std::string cutName= "";    
		line >> cutName;
    if(cutName=="\n" || cutName=="") continue;
		#ifdef LOGGING_ENABLED
			INFO_LOG("Cut: "<<cutName);
		#endif

		//Get cut enabled flag
		std::string cutEnabledFlag= "";
		bool cutEnabled= false;
		line >> cutEnabledFlag;
		if(cutEnabledFlag=="T") cutEnabled= true;
		else if(cutEnabledFlag=="F") cutEnabled= false;
		else{
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to parse cut enabled flag (expected T/F, got "<<cutEnabledFlag<<"!");
			#endif
			return -1;
		}
	
		#ifdef LOGGING_ENABLED
			INFO_LOG("cutEnabledFlag: "<<cutEnabledFlag<<", cutEnabled? "<<cutEnabled);
		#endif

		//Get cut type
		std::string cutTypeStr= "";
		line >> cutTypeStr;
		bool cutReversed= false;
		int cutType= eUNKNOWN_CUT;

		if(cutTypeStr=="=="){//equality cut
			cutType= eEQUALITY_CUT;
			cutReversed= false;
		}
		else if(cutTypeStr=="!=") {//inequality cut
			cutType= eEQUALITY_CUT;
			cutReversed= true;
		}	
		else if(cutTypeStr=="<>"){//bound cut
			cutType= eBOUND_CUT;
			cutReversed= false;
		}
		else if(cutTypeStr=="><"){//bound cut reversed
			cutType= eBOUND_CUT;
			cutReversed= true;
		}
		else if(cutTypeStr==">="){//single bound cut
			cutType= eSINGLE_BOUND_CUT;
			cutReversed= false;
		}
		else if(cutTypeStr=="<="){//single bound cut reversed
			cutType= eSINGLE_BOUND_CUT;
			cutReversed= true;
		}
		else{
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Unknown cut type found ("<<cutTypeStr<<")!");
			#endif
			return -1;
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("cutTypeStr: "<<cutTypeStr<<", cutType="<<cutType<<", cutReversed? "<<cutReversed);
		#endif

		//Get cut AND/OR
		std::string cutCombineFlag= "";
		line >> cutCombineFlag;
		bool cutCombineInOR= false;
		if(cutCombineFlag=="OR") cutCombineInOR= true;
		else if(cutCombineFlag=="AND") cutCombineInOR= false;
		else{
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Unknown/invalid cut combine flag found ("<<cutCombineFlag<<")!");
			#endif
			return -1;
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("cutCombineFlag: "<<cutCombineFlag<<", cutCombineInOR? "<<cutCombineInOR);
		#endif

		//Register cut parameters
		if(cutType==eEQUALITY_CUT)
		{
			std::string cutParamsField= "";
			line >> cutParamsField;
			if(cutParamsField=="\n" || cutParamsField=="") {
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Missing cut parameters in equality cut parsed!");
				#endif
				return -1;
			}
			CodeUtils::RemovePatternInString(cutParamsField,"{");
			CodeUtils::RemovePatternInString(cutParamsField,"}");
			std::vector<std::string> cutParams_str= CodeUtils::SplitStringOnPattern(cutParamsField,',');
			std::vector<double> cutParams= CodeUtils::StringVecToTypedVec<double>(cutParams_str);	
			CutFactory::Instance().RegisterEqualityCut(cutName,cutParams,cutEnabled,cutReversed,cutCombineInOR);
		}
		else if(cutType==eBOUND_CUT)
		{
			std::string cutParamsField_min= "";
			std::string cutParamsField_max= "";
			line >> cutParamsField_min;
			line >> cutParamsField_max;
			if(cutParamsField_min=="\n" || cutParamsField_min=="") {
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Missing min cut parameters in bound cut parsed!");
				#endif
				return -1;
			}	
			if(cutParamsField_max=="\n" || cutParamsField_max=="") {
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Missing max cut parameters in bound cut parsed!");
				#endif
				return -1;
			}
			CodeUtils::RemovePatternInString(cutParamsField_min,"{");
			CodeUtils::RemovePatternInString(cutParamsField_min,"}");
			CodeUtils::RemovePatternInString(cutParamsField_max,"{");
			CodeUtils::RemovePatternInString(cutParamsField_max,"}");
			std::vector<std::string> cutMinParams_str= CodeUtils::SplitStringOnPattern(cutParamsField_min,',');
			std::vector<std::string> cutMaxParams_str= CodeUtils::SplitStringOnPattern(cutParamsField_max,',');
			std::vector<double> cutMinParams= CodeUtils::StringVecToTypedVec<double>(cutMinParams_str);	
			std::vector<double> cutMaxParams= CodeUtils::StringVecToTypedVec<double>(cutMaxParams_str);
			if(cutMaxParams.size()!=cutMinParams.size()){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Parsed min/max bound parameter vectors have different size!");
				#endif
				return -1;
			}
			std::vector<std::pair<double,double>>	cutParams;
			for(size_t i=0;i<cutMaxParams.size();i++) cutParams.push_back( std::make_pair(cutMinParams[i],cutMaxParams[i]) );
			CutFactory::Instance().RegisterBoundCut(cutName,cutParams,cutEnabled,cutReversed,cutCombineInOR);			
		}
		else if(cutType==eSINGLE_BOUND_CUT){
			std::string cutParamsField= "";
			line >> cutParamsField;
			if(cutParamsField=="\n" || cutParamsField=="") {
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Missing cut parameters in bound cut parsed!");
				#endif
				return -1;
			}	
			CodeUtils::RemovePatternInString(cutParamsField,"{");
			CodeUtils::RemovePatternInString(cutParamsField,"}");
			std::vector<std::string> cutParams_str= CodeUtils::SplitStringOnPattern(cutParamsField,',');
			std::vector<double> cutParams= CodeUtils::StringVecToTypedVec<double>(cutParams_str);	
			CutFactory::Instance().RegisterSingleBoundCut(cutName,cutParams,cutEnabled,cutReversed,cutCombineInOR);			
		}
		else{
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Unknown/invalid cut type found ("<<cutType<<")!");
			#endif
			return -1;
		}
		
		if (!in.good()) break;
		
	}//close while

	//## Print parsed cuts
	#ifdef LOGGING_ENABLED
		INFO_LOG("Printing parsed cuts...");
	#endif
	PrintCuts();

	//## Close file
	in.close();

	return 0;

}//close Parse()


int CutParser::RegisterPredefinedCuts()
{	
	//Register pre-defined cuts
	try {
		//=============
		//==  CUTS   ==
		//=============
		//REGISTER_OPTION(inputFile,std::string,"","","");
		
		
		//Set has_registered flag (otherwise cuts are re-built)
		m_hasRegisteredCuts= true;

	}//close try block
	catch(...)
	{
		cerr<<"CutParser::RegisterPredefinedCuts(): ERROR: Failed to load predefined cuts!"<<endl;
		return -1;
	}

	return 0;

}//close RegisterPredefinedCuts()

}//close namespace

