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
* @file Option.h
* @class Option
* @brief Define configuration option
* 
* @author S. Riggi
* @date 5/6/2016
*/

#ifndef _OPTION_h
#define _OPTION_h 1

#include <Logger.h>

#include <json/json.h>

#include <TObject.h>
#include <TBranch.h>
#include <TTree.h>

#include <vector>
#include <string>

#include <string>

using namespace std;

namespace Caesar{

class OptionBase : public TObject {
	public:
    virtual ~OptionBase() {}
	public:
		virtual void Print() = 0;
		virtual int SetValueFromStream(std::stringstream& sstream) = 0;
		virtual int SetValueFromString(std::string const& s) = 0;
		virtual int AddBranch(TTree* tree) = 0;
		virtual int GetJson(Json::Value&) = 0;
		virtual int GetJsonString(std::string&,bool) = 0;
		template <typename T>
    	int GetValue(T& t);
		template <typename T>
    	int GetDefaultValue(T& t);
		template <typename T>
    	int GetRange(T& xmin,T& xmax);

	ClassDef(OptionBase,1)	
};
typedef OptionBase* OptionBasePtr;
typedef std::map<std::string, OptionBasePtr> OptionMap;



//== OPTION HELPERS ==
template <typename T>
struct OptionHelper  {
	
	public:
		OptionHelper(){};
		virtual ~OptionHelper(){};
	
	public:
	static bool CheckRange(T& value,T& minValue,T& maxValue){
		if(typeid(T) == typeid(std::string)){
			return true;
		}
		if( (typeid(T)==typeid(char)) || (typeid(T)==typeid(char*)) ) {
			return true;
		}
		return (value>minValue && value<maxValue);
	}//close CheckRange()

	static int SetValueFromStream(std::stringstream& sstream,T& value,T& minValue,T& maxValue) {
		T tmpValue= T(0);	
		sstream>>tmpValue;	
		if(!CheckRange(tmpValue,minValue,maxValue)) return -1;
		value= tmpValue;
		return 0;
	}
	static int GetJson(Json::Value& root,std::string& name,T& value,T& defaultValue,T& minValue,T& maxValue) {
		root["name"]= name;
		root["value"]= value;	
		root["default"]= defaultValue;
		root["min"]= minValue;
		root["max"]= maxValue;
		return 0;
	}
	static std::string GetPrintable(std::string& name,T& value,T& defaultValue,T& minValue,T& maxValue){
		std::stringstream ss;
		ss<<"Option: {\"name\": \""<<name<<"\", \"value\": "<<value<<", \"default\": "<<defaultValue<<", \"min\": "<<minValue<<" \"max\": "<<maxValue<<"}";
		return ss.str();
	}
	
	ClassDef(OptionHelper,1)

};//close OptionHelper()


template <>
struct OptionHelper<std::string> {

	static bool CheckRange(std::string& value,std::string& minValue,std::string& maxValue){
		return true;
	}
	static int SetValueFromStream(std::stringstream& sstream,std::string& value,std::string& minValue,std::string& maxValue) {
		std::string tmpValue= "";	
		sstream>>tmpValue;	
		if(!CheckRange(tmpValue,minValue,maxValue)) return -1;
		value= tmpValue;
		return 0;
	}
	static int GetJson(Json::Value& root,std::string& name,std::string& value,std::string& defaultValue,std::string& minValue,std::string& maxValue) {
		root["name"]= name;
		root["value"]= value;	
		root["default"]= defaultValue;
		root["min"]= minValue;
		root["max"]= maxValue;
		return 0;
	}
	static std::string GetPrintable(std::string& name,std::string& value,std::string& defaultValue,std::string& minValue,std::string& maxValue){
		std::stringstream ss;
		ss<<"Option: {\"name\": \""<<name<<"\", \"value\": "<<value<<", \"default\": "<<defaultValue<<", \"min\": "<<minValue<<" \"max\": "<<maxValue<<"}";
		return ss.str();
	}

};//close template specialization for string

template <>
struct OptionHelper<long int> {
	
	static bool CheckRange(long int& value,long int& minValue,long int& maxValue){
		return (value>minValue && value<maxValue);
	}
	static int SetValueFromStream(std::stringstream& sstream,long int& value,long int& minValue,long int& maxValue) {
		long int tmpValue;	
		sstream>>tmpValue;	
		if(!CheckRange(tmpValue,minValue,maxValue)) return -1;
		value= tmpValue;
		return 0;
	}
	static int GetJson(Json::Value& root,std::string& name,long int& value,long int& defaultValue,long int& minValue,long int& maxValue) {
		root["name"]= name;
		root["value"]= Json::Value( static_cast<Json::Int64>(value) );
		root["default"]= Json::Value( static_cast<Json::Int64>(defaultValue) );
		root["min"]= Json::Value( static_cast<Json::Int64>(minValue) );
		root["max"]= Json::Value( static_cast<Json::Int64>(maxValue) );
		return 0;
	}	
	static std::string GetPrintable(std::string& name,long int& value,long int& defaultValue,long int& minValue,long int& maxValue){
		std::stringstream ss;
		ss<<"Option: {\"name\": \""<<name<<"\", \"value\": "<<value<<", \"default\": "<<defaultValue<<", \"min\": "<<minValue<<" \"max\": "<<maxValue<<"}";
		return ss.str();
	}
	//ClassDef(OptionHelper<long int>,1)

};//close template specialization for long int


template <>
struct OptionHelper<bool>  {
	static bool CheckRange(bool& value,bool& minValue,bool& maxValue){
		return true;
	}
	static int SetValueFromStream(std::stringstream& sstream,bool& value,bool& minValue,bool& maxValue) {
		sstream>> std::boolalpha >> value;
		return 0;
	}
	static int GetJson(Json::Value& root,std::string& name,bool& value,bool& defaultValue,bool& minValue,bool& maxValue) {
		root["name"]= name;
		root["value"]= value;	
		root["default"]= defaultValue;
		root["min"]= minValue;
		root["max"]= maxValue;
		return 0;
	}
	static std::string GetPrintable(std::string& name,bool& value,bool& defaultValue,bool& minValue,bool& maxValue){
		std::stringstream ss;
		ss<<std::boolalpha<<"Option: {\"name\": \""<<name<<"\", \"value\": "<<value<<", \"default\": "<<defaultValue<<", \"min\": "<<minValue<<" \"max\": "<<maxValue<<"}";
		return ss.str();
	}
	//ClassDef(OptionHelper<bool>,1)
};//close template specialization for long int


template <typename T>
struct OptionHelper< std::vector<T> > {
	static bool CheckRange(std::vector<T>& value,std::vector<T>& minValue,std::vector<T>& maxValue){	
		if(typeid(T) == typeid(std::string)){
			return true;
		}
		if( (typeid(T)==typeid(char)) || (typeid(T)==typeid(char*)) ) {
			return true;
		}
		for(unsigned int i=0;i<value.size();i++){
			if(value[i]<minValue[i] || value[i]>maxValue[i]) return false;
		}
		return true;
	}

	static int SetValueFromStream(std::stringstream& sstream,std::vector<T>& value,std::vector<T>& minValue,std::vector<T>& maxValue) {
		std::vector<T> tmpValue;
		do {
    	for (T thisVal; sstream >> thisVal;) {
      	tmpValue.push_back(thisVal);
      }
      if (sstream.fail()) {
      	sstream.clear();
       	std::string token;
        sstream >> token;
      }
    }
    while (!sstream.eof());
		if(!CheckRange(tmpValue,minValue,maxValue)) return -1;
		value.insert(value.end(),tmpValue.begin(),tmpValue.end());
		return 0;
	}
	
	static int GetJson(Json::Value& root,std::string& name,std::vector<T>& value,std::vector<T>& defaultValue,std::vector<T>& minValue,std::vector<T>& maxValue) {
		for(unsigned int i=0;i<value.size();i++){	
			root["value"].append( Json::Value(value[i]) );
		}
		for(unsigned int i=0;i<defaultValue.size();i++){	
			root["default"].append( Json::Value(defaultValue[i]) );
		}
		for(unsigned int i=0;i<minValue.size();i++){	
			root["min"].append( Json::Value(minValue[i]) );
		}
		for(unsigned int i=0;i<maxValue.size();i++){	
			root["max"].append( Json::Value(maxValue[i]) );
		}
		return 0;
	}
	static std::string GetPrintable(std::string& name,std::vector<T>& value,std::vector<T>& defaultValue,std::vector<T>& minValue,std::vector<T>& maxValue){
		std::stringstream ss;
		ss<<"Option: {\"name\": \""<<name<<"\", \"value\": [";
		for(unsigned int i=0;i<value.size()-1;i++) ss<<value[i]<<",";
		ss<<value[value.size()-1]<<"]";

		ss<<", \"default\": [";
		for(unsigned int i=0;i<defaultValue.size()-1;i++) ss<<defaultValue[i]<<",";
		ss<<defaultValue[defaultValue.size()-1]<<"]";

		ss<<", \"min\": [";
		for(unsigned int i=0;i<minValue.size()-1;i++) ss<<minValue[i]<<",";
		ss<<minValue[minValue.size()-1]<<"]";

		ss<<" \"max\": [";
		for(unsigned int i=0;i<maxValue.size()-1;i++) ss<<maxValue[i]<<",";
		ss<<maxValue[maxValue.size()-1]<<"]";
		ss<<"}";
		return ss.str();
	}
	//ClassDef(OptionHelper<std::vector<T>>,1)
};//close template specialization for std::vector<T>


template <>
struct OptionHelper< std::vector<long int> >  {
	static bool CheckRange(std::vector<long int>& value,std::vector<long int>& minValue,std::vector<long int>& maxValue){	
		for(unsigned int i=0;i<value.size();i++){
			if(value[i]<minValue[i] || value[i]>maxValue[i]) return false;
		}
		return true;
	}
	static int SetValueFromStream(std::stringstream& sstream,std::vector<long int>& value,std::vector<long int>& minValue,std::vector<long int>& maxValue) {
		std::vector<long int> tmpValue;
		do {
    	for (long int thisVal; sstream >> thisVal;) {
      	tmpValue.push_back(thisVal);
      }
      if (sstream.fail()) {
      	sstream.clear();
       	std::string token;
        sstream >> token;
      }
    }
    while (!sstream.eof());
		if(!CheckRange(tmpValue,minValue,maxValue)) return -1;
		value.insert(value.end(),tmpValue.begin(),tmpValue.end());
		return 0;
	}

	static int GetJson(Json::Value& root,std::string& name,std::vector<long int>& value,std::vector<long int>& defaultValue,std::vector<long int>& minValue,std::vector<long int>& maxValue) {
		for(unsigned int i=0;i<value.size();i++){	
			root["value"].append( Json::Value( static_cast<Json::Int64>(value[i]) ) );
		}
		for(unsigned int i=0;i<defaultValue.size();i++){	
			root["default"].append( Json::Value( static_cast<Json::Int64>(defaultValue[i]) ) );
		}
		for(unsigned int i=0;i<minValue.size();i++){	
			root["min"].append( Json::Value( static_cast<Json::Int64>(minValue[i]) ) );
		}
		for(unsigned int i=0;i<maxValue.size();i++){	
			root["max"].append( Json::Value( static_cast<Json::Int64>(maxValue[i]) ) );
		}
		return 0;
	}
	static std::string GetPrintable(std::string& name,std::vector<long int>& value,std::vector<long int>& defaultValue,std::vector<long int>& minValue,std::vector<long int>& maxValue){
		std::stringstream ss;
		ss<<"Option: {\"name\": \""<<name<<"\", \"value\": [";
		for(unsigned int i=0;i<value.size()-1;i++) ss<<value[i]<<",";
		ss<<value[value.size()-1]<<"]";

		ss<<", \"default\": [";
		for(unsigned int i=0;i<defaultValue.size()-1;i++) ss<<defaultValue[i]<<",";
		ss<<defaultValue[defaultValue.size()-1]<<"]";

		ss<<", \"min\": [";
		for(unsigned int i=0;i<minValue.size()-1;i++) ss<<minValue[i]<<",";
		ss<<minValue[minValue.size()-1]<<"]";

		ss<<" \"max\": [";
		for(unsigned int i=0;i<maxValue.size()-1;i++) ss<<maxValue[i]<<",";
		ss<<maxValue[maxValue.size()-1]<<"]";
		ss<<"}";
		return ss.str();
	}
	//ClassDef(OptionHelper<std::vector<long int>>,1)
};//close template specialization for std::vector<long int>


template <>
struct OptionHelper< std::vector<bool> > {
	
	static int SetValueFromStream(std::stringstream& sstream,std::vector<bool>& value,std::vector<bool>& minValue,std::vector<bool>& maxValue) {
		std::vector<bool> tmpValue;
		do {
    	for (bool thisVal; sstream >> std::boolalpha >> thisVal;) {
      	tmpValue.push_back(thisVal);
      }
      if (sstream.fail()) {
      	sstream.clear();
       	std::string token;
        sstream >> token;
      }
    }
    while (!sstream.eof());
		value.insert(value.end(),tmpValue.begin(),tmpValue.end());
		return 0;
	}
	static int GetJson(Json::Value& root,std::string& name,std::vector<bool>& value,std::vector<bool>& defaultValue,std::vector<bool>& minValue,std::vector<bool>& maxValue) {
		for(unsigned int i=0;i<value.size();i++){	
			root["value"].append( Json::Value(value[i]) );
		}
		for(unsigned int i=0;i<defaultValue.size();i++){	
			root["default"].append( Json::Value(defaultValue[i]) );
		}
		for(unsigned int i=0;i<minValue.size();i++){	
			root["min"].append( Json::Value( minValue[i]) );
		}
		for(unsigned int i=0;i<maxValue.size();i++){	
			root["max"].append( Json::Value( maxValue[i]) );
		}
		return 0;
	}
	static std::string GetPrintable(std::string& name,std::vector<bool>& value,std::vector<bool>& defaultValue,std::vector<bool>& minValue,std::vector<bool>& maxValue){
		std::stringstream ss;
		ss<<"Option: {\"name\": \""<<name<<"\", \"value\": [";
		for(unsigned int i=0;i<value.size()-1;i++) ss<<value[i]<<",";
		ss<<value[value.size()-1]<<"]";

		ss<<", \"default\": [";
		for(unsigned int i=0;i<defaultValue.size()-1;i++) ss<<defaultValue[i]<<",";
		ss<<defaultValue[defaultValue.size()-1]<<"]";

		ss<<", \"min\": [";
		for(unsigned int i=0;i<minValue.size()-1;i++) ss<<minValue[i]<<",";
		ss<<minValue[minValue.size()-1]<<"]";

		ss<<" \"max\": [";
		for(unsigned int i=0;i<maxValue.size()-1;i++) ss<<maxValue[i]<<",";
		ss<<maxValue[maxValue.size()-1]<<"]";
		ss<<"}";
		return ss.str();
	}
	//ClassDef(OptionHelper<std::vector<bool>>,1)
};//close template specialization for std::vector<bool>
//=====================================================================





template <typename T> 
class Option : public OptionBase {

	public:
		/** 
		\brief Default Class constructor (needed by ROOT IO)
 		*/
		Option(){};
		/** 
		\brief Class constructor
 		*/
  	Option(std::string name,T defaultValue,T minValue,T maxValue) : 
			m_name(name), m_defaultValue(defaultValue), m_minValue(minValue), m_maxValue(maxValue)
		{
			m_value= m_defaultValue;
			m_hasSetValue= false;
		};
		/** 
		\brief Class destructor
 		*/
  	virtual ~Option(){};

		typedef T OptionType;

	public:
		/** 
		\brief Set option value from stringstream
 		*/
		int SetValueFromStream(std::stringstream& sstream){
			//Set value (range check performed internally)
			if(OptionHelper<T>::SetValueFromStream(sstream,m_value,m_minValue,m_maxValue)<0){
				return -1;
			}	
			m_hasSetValue= true;
			return 0;	
		}

		/*
		int SetValueFromStream(std::stringstream& sstream){
			T tmpValue;
			if(typeid(T) == typeid(bool)){
				sstream >> std::boolalpha >> tmpValue;	
			}
			else{
				sstream>>tmpValue;
			}
			
			if(CheckRange() && (tmpValue<m_minValue || tmpValue>m_maxValue) ) {
				return -1;
			}
			m_value= tmpValue;
			m_hasSetValue= true;
			return 0;
		}
		*/

		/** 
		\brief Set option value from string
 		*/
		int SetValueFromString(std::string const& s){
			std::stringstream sstream(s);
			return SetValueFromStream(sstream);
		}
		/** 
		\brief Set option value
 		*/
		int SetValue(T value){
			if(!OptionHelper<T>::CheckRange(value,m_minValue,m_maxValue)) return -1;
			m_value= value;
			m_hasSetValue= true;
			return 0;
		}
		/** 
		\brief Get option value
 		*/
		T GetValue(){	
			if(!m_hasSetValue) return GetDefaultValue();//return default if no value has been set
			return m_value;
		}
		/** 
		\brief Get default value
 		*/	
		T GetDefaultValue(){
			return m_defaultValue;
		}

		/** 
		\brief Get json object from option
 		*/
		int GetJson(Json::Value& jsonObj){
			return OptionHelper<T>::GetJson(jsonObj,m_name,m_value,m_defaultValue,m_minValue,m_maxValue);	
		}
		
		/** 
		\brief Get json string
 		*/
		int GetJsonString(std::string& jsonString,bool isMinified=true){
			Json::Value jsonObj;
			if(GetJson(jsonObj)<0) {
				ERROR_LOG("Failed to encode option to json object!");
				return -1;
			}
			try {
				if(isMinified){// write in a minified way
					Json::FastWriter fastWriter;
					jsonString= fastWriter.write(jsonObj);
				}
				else{	// write in a nice readible way
					Json::StyledWriter formattedWriter;
					jsonString= formattedWriter.write(jsonObj);
				}
			}
			catch(...){
				ERROR_LOG("Failed to encode argument to json string!");
				return -1;
			}
			return 0;
		};

		/** 
		\brief Check given value against registered range (non sense for bool, string/char)?
 		*/	
		bool CheckRange(T& value){
			return OptionHelper<T>::CheckRange(value,m_minValue,m_maxValue);	
		}
		bool CheckRange(T* value){
			return OptionHelper<T>::CheckRange(*value,m_minValue,m_maxValue);	
		}

		/*
		bool CheckRange(){
			if(typeid(T) == typeid(std::string)){
				//return false;
				return true;
			}
			if( (typeid(T)==typeid(char)) || (typeid(T)==typeid(char*)) ) {
				//return false;	
				return true;
			}
			if(typeid(T) == typeid(bool)){
				//return false;
				return true;
			}
			return (m_value>m_minValue && m_value<m_maxValue);
			//return true;
		}//close CheckRange()
		*/

		/** 
		\brief Print option
 		*/
		void Print(){
			cout<<GetPrintable()<<endl;
			//cout<<"Option: {\"name\": \""<<m_name<<"\", \"value\": "<<m_value<<", \"default\": "<<m_defaultValue<<", \"min\": "<<m_minValue<<" \"max\": "<<m_maxValue<<"}"<<endl;
		}

		std::string GetPrintable(){
			return OptionHelper<T>::GetPrintable(m_name,m_value,m_defaultValue,m_minValue,m_maxValue);
		}

		/** 
		\brief Add a TBranch to input TTree
 		*/
		int AddBranch(TTree* tree){
			if(!tree) {
				ERROR_LOG("Null tree ptr given!");
				return -1;
			}
			TBranch* branch= tree->Branch(m_name.c_str(),this);
			if(!branch) {
				WARN_LOG("Failed to create a branch for this option...");
				return -1;
			}
			return 0;
		}

	protected:
		
		std::string m_name;
		T m_value;
		T m_defaultValue;
		T m_minValue;
		T m_maxValue;
		bool m_hasSetValue;

	ClassDef(Option,1)	
};


//== OPTION SPECIALIZATION ==
//*** GetJson()  ********
/*
template <typename T> inline
int Option<T>::GetJson(Json::Value& jsonObj){	
	jsonObj["name"]= m_name;
	jsonObj["value"]= m_value;	
	jsonObj["default"]= m_defaultValue;
	jsonObj["min"]= m_minValue;
	jsonObj["max"]= m_maxValue;
	return 0;
}
template <> inline
int Option<long int>::GetJson(Json::Value& jsonObj) {
 	jsonObj["name"]= m_name;
	jsonObj["value"]= Json::Value( static_cast<Json::Int64>(m_value) );
	jsonObj["default"]= Json::Value( static_cast<Json::Int64>(m_defaultValue) );
	jsonObj["min"]= Json::Value( static_cast<Json::Int64>(m_minValue) );
	jsonObj["max"]= Json::Value( static_cast<Json::Int64>(m_maxValue) );
	return 0;
}
*/
//****************************

//*** SetValueFromStream()   *****
/*
template <typename T> inline
int Option<T>::SetValueFromStream(std::stringstream& sstream) {
	T tmpValue;
	sstream>>tmpValue;
	if( !CheckRange(tmpValue) ) return -1;
	m_value= tmpValue;
	m_hasSetValue= true;
	return 0;
}
*/




template <typename T>
int OptionBase::GetValue(T& val) {
	try{
		auto castThis = dynamic_cast<Option<T>*>(this);
		val= castThis->m_value;
	}
	catch(std::bad_cast& e){
		cerr<<"Option::GetValue(): ERROR: Failed to get option value (invalid type cast)!"<<endl;
		return -1;
	}
	return 0;		
}//close GetValue()

template <typename T>
int OptionBase::GetDefaultValue(T& val) {
	try{
		auto castThis = dynamic_cast<Option<T>*>(this);
		val= castThis->m_defaultValue;
	}
	catch(std::bad_cast& e){
		cerr<<"Option::GetDefaultValue(): ERROR: Failed to get default option value (invalid type cast)!"<<endl;
		return -1;
	}
	return 0;		
}

template <typename T>
int OptionBase::GetRange(T& minVal,T& maxVal) {
	try{
		auto castThis = dynamic_cast<Option<T>*>(this);
		minVal= castThis->m_minValue;
		maxVal= castThis->m_maxValue;
	}
	catch(std::bad_cast& e){
		cerr<<"Option::GetRange(): ERROR: Failed to get default option value (invalid type cast)!"<<endl;
		return -1;
	}
	return 0;		
}






//== OPTION FACTORY ==
class OptionFactory : public TObject {

	public:
		static OptionFactory& Instance() {
    	// Since it's a static variable, if the class has already been created,
      // It won't be created again.
      // And it is thread-safe in C++11.
      static OptionFactory myInstance;
 
      // Return a reference to our instance.
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    OptionFactory(OptionFactory const&) = delete;             // Copy construct
    OptionFactory(OptionFactory&&) = delete;                  // Move construct
    OptionFactory& operator=(OptionFactory const&) = delete;  // Copy assign
    OptionFactory& operator=(OptionFactory &&) = delete;      // Move assign


	public:

		/** 
		\brief Create an option
 		*/
		template<typename T>
		static OptionBasePtr Create(std::string name,T defaultValue,T minValue,T maxValue) {	
			return OptionBasePtr(new Option<T>(name,defaultValue,minValue,maxValue));
		}//close Create()

		/** 
		\brief Register option
 		*/
		template <typename T>
  		int Register(std::string name,T defaultValue,T minValue,T maxValue) {
				//Check args
				if(name=="") {
					cerr<<"OptionFactory()::Register(): WARN: Invalid option name given!"<<endl;					
					return -1;
				}
				if(minValue>maxValue) {
					cerr<<"OptionFactory()::Register(): WARN: Invalid range parameter given!"<<endl;					
					return -1;
				}

				//Check if option with given name already exists
				OptionBasePtr thisOption= GetOptionBase(name);
				if(thisOption){
					cerr<<"OptionFactory()::Register(): WARN: Option "<<name<<" already registered, skip it!"<<endl;
					return 0;
				}
	
				//Option does not exist, create one and add to the map
				OptionBasePtr aNewOption= Create<T>(name,defaultValue,minValue,maxValue);
				m_RegisteredOptions.insert( std::pair<std::string,OptionBasePtr>(name,aNewOption) );
		
				return 0;
  		}//close RegisterOption()
		
		/** 
		\brief Has option?
 		*/
		bool HasOption(std::string name){
			OptionMap::iterator it= m_RegisteredOptions.find(name);
			if ( m_RegisteredOptions.empty() || it==m_RegisteredOptions.end() )
				return false;
			return true;
		}

		/** 
		\brief Get option base pointer
 		*/
		OptionBasePtr GetOptionBase(std::string name){
			OptionMap::iterator it= m_RegisteredOptions.find(name);
			if ( m_RegisteredOptions.empty() || it==m_RegisteredOptions.end() ) return nullptr;
			return it->second;
		}

		/** 
		\brief Get option impl pointer
 		*/
		template <typename T>
		//std::shared_ptr<Option<T>> GetOption(std::string name){
		Option<T>* GetOption(std::string name){
			//Try to dynamical cast to template type
			OptionBasePtr thisOption= GetOptionBase(name);
			if(!thisOption){
				cerr<<"OptionFactory()::GetOption(): INFO: Option not found!"<<endl;
				return 0;
			}
			//std::shared_ptr<Option<T>> thisCastedOption= std::dynamic_pointer_cast<Option<T>>(thisOption);
			Option<T>* thisCastedOption= dynamic_cast<Option<T>*>(thisOption);
			
			if(!thisCastedOption){
				cerr<<"OptionFactory()::GetOption(): ERROR: Failed to cast option to given data type (option registered with different type?)!"<<endl;
				return 0;
			}
			return thisCastedOption;
		}

		/** 
		\brief Get option value
 		*/
		template <typename T>
			int GetOptionValue(std::string name,T& value){
				Option<T>* thisOption= GetOption<T>(name);
				if(!thisOption) return -1;
				value= thisOption->GetValue();
				return 0;
			}//close GetOptionValue()

		/** 
		\brief Set option value
 		*/
		template <typename T>
			int SetOption(std::string name,T value){
				//std::shared_ptr<Option<T>> thisOption= GetOption<T>(name);	
				Option<T>* thisOption= GetOption<T>(name);
				if(!thisOption){
					cerr<<"OptionFactory()::SetOption(): ERROR: Cannot set value ("<<value<<") in option "<<name<<" (option not registered)!"<<endl;	
					return -1;
				}
				thisOption->SetValue(value);
				return 0;
			}//close SetOption()

		/** 
		\brief Set option value from string
 		*/
		int SetOptionFromString(std::string name,std::string stringified_value){
				OptionBasePtr thisOption= GetOptionBase(name);
				if(!thisOption){
					cerr<<"OptionFactory()::SetOptionFromString(): ERROR: Cannot set value ("<<stringified_value<<") in option "<<name<<" (option not registered)!"<<endl;	
					return -1;
				}
				thisOption->SetValueFromString(stringified_value);
				return 0;
			}//close SetOption()
		

		/** 
		\brief Print registered options
 		*/
		void PrintOptions(){
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
    		(it->second)->Print();
			}
		}//close PrintOptions()

		/** 
		\brief Get json
 		*/
		int GetJson(Json::Value& root){
			Json::Value optionArray;
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
				Json::Value jsonObj;
    		int status= (it->second)->GetJson(jsonObj);
				if(status<0) return -1;
				optionArray.append(jsonObj);
			}
			root["options"]= optionArray;
			return 0;
		}

		/** 
		\brief Get json string
 		*/
		int GetJsonString(std::string& jsonString,bool isMinified=true){
			Json::Value jsonObj;
			if(GetJson(jsonObj)<0) {
				cerr<<"OptionFactory::GetJsonString(): ERROR: Failed to encode option to json object!"<<endl;
				return -1;
			}
			try {
				if(isMinified){// write in a minified way
					Json::FastWriter fastWriter;
					jsonString= fastWriter.write(jsonObj);
				}
				else{	// write in a nice readible way
					Json::StyledWriter formattedWriter;
					jsonString= formattedWriter.write(jsonObj);
				}
			}
			catch(...){
				cerr<<"OptionFactory::GetJsonString(): ERROR: Failed to encode argument to json string!"<<endl;
				return -1;
			}
			return 0;
		};

		/** 
		\brief Make a TTree with all options
 		*/	
		TTree* MakeTTree(std::string treeName="CfgInfo"){

			//Create tree
			TTree* tree= new TTree(treeName.c_str(),treeName.c_str());

			//Adding branches
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
    		if( (it->second)->AddBranch(tree)<0 ) {
					cerr<<"OptionFactory::MakeTTree(): WARN: Failed to add branch for current option... "<<endl;
					continue;
				}
			}

			//Fill tree
			tree->Fill();

			return tree;
		}//close MakeTTree()
		
	protected:
		OptionFactory(){};
		virtual ~OptionFactory(){
			//Deleting map with options
			cout<<"~OptionFactory(): INFO: Deleting registered options..."<<endl;
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
				if(it->second) delete it->second;
			}//end loop map
			m_RegisteredOptions.clear();
		};

	private:
		OptionMap m_RegisteredOptions;		
	
	ClassDef(OptionFactory,1)	
};



#define REGISTER_OPTION(name,type,default_value,min_value,max_value) (OptionFactory::Instance().Register<type>(#name,default_value,min_value,max_value))
#define GET_OPTION(name,type) (OptionFactory::Instance().GetOption<type>(#name))
#define SET_OPTION(name,type,value) (OptionFactory::Instance().SetOption<type>(#name,value))
#define PRINT_OPTIONS() (OptionFactory::Instance().PrintOptions())
#define GET_OPTION_VALUE(name,value) (OptionFactory::Instance().GetOptionValue(#name,value))


#ifdef __MAKECINT__
#pragma link C++ class OptionBase+;
#pragma link C++ class OptionBase*+;
#pragma link C++ class std::map<std::string,Caesar::OptionBase*>+;
#pragma link C++ class Option<int>+;
#pragma link C++ class Option<long int>+;
#pragma link C++ class Option<float>+;
#pragma link C++ class Option<double>+;
#pragma link C++ class Option<std::string>+;
#pragma link C++ class Option<bool>+;
#pragma link C++ class OptionFactory+;
#pragma link C++ MACRO REGISTER_OPTION;
#pragma link C++ MACRO GET_OPTION;
#pragma link C++ MACRO SET_OPTION;
#pragma link C++ MACRO PRINT_OPTIONS;
#pragma link C++ MACRO GET_OPTION_VALUE;
#endif


}//close namespace 

#endif
 
