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
* @file CodeUtils.h
* @class CodeUtils
* @brief Utility functions for programming shortcut tasks
*
* Utility functions for programming shortcut tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef _CODE_UTILS_h
#define _CODE_UTILS_h 1

#include <Logger.h>

#include <TObject.h>

#include <json/json.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <time.h>
#include <ctime>


using namespace std;

namespace Caesar {


template <typename T>
struct MatchJsonValue {
	MatchJsonValue(const T& value,const std::string& key) 
		: m_value(value),m_key(key)
	{}
 			
	//friend bool operator()(const Json::Value& obj) const {
	bool operator()(const Json::Value& obj) const {		
		if(obj[m_key.c_str()].isNull()) return false;
		return (obj[m_key.c_str()].asString()==m_value);
 	}

 	private:
  	const T& m_value;
		const std::string& m_key;

};//close MatchJsonValue()

template <>
inline bool MatchJsonValue<float>::operator()(const Json::Value& obj) const {	
	if(obj[m_key.c_str()].isNull()) return false;
	return (obj[m_key.c_str()].asFloat()==m_value);
}


template <>
inline bool MatchJsonValue<double>::operator()(const Json::Value& obj) const {	
	if(obj[m_key.c_str()].isNull()) return false;
	return (obj[m_key.c_str()].asDouble()==m_value);
}

template <>
inline bool MatchJsonValue<std::string>::operator()(const Json::Value& obj) const {	
 	if(obj[m_key.c_str()].isNull()) return false;
	return (obj[m_key.c_str()].asString()==m_value);
}

template <>
inline bool MatchJsonValue<bool>::operator()(const Json::Value& obj) const {	
	if(obj[m_key.c_str()].isNull()) return false;
	return (obj[m_key.c_str()].asBool()==m_value);
}
		
template <>
inline bool MatchJsonValue<int>::operator()(const Json::Value& obj) const {	
	if(obj[m_key.c_str()].isNull()) return false;
	return (obj[m_key.c_str()].asInt()==m_value);
}
		
template <>
inline bool MatchJsonValue<long int>::operator()(const Json::Value& obj) const {	
 	if(obj[m_key.c_str()].isNull()) return false;
	return (obj[m_key.c_str()].asInt64()==m_value);
}

template <>
inline bool MatchJsonValue<char*>::operator()(const Json::Value& obj) const {	
	if(obj[m_key.c_str()].isNull()) return false;
	return (strcmp(obj[m_key.c_str()].asCString(),m_value)==0);
}



//===============================================
//==          CODE UTILS CLASS
//===============================================
class CodeUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    CodeUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~CodeUtils();

		
	public:

		static std::string GenerateUUID(){
			boost::uuids::uuid random_uuid = boost::uuids::random_generator()();
			return boost::lexical_cast<std::string>(random_uuid);
		}

		static int JsonToString(std::string& jsonString,Json::Value& jsonObj,bool isMinified=true){
			//Encode to string
			try {
				if(isMinified){// write in a minified way
					Json::FastWriter fastWriter;
					jsonString= fastWriter.write(jsonObj);
				}
				else{	// write in a nice readible way
					Json::StyledWriter formattedWriter;
					jsonString= formattedWriter.write(jsonObj);
				}
			}//close try
			catch(...){
				ERROR_LOG("Failed to encode to json string!");
				return -1;
			}		
			return 0;
		}//close JsonToString()

		static int StringToJson(Json::Value& root,std::string& jsonString){
			
			Json::Reader reader;
			if(!reader.parse(jsonString, root)) {
				ERROR_LOG("Failed to encode string to json ("<<reader.getFormattedErrorMessages()<<")");
				return -1;
			}
			return 0;
		}//close StringToJson()

		/**
		* \brief Find json option in array by name
		*/
		template<typename T>
		static int FindJsonValue(int& pos,Json::Value& root,T value,std::string key="name"){
			pos= -1;
			
			//Check option name
			if(key=="") return -1;
			
			//Check json (it shall be an array)
			if(root.isNull() || root.empty() || !root.isArray() ){
				return -1;
			}
				
			//Find option with given key
			Json::ValueIterator it= std::find_if(root.begin(),root.end(),MatchJsonValue<T>(value,key));
			if(it==root.end()){
				pos= -1;
				return -1;
			}	
			pos= (int)(it-root.begin());
			return 0;
		}//close FindJsonOption()

		/**
		* \brief Find item in a vector and returns item position
		*/
		template<class T>
			static bool FindItem(std::vector<T>& v, T itemValue,int& pos) {
				if(v.size()<=0){
					pos= -1;
					return false;
				}	
				typename std::vector<T>::iterator it= std::find(v.begin(),v.end(),itemValue);
				if(it==v.end()){
					pos= -1;
					return false;
				}	
				pos= (int)(it-v.begin());
				return true;
			}//close FindItem()


		/**
		* \brief Delete selected items from a vector
		*/
		template<class T,typename K>
			static void DeleteItems(std::vector<T>& data, const std::vector<K>& deleteIndices) {
    	std::vector<bool> markedElements(data.size(), false);
    	std::vector<T> tempBuffer;
    	tempBuffer.reserve(data.size()-deleteIndices.size());

			typedef typename std::vector<K>::const_iterator Iter;
   	 	for (Iter itDel = deleteIndices.begin(); itDel != deleteIndices.end(); itDel++)
      	markedElements[*itDel] = true;

    	for (size_t i=0; i<data.size(); i++) {
				if (markedElements[i]) {//Free memory!
					if(data[i]){
						delete data[i];
						data[i]= 0;
					}			
				}
				else{
					tempBuffer.push_back(data[i]);
				}
				/*
      	if (!markedElements[i]) {
					tempBuffer.push_back(data[i]);
				}
				*/
    	}
    	data = tempBuffer;
		}//close DeleteItems()


		/**
		* \brief Order vectors and get ordering index
		*/
		template<class T> struct index_cmp{

  		index_cmp(const T arr) : arr(arr) {}
  		bool operator()(const size_t a, const size_t b) const
 			{
    		return arr[a] < arr[b];
  		}
  		const T arr;
		};
		template<class T> struct descending_index_cmp{

  		descending_index_cmp(const T arr) : arr(arr) {}
  		bool operator()(const size_t a, const size_t b) const
 			{
    		return arr[a] > arr[b];
  		}
  		const T arr;
		};

		template< class T >
			static void reorder(std::vector<T> & unordered,std::vector<size_t> const & index_map,std::vector<T> & ordered){
  			// copy for the reorder according to index_map, because unsorted may also be
  			// sorted
  			std::vector<T> copy = unordered;
  			ordered.resize(index_map.size());
  			for(unsigned int i = 0; i<index_map.size();i++)
					ordered[i] = copy[index_map[i]];
			}

		template <class T>
			static void sort(std::vector<T> & unsorted,std::vector<T> & sorted,std::vector<size_t> & index_map){
  			// Original unsorted index map
  			index_map.resize(unsorted.size());
 				for(size_t i=0;i<unsorted.size();i++)
					index_map[i] = i;
  
  			// Sort the index map, using unsorted for comparison
  			std::sort(index_map.begin(),index_map.end(),index_cmp<std::vector<T>& >(unsorted));
  			sorted.resize(unsorted.size());
  			reorder(unsorted,index_map,sorted);
			}
	
		template <class T>
			static void sort_descending(std::vector<T> & unsorted,std::vector<T> & sorted,std::vector<size_t> & index_map){
  			// Original unsorted index map
  			index_map.resize(unsorted.size());
 				for(size_t i=0;i<unsorted.size();i++)
					index_map[i] = i;
  
  			// Sort the index map, using unsorted for comparison
  			std::sort(index_map.begin(),index_map.end(),descending_index_cmp<std::vector<T>& >(unsorted));
  			sorted.resize(unsorted.size());
  			reorder(unsorted,index_map,sorted);
			}

		/**
		* \brief String find and replace
		*/
		static void StringFindAndReplace(std::string& str, const std::string& oldstr, const std::string& newstr){
  		size_t pos = 0;
  		while((pos = str.find(oldstr, pos)) != std::string::npos){
    		str.replace(pos, oldstr.length(), newstr);
     		pos += newstr.length();
  		}	
		}//close StringFindAndReplace()

		/**
		* \brief Remove substring
		*/
		static std::string ExtractSubString(const std::string& s, const std::string& pattern, bool extractleft=true){
			size_t pos = s.find(pattern);
			if(pos==std::string::npos) return s;
			if(extractleft) return s.substr(0,pos);
			else return s.substr(pos+1,s.length()-pos);
		}

		// return an evenly spaced 1-d grid of doubles.
		template <typename T>
			static std::vector<T> linspace(T first,T last,int len){
				std::vector<T> result(len);
  			double step = (last-first) / (len - 1);
  			for (int i=0; i<len; i++) { result[i] = first + i*step; }
  			return result;
			}


		template< class InputIterator, class Function, class Predicate >
    static Function for_each_if(InputIterator first, 
    					 InputIterator last, 
    					 Predicate pred, 
    					 Function f)
    {
    	for( ; first != last; ++first)
    	{
    		if( pred(*first) )
    			f(*first);
    	}
    	return f;
    };

		/** 
		\brief Check if map has key
 		*/
		template <typename T,typename K>
		static bool HasMapKey(std::map<T,K> const& m, T key){
			if(m.empty() || m.find(key)==m.end()) return false;
			return true;
		}

	private:
	
		ClassDef(CodeUtils,1)
};



#ifdef __MAKECINT__
#pragma link C++ class CodeUtils+;
#endif	

}//close namespace


#endif 
