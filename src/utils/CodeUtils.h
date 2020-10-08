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

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <WCSUtils.h>

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
#include <random>
#include <numeric>

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

		/**
		* \brief Generate uuid string
		*/
		static std::string GenerateUUID(){
			boost::uuids::uuid random_uuid = boost::uuids::random_generator()();
			return boost::lexical_cast<std::string>(random_uuid);
		}

		/**
		* \brief Convert json to string (DEPRECATED JSONCPP API)
		*/
		/*
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
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to encode to json string!");
				#endif
				return -1;
			}		
			return 0;
		}//close JsonToString()
		*/

		/**
		* \brief Concatenate integers and create string code given digits
		*/
		static std::string GetStringCodeFromIntegers(std::vector<int> ids,size_t ndigits)
		{
			std::string code= "";
			if(ndigits<=0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("ndigits must be >0");
				#endif
				return code;
			}

			for(size_t i=0;i<ids.size();i++)
			{
				std::stringstream ss;
				ss<<std::setfill('0')<<std::setw(ndigits)<<ids[i];
				code+= ss.str();
			}

			return code;
		}

		/**
		* \brief Concatenate integers and create int code given digits
		*/
		static long int GetCodeFromIntegers(std::vector<int> ids,size_t ndigits)
		{
			std::string code_str= GetStringCodeFromIntegers(ids,ndigits);
			long int code= stol(code_str.c_str());
			return code;
		}

		/**
		* \brief Split a string in equal parts
		*/
		static std::vector<std::string> SplitStringInEqualParts(const std::string& str, int splitLength)
		{
   		int NumSubstrings = str.length()/splitLength;
   		std::vector<std::string> ret;

   		for (auto i = 0; i < NumSubstrings; i++)
   		{
      	ret.push_back(str.substr(i * splitLength, splitLength));
   		}

   		// If there are leftover characters, create a shorter item at the end.
   		if (str.length() % splitLength != 0)
   		{
      	ret.push_back(str.substr(splitLength * NumSubstrings));
   		}

   		return ret;
		}

		/**
		* \brief Get integer codes from string codes
		*/
		static int DecodeIntCodes(std::vector<int>& ids,long int code,size_t ndigits)
		{
			ids.clear();
			if(ndigits<=0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("ndigits must be >0");
				#endif
				return -1;
			}
			std::string code_str= to_string(code);
			int n= code_str.size();
			int ntrailing_zeros= ndigits - n%ndigits;
			std::stringstream ss;
			if(n>0) code_str= std::string(ntrailing_zeros,'0').append(code_str);
			
			std::vector<std::string> codes_str= SplitStringInEqualParts(code_str,ndigits);
			
			for(size_t i=0;i<codes_str.size();i++){
				int thisCode= stoi(codes_str[i].c_str());	
				ids.push_back(thisCode);
			}
	
			return 0;
		}

		/**
		* \brief Convert json to string (NEW JSONCPP API)
		*/
		static int JsonToString(std::string& jsonString,Json::Value& jsonObj,bool isMinified=true)
		{
			//Encode to string
			Json::StreamWriterBuilder builder;
			builder["commentStyle"] = "None";
			if(isMinified) builder["indentation"] = "";//written in a single line			
			else builder["indentation"] = "  ";

			try {
				jsonString= Json::writeString(builder, jsonObj);
			}
			catch(...){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to encode to json string!");
				#endif
				return -1;
			}
			return 0;
		}//close JsonToString()

		/**
		* \brief Convert string to json
		*/
		static int StringToJson(Json::Value& root,std::string& jsonString)
		{		
			//DEPRECATED JSONCPP API
			/*	
			Json::Reader reader;
			if(!reader.parse(jsonString, root)) {
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to encode string to json ("<<reader.getFormattedErrorMessages()<<")");
				#endif
				return -1;
			}
			*/

			//NEW JSONCPP API
			Json::CharReaderBuilder builder;
			Json::CharReader* reader = builder.newCharReader();
			Json::Value output;
			std::string errors;

			bool status= reader->parse(jsonString.c_str(), jsonString.c_str() + jsonString.length(), &root, &errors);
			if(reader){
				delete reader;
				reader= 0;
			}
			if(!status){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to encode string to json (err="<<errors<<")");
				#endif
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
				if(!markedElements[i]) tempBuffer.push_back(data[i]);
    	}
    	data = tempBuffer;

		}//close DeleteItems()


		/**
		* \brief Delete selected items from a vector
		*/
		template<class T,typename K>
			static void DeletePtrItems(std::vector<T>& data, const std::vector<K>& deleteIndices) {
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
    	}
    	data = tempBuffer;

		}//close DeletePtrItems()

		/**
		* \brief Delete object pointer collection
		*/
		template <class T>
		static void DeletePtrCollection(std::vector<T*>& data){
			for(size_t i=0;i<data.size();i++){
				if(data[i]){
					delete data[i];
					data[i]= 0;
				}
			}
			data.clear();	
		}//close DeletePtrCollection()

		/**
		* \brief Delete object pointer collection
		*/
		template <class T>
		static void DeletePtrCollection(std::initializer_list<T*> data){
			for (auto it : data){
				T* item= it;
				if(item){
					delete item;
					item= 0;
				}
			}
		}//close DeletePtrCollection()

		/**
		* \brief Delete object pointer 
		*/
		template <class T>
		static void DeletePtr(T* data){
			if(data){
				delete data;
				data= 0;
			}
		}//close DeletePtr()
		
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

		
		/**
		* \brief Reorder vector
		*/
		template< class T >
			static void reorder(std::vector<T> & unordered,std::vector<size_t> const & index_map,std::vector<T> & ordered){
  			// copy for the reorder according to index_map, because unsorted may also be
  			// sorted
  			std::vector<T> copy = unordered;
  			ordered.resize(index_map.size());
  			for(unsigned int i = 0; i<index_map.size();i++)
					ordered[i] = copy[index_map[i]];
			}

		/**
		* \brief Sort vector
		*/
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
	
		/**
		* \brief Sort vector in descending order
		*/
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

		typedef std::vector< std::pair<long int,long int> > IndexPairs;

		/**
		* \brief Find index of equal elements in two vectors. 
		*/
		template <class Iterator,class Comparator>
		static IndexPairs FindIntersectionIndexes(Iterator first1,Iterator last1,Iterator first2,Iterator last2,Comparator comp,bool sorted=false)
		{			
			//Sort vectors if not sorted
			if(!sorted){
				std::sort( first1, last1, comp );
				std::sort( first2, last2, comp );
			}

			//Find intersection indexes
			IndexPairs indexes;
			Iterator begin1= first1;
			Iterator begin2= first2;
	
			while (first1 != last1 && first2 != last2) {
  			if (comp(*first1, *first2)) {
    			++first1;
    		} 
				else if(comp(*first2, *first1)){
					++first2;
				}
				else {
					long int index1= std::distance(begin1,first1);
					long int index2= std::distance(begin2,first2);
					indexes.push_back( std::make_pair(index1,index2) );
					++first1;
					++first2;
   			}
 			}//end while loop

			return indexes;

		}//close FindIntersectionIndexes()
		
		/**
		* \brief Extract a subset of random index without repetitions from a container
		*/
		template<class Iterator>
		static Iterator random_unique(Iterator begin, Iterator end, size_t num_random) 
		{
			//Fisher-Yates shuffle (https://stackoverflow.com/questions/9345087/choose-m-elements-randomly-from-a-vector-containing-n-elements)
    	size_t left = std::distance(begin, end);
    	while (num_random--) {
        Iterator r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
   	 	}
    	return begin;
		}

		/** 
		\brief Extract sample from data vector with given sample size, with/without repetitions and uniform weights
 		*/
		template<typename T>	
		static int ExtractVectorRandomSample(std::vector<T>& sample_data,const std::vector<T>& data,long int n=-1,bool repeate=true)
		{
			//Init sample data
			sample_data.clear();

			//Check data size
			if(data.empty()){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Input data are empty, nothing to be sampled!");
				#endif
				return -1;
			}
			
			//Check if all data are to be used
			long int data_size= static_cast<long int>(data.size());
			if(n<=0 || n>data_size){
				n= data_size;
			}

			//Init random generator
			std::random_device rd; // obtain a random number from hardware
    	std::mt19937 gen(rd()); // seed the generator
    	std::uniform_int_distribution<long int> dis(0,n-1);

			//Check if repetitions are to be used
			long int sample_size= static_cast<long int>(sample_data.size());

			if(repeate){
				for(long int i=0;i<n;i++){
					long int rand_index= dis(gen);
					T dataValue= data[rand_index];
					sample_data.push_back(dataValue);
				}
			}//close if
			else{
				if(n==data_size){//return input data in case no repetition with all data are requested
					sample_data= data;
				}
				else{
					//Copy input data and extract n random unique samples
					std::vector<T> tmp= data;
					random_unique(tmp.begin(),tmp.end(),n);
					for(long int i=0;i<n;i++) sample_data.push_back(tmp[i]);
				}//close else
			}//close else NO REPETITIONS

			return 0;

		}//close ExtractVectorRandomSample()

		/** 
		\brief Extract a number of random samples from data vector with given sample size, with/without repetitions and uniform weights
 		*/
		template<typename T>
		static int ExtractVectorRandomSamples(std::vector<std::vector<T>>& data_samples,const std::vector<T>& data,int nSamples,long int n=-1,bool repeate=true)
		{
			//Init samples
			data_samples.clear();

			//Check data size
			if(data.empty()){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Input data are empty, nothing to be sampled!");
				#endif
				return -1;
			}

			//Generate samples
			for(int i=0;i<nSamples;i++){
				std::vector<T> sample_data;
				int status= ExtractVectorRandomSample(sample_data,data,n,repeate);
				if(status<0){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Failed to generate random sample no. "<<i+1<<", exit generation!");
					#endif
					data_samples.clear();
					return -1;
				}

				//Add sample to collection	
				data_samples.push_back(sample_data);

			}//end loop samples
		
			return 0;

		}//close ExtractVectorRandomSamples()

		/**
		* \brief Find vector index at which the cumulative sum is smaller then given value (comparator version)
		*/
		template <class T,class Comparator>
		static T FindCumulativeSumFractionThr(std::vector<T>& data,Comparator comp,double thr, bool sorted=false)
		{			
			//Sort vector if not sorted
			if(!sorted){
				std::sort( data.begin(), data.end(), comp );
			}

			//Fill vector of cumulative sum
			std::vector<T> cdf;
			std::partial_sum(data.begin(),data.end(),std::back_inserter(cdf));
			double sum= cdf[cdf.size()-1];
	
			//Find element at which the cumulative sum is smaller than the given thr 
			size_t pos= std::lower_bound (cdf.begin(), cdf.end(),thr*sum, comp ) - cdf.begin();

			return data[pos];

		}//close FindCumulativeSumFractionThr()
		
		/**
		* \brief Find vector index at which the cumulative sum is smaller then given value 
		*/
		template <class T>
		static T FindCumulativeSumFractionThr(std::vector<T>& data,double thr, bool sorted=false)
		{			
			//Sort vector if not sorted
			if(!sorted){
				std::sort( data.begin(), data.end() );
			}

			//Fill vector of cumulative sum
			std::vector<T> cdf;
			std::partial_sum(data.begin(),data.end(),std::back_inserter(cdf));
			double sum= cdf[cdf.size()-1];
	
			//Find element at which the cumulative sum is smaller than the given thr 
			size_t pos= std::lower_bound (cdf.begin(), cdf.end(),thr*sum ) - cdf.begin();

			return data[pos];

		}//close FindCumulativeSumFractionThr()
		

		/**
		* \brief Find vector index at which the cumulative sum is smaller then given value 
		*/
		template<class InputIt, class T = typename std::iterator_traits<InputIt>::value_type>
		static T FindVectorMode(InputIt begin, InputIt end,int& nmodes)
		{
    	std::map<T, int> counts;
    	for (InputIt it = begin; it != end; ++it) {
      	if (counts.find(*it) != counts.end()) {
        	++counts[*it];
        }
        else {
        	counts[*it] = 1;
        }
    	}
   	 	auto modeItem= std::max_element(counts.begin(), counts.end(),
      	[] (const std::pair<T, int>& pair1, const std::pair<T, int>& pair2) 
				{return pair1.second < pair2.second;});
		
			T mode= modeItem->first;
			int modeFreq= modeItem->second;
	
			nmodes= std::count_if(
				counts.begin(),
				counts.end(),
				[&modeFreq](std::pair<T,int> p){
					return p.second==modeFreq;
				}
			);
		
			return mode;

		}//close FindVectorMode()

		

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
		* \brief Remove pattern in string
		*/
		static void RemovePatternInString(std::string& str,const std::string pattern)
		{
			StringFindAndReplace(str,pattern,"");
		}

		/**
		* \brief Extract substring
		*/
		static std::string ExtractSubString(const std::string& s, const std::string& pattern, bool extractleft=true){
			size_t pos = s.find(pattern);
			if(pos==std::string::npos) return s;
			if(extractleft) return s.substr(0,pos);
			else return s.substr(pos+1,s.length()-pos);
		}

		/**
		* \brief Split string on whitespaces
		*/
		static std::vector<std::string> SplitStringOnWhitespaces(const std::string& s)
		{
			std::vector<std::string> result;
			std::istringstream iss(s);
			for(std::string s; iss >> s;) result.push_back(s);
			return result;
		}

		/**
		* \brief Split string on pattern
		*/
		static std::vector<std::string> SplitStringOnPattern(const std::string& s,char delim)
		{			
			std::vector<std::string> result;
			std::stringstream ss(s);
    	std::string item;
    	while (std::getline(ss, item, delim)) {
      	if (item.length() > 0) {
        	result.push_back(item);  
        }
    	}
    	return result;		
		}

		/**
		* \brief Convert string vector to typed vector
		*/
		template<typename T>
		static std::vector<T> StringVecToTypedVec(const std::vector<std::string>& vec)
		{
			std::vector<T> vec_typed;
			if(vec.empty()) return vec_typed;
			for(size_t i=0;i<vec.size();i++){
				std::stringstream ss(vec[i]);
				T val= T(0);
				ss >> val;
				vec_typed.push_back(val);
			}
			return vec_typed;
		}

		
		
		/**
		* \brief Find pattern in string
		*/
		static bool HasPatternInString(std::string str,std::string pattern)
		{
			std::size_t found = str.find(pattern);
  		if (found!=std::string::npos) return true;
			return false;
		}

		/**
		* \brief Strip blank spaces from string
		*/
		static int StripBlankSpaces(std::string& s){
			try {
				s.erase( 
					std::remove_if(
						s.begin(),
						s.end(),
						[](char c){ 
      				return std::isspace(static_cast<unsigned char>(c));
      			}
					),
					s.end()
				);
			}//close try
			catch(std::exception& e){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("C++ exception ("<<e.what()<<") occurred while stripping blank spaces from string!");
				#endif
				return -1;
			}
			return 0;
		}//close StripBlankSpaces()

		/**
		* \brief Collapse a collection in a string (equivalent of python join)
		*/
		template <typename Iter>
		static std::string JoinCollection(Iter begin, Iter end, std::string separator="")
		{
  		std::ostringstream result;
  		if (begin != end) result << *begin++;
  		while (begin != end) result << separator << *begin++;
  		return result.str();
		}

		/**
		* \brief Join vectors
		*/
		template <typename T>
		static std::string JoinVec(const std::vector<T>& data,std::string separator="")
		{
			return JoinCollection(data.begin(),data.end(),separator);
		}

		/**
		* \brief Perform logical OR among all vector elements
		*/
		static bool GetVecLogicalOr(const std::vector<bool>& v)
		{
			bool res= std::any_of(v.begin(),v.end(), [](bool elem) { return elem; });
			return res;
		}

		/**
		* \brief Perform logical AND among all vector elements
		*/
		static bool GetVecLogicalAnd(const std::vector<bool>& v)
		{
			bool res= std::all_of(v.begin(),v.end(), [](bool elem) { return elem; });
			return res;
		}
		
		/**
		* \brief Compare string case insensitive 
		*/
		static bool AreEqualStringNoCase(std::string str1,std::string str2)
		{
			std::string s1= str1;
			std::string s2= str2;
			bool areEquals= ( strcasecmp(s1.c_str(),s2.c_str())==0 );
			return areEquals;
		}

		/**
		* \brief Extract filename from path
		*/
		static std::string ExtractFileNameFromPath(const std::string& s,bool strip_extension=false){
			char sep = '/';
   		size_t pos = s.rfind(sep, s.length());
			std::string filename= s;
   		if (pos != string::npos) {
      	filename= s.substr(pos+1, s.length() - pos);
   		}
			if(strip_extension) {
				std::string filename_noext= ExtractSubString(filename,".");
				filename= filename_noext;
			}
   		return filename;
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


/**
* \brief Convert string vector to bool vector (specialization of general method)
*/
template<>
inline std::vector<bool> CodeUtils::StringVecToTypedVec<bool>(const std::vector<std::string>& vec)
{
	std::vector<bool> vec_typed;
	if(vec.empty()) return vec_typed;
	for(size_t i=0;i<vec.size();i++){
		std::stringstream ss(vec[i]);
		bool val= false;
		ss >> std::boolalpha >> val;
		vec_typed.push_back(val);
	}
	return vec_typed;
}


#ifdef __MAKECINT__
#pragma link C++ class CodeUtils+;
#endif	

}//close namespace





#endif 
