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
* @file Cut.h
* @class Cut
* @brief Define cuts
* 
* @author S. Riggi
* @date 20/6/2019
*/


#ifndef _CUT_h
#define _CUT_h 1

#include <CodeUtils.h>

#include <TObject.h>

#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>

using namespace std;


namespace Caesar {

//=======================================
//==           CUT CLASS
//=======================================
class Cut : public TObject 
{
	public:
		
		/** 
		\brief Cut constructor
 		*/
		Cut(std::string name,bool enabled=false,bool reverse=false,bool combineInOR=false):
			m_name(name), m_enabled(enabled), m_reverse(reverse), m_combineInOR(combineInOR)
		{}
		/** 
		\brief Cut destructor
 		*/
    virtual ~Cut() {}		
		
	public:
		/** 
		\brief Check if cut is passed
 		*/
		virtual bool isPassed(double val){return false;}
		virtual bool isPassed(float val){return false;}
		virtual bool isPassed(int val){return false;}
		virtual bool isPassed(long int val){return false;}	
		virtual bool isPassed(bool val){return false;}	
		/** 
		\brief Enable cut
 		*/
		virtual void enable(){m_enabled=true;}
		/** 
		\brief Disable cut
 		*/
		virtual void disable(){m_enabled=false;}
		/** 
		\brief Is cut enabled?
 		*/
		virtual bool isEnabled(){return m_enabled;}
		/** 
		\brief Print info
 		*/
		virtual void PrintInfo(){}

	public:
		std::string m_name;
		bool m_enabled;	
		bool m_reverse;
		bool m_combineInOR;
		

	ClassDef(Cut,1)	
};
typedef Cut* CutPtr;
typedef std::map<std::string, CutPtr> CutMap;


//=======================================
//==       EQUALITY CUT CLASS
//=======================================
template<typename T>
class EqualityCut : public Cut 
{
	public:
		/** 
		\brief Equality cut constructor
 		*/
		EqualityCut(std::string name,T cutVal,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			m_cutValues.insert(cutVal);
		}
		/** 
		\brief Equality cut constructor
 		*/
		EqualityCut(std::string name,std::initializer_list<T> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false):
			Cut(name,enabled,reverse,combineInOR)
		{
			for(auto item: cutValues) m_cutValues.insert(item);
		}
		/** 
		\brief Equality cut constructor
 		*/
		EqualityCut(std::string name,std::vector<T>& cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false):
			Cut(name,enabled,reverse,combineInOR)
		{
			for(auto item: cutValues) m_cutValues.insert(item);
		}
		/** 
		\brief Equality cut destructor
 		*/
    virtual ~EqualityCut() 
		{
			m_cutValues.clear();
		}
		/** 
		\brief Print info
 		*/
		virtual void PrintInfo()
		{
			std::string equalitySign= "==";
			if(m_reverse) equalitySign= "!=";
			cout<<"== CUT "<<m_name<<" (enabled? "<<m_enabled<<") =="<<endl;
			size_t iter= 0;
			for(auto cutValue: m_cutValues){
				if(m_combineInOR) {
					if(iter==m_cutValues.size()-1) cout<<"x"<<equalitySign<<cutValue<<endl;
					else cout<<"x"<<equalitySign<<cutValue<<" || ";
				}
				else{
					if(iter==m_cutValues.size()-1) cout<<"x"<<equalitySign<<cutValue<<endl;
					else cout<<"x"<<equalitySign<<cutValue<<" && ";
				}
				iter++;
			}
			cout<<"==========================================="<<endl;
		}

		
			
	public:
		/** 
		\brief Check if cut is passed
 		*/
		virtual bool isPassed(double val){return isEqualityCutPassed(val);}
		virtual bool isPassed(float val){return isEqualityCutPassed(val);}
		virtual bool isPassed(int val){return isEqualityCutPassed(val);}
		virtual bool isPassed(long int val){return isEqualityCutPassed(val);}
		virtual bool isPassed(bool val){return isEqualityCutPassed(val);}
		
	private:
		/** 
		\brief Check if cut is passed
 		*/
		template<typename K>
		bool isEqualityCutPassed(K& val)
		{
			//If not enabled pass
			if(!m_enabled) return true;

			//Compute cut pass results
			std::vector<bool> cutResults;
			for (auto cutValue : m_cutValues){
				K cutValue_casted= static_cast<K>(cutValue);
				bool passed= (val==cutValue_casted);
				if(m_reverse) passed= (val!=cutValue_casted);
				cutResults.push_back(passed);
			}
	
			//Perform OR/AND among all cut results
			bool cutPassed= false;
			if(m_combineInOR) cutPassed= CodeUtils::GetVecLogicalOr(cutResults); 
			else cutPassed= CodeUtils::GetVecLogicalAnd(cutResults); 
				
			return cutPassed;
		}

	public:
		std::set<T> m_cutValues;

	ClassDef(EqualityCut,1)	
};


//=======================================
//==       BOUND CUT CLASS
//=======================================
template<typename T>
class BoundCut : public Cut 
{
	public:
		/** 
		\brief Bound cut constructor
 		*/
		BoundCut(std::string name,T cutVal_min,T cutVal_max,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			m_cutValues.insert(std::make_pair<cutVal_min,cutVal_max>);
		}
		/** 
		\brief Bound cut constructor
 		*/
		BoundCut(std::string name,std::initializer_list<std::pair<T,T>> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			for(auto item: cutValues) m_cutValues.insert(item);
		}
		/** 
		\brief Bound cut constructor
 		*/
		BoundCut(std::string name,std::vector<std::pair<T,T>> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			for(auto item: cutValues) m_cutValues.insert(item);
		}
		/** 
		\brief Bound cut destructor
 		*/
    virtual ~BoundCut() 
		{
			m_cutValues.clear();
		}
	
	public:
		/** 
		\brief Check if cut is passed
 		*/
		virtual bool isPassed(double val){return isBoundCutPassed(val);}
		virtual bool isPassed(float val){return isBoundCutPassed(val);}
		virtual bool isPassed(int val){return isBoundCutPassed(val);}
		virtual bool isPassed(long int val){return isBoundCutPassed(val);}
		virtual bool isPassed(bool val){return isBoundCutPassed(val);}

	private:
		/** 
		\brief Check if cut is passed
 		*/
		template<typename K>
		bool isBoundCutPassed(K& val)
		{
			//If not enabled pass
			if(!m_enabled) return true;

			//Compute cut pass results
			std::vector<bool> cutResults;
			for (auto cutValue : m_cutValues)
			{
				K cutMinValue_casted= static_cast<K>(cutValue.first);
				K cutMaxValue_casted= static_cast<K>(cutValue.second);
				bool passed= (val>=cutMinValue_casted && val<=cutMaxValue_casted);	
				if(m_reverse) passed= (val<=cutMinValue_casted && val>=cutMaxValue_casted);
				cutResults.push_back(passed);
			}
			
			//Perform OR/AND among all cut results
			bool cutPassed= false;
			if(m_combineInOR) cutPassed= CodeUtils::GetVecLogicalOr(cutResults);
			else cutPassed= CodeUtils::GetVecLogicalAnd(cutResults);
			
			return cutPassed;
		}

		/** 
		\brief Print info
 		*/
		virtual void PrintInfo()
		{
			std::string boundSign1= "<=";
			std::string boundSign2= "<=";
			if(m_reverse) {
				boundSign1= ">=";
				boundSign2= ">=";
			}

			cout<<"== CUT "<<m_name<<" (enabled? "<<m_enabled<<") =="<<endl;
			size_t iter= 0;
			for(auto cutValue: m_cutValues){
				if(m_combineInOR) {
					if(iter==m_cutValues.size()-1) cout<<cutValue.first<<boundSign1<<"x"<<boundSign2<<cutValue.second<<endl;
					else cout<<cutValue.first<<boundSign1<<"x"<<boundSign2<<cutValue.second<<" || ";
				}
				else{
					if(iter==m_cutValues.size()-1) cout<<cutValue.first<<boundSign1<<"x"<<boundSign2<<cutValue.second<<endl;
					else cout<<cutValue.first<<boundSign1<<"x"<<boundSign2<<cutValue.second<<" && ";
				}
				iter++;
			}
			cout<<"==========================================="<<endl;
		}

	public:
		std::set<std::pair<T,T>> m_cutValues;
		

	ClassDef(BoundCut,1)	
};


//=======================================
//==       SINGLE BOUND CUT CLASS
//=======================================
template<typename T>
class SingleBoundCut : public Cut 
{
	public:
		/** 
		\brief Single bound cut constructor
 		*/
		SingleBoundCut(std::string name,T cutVal,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			m_cutValues.insert(cutVal);
		}
		/** 
		\brief Single bound cut constructor
 		*/
		SingleBoundCut(std::string name,std::initializer_list<T> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			for(auto item: cutValues) m_cutValues.insert(item);
		}
		/** 
		\brief Single bound cut constructor
 		*/
		SingleBoundCut(std::string name,std::vector<T> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false):	
			Cut(name,enabled,reverse,combineInOR)
		{
			for(auto item: cutValues) m_cutValues.insert(item);
		}
		/** 
		\brief Single bound cut destructor
 		*/
    virtual ~SingleBoundCut() 
		{
			m_cutValues.clear();
		}
	
	public:
		/** 
		\brief Check if cut is passed
 		*/
		virtual bool isPassed(double val){return isSingleBoundCutPassed(val);}
		virtual bool isPassed(float val){return isSingleBoundCutPassed(val);}
		virtual bool isPassed(int val){return isSingleBoundCutPassed(val);}
		virtual bool isPassed(long int val){return isSingleBoundCutPassed(val);}
		virtual bool isPassed(bool val){return isSingleBoundCutPassed(val);}
		
	private:
		/** 
		\brief Check if cut is passed
 		*/
		template<typename K>
		bool isSingleBoundCutPassed(K& val)
		{
			//If not enabled pass
			if(!m_enabled) return true;

			//Compute cut pass results
			std::vector<bool> cutResults;
			for (auto cutValue : m_cutValues)
			{
				K cutValue_casted= static_cast<K>(cutValue);
				bool passed= (val>=cutValue_casted);
				if(m_reverse) passed= (val<=cutValue_casted);
				cutResults.push_back(passed);
			}
			
			//Perform OR/AND among all cut results
			bool cutPassed= false;
			if(m_combineInOR) cutPassed= CodeUtils::GetVecLogicalOr(cutResults); 
			else cutPassed= CodeUtils::GetVecLogicalAnd(cutResults); 
				
			return cutPassed;
		}

		/** 
		\brief Print info
 		*/
		virtual void PrintInfo()
		{
			std::string boundSign= ">=";
			if(m_reverse) {
				boundSign= "<=";
			}

			cout<<"== CUT "<<m_name<<" (enabled? "<<m_enabled<<") =="<<endl;
			size_t iter= 0;
			for(auto cutValue: m_cutValues){
				if(m_combineInOR) {
					if(iter==m_cutValues.size()-1) cout<<"x"<<boundSign<<cutValue<<endl;
					else cout<<"x"<<boundSign<<cutValue<<" || ";
				}
				else{
					if(iter==m_cutValues.size()-1) cout<<"x"<<boundSign<<cutValue<<endl;
					else cout<<"x"<<boundSign<<cutValue<<" && ";
				}
				iter++;
			}
			cout<<"==========================================="<<endl;
		}

	public:
		std::set<T> m_cutValues;
		

	ClassDef(SingleBoundCut,1)	
};


//=======================================
//==       CUT SPEC CLASS
//=======================================
/*
class CutSpec : public TObject 
{
	public:
 	 typedef bool(*CutFunctionPtr)(Cut &);

  	// A cut specification with name, implementation, and number of parameters
  	CutSpec(const std::string cutName, CutFunctionPtr cutFunction, int nParams);
  /// A cut specification with name and implementation
  CutSpec(const std::string cutName, CutFunctionPtr cutFunction);

  /// A cut specification with only a name is useful for the END marker.
  CutSpec(const std::string cutName);
  /// An empty cut specification is only useful for placeholders
  CutSpec();


  // Access the name of the cut
  const std::string & GetCutName() const { return fCutName; }
  // Access the cut implementation function
  CutFunctionPtr GetCutFunction() const { return fCutFunction; }
  
	// Access the number of cut parameters (-1 for 'not defined')
  int GetNParams() const { return fNParams; }

private:
  std::string fCutName;
  CutFunctionPtr fCutFunction;
  int fNParams;
};
*/


//=======================================
//==        PARAMETRIC CUT CLASS
//=======================================
/*
class ParametricCut : public TObject 
{
	public:
		
		
		ParametricCut(std::string name,bool enabled=false,bool reverse=false,bool combineInOR=false):
			m_name(name), m_enabled(enabled), m_reverse(reverse), m_combineInOR(combineInOR)
		{}
		
    virtual ~ParametricCut() {}		
		
	public:
		
		virtual bool isPassed(double val){return false;}
		virtual bool isPassed(float val){return false;}
		virtual bool isPassed(int val){return false;}
		virtual bool isPassed(long int val){return false;}	
		virtual bool isPassed(bool val){return false;}	
		
		virtual void enable(){m_enabled=true;}
		
		virtual void disable(){m_enabled=false;}
		
		virtual bool isEnabled(){return m_enabled;}
		
		virtual void PrintInfo(){}

	public:
		std::string m_name;
		bool m_enabled;	
		bool m_reverse;
		bool m_combineInOR;
		

	ClassDef(ParametricCut,1)	
};
typedef ParametricCut* ParametricCutPtr;
typedef std::map<std::string, ParametricCutPtr> ParametricCutMap;
*/

//=======================================
//==       CUT FACTORY CLASS
//=======================================
class CutFactory : public TObject 
{
	public:
		/** 
		\brief Return cut factory instance
 		*/
		static CutFactory& Instance() 
		{
    	static CutFactory myInstance;
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    CutFactory(CutFactory const&) = delete;             // Copy construct
    CutFactory(CutFactory&&) = delete;                  // Move construct
    CutFactory& operator=(CutFactory const&) = delete;  // Copy assign
    CutFactory& operator=(CutFactory &&) = delete;      // Move assign


	public:

		/** 
		\brief Create an equality cut
 		*/
		template<typename T>
		static CutPtr CreateEqCut(std::string name,T cutValue,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new EqualityCut<T>(name,cutValue,enabled,reverse,combineInOR));
		}
		/** 
		\brief Create an equality cut
 		*/
		template<typename T>
		static CutPtr CreateEqCut(std::string name,std::vector<T>& cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new EqualityCut<T>(name,cutValues,enabled,reverse,combineInOR));
		}
		/** 
		\brief Create an equality cut
 		*/
		template<typename T>
		static CutPtr CreateEqCut(std::string name,std::initializer_list<T> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new EqualityCut<T>(name,cutValues,enabled,reverse,combineInOR));
		}

		/** 
		\brief Register equality cut
 		*/
		template <typename T>
  	int RegisterEqualityCut(std::string name,T cutValue,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{				
			//Check args
			if(name=="") {
				cerr<<"CutFactory()::RegisterEqualityCut(): WARN: Invalid cut name given!"<<endl;					
				return -1;
			}
			
			//Check if option with given name already exists
			CutPtr cut= GetCut(name);
			if(cut){
				cerr<<"CutFactory()::RegisterEqualityCut(): WARN: Cut "<<name<<" already registered, skip it!"<<endl;
				return 0;
			}
	
			//Cut does not exist, create one and add to the map
			cut= CreateEqCut<T>(name,cutValue,enabled,reverse,combineInOR);
			m_registeredCuts.insert( std::pair<std::string,CutPtr>(name,cut) );
		
			return 0;

  	}//close RegisterEqualityCut()

		/** 
		\brief Register equality cut
 		*/
		template <typename T>
  	int RegisterEqualityCut(std::string name,std::vector<T>& cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{				
			//Check args
			if(name=="") {
				cerr<<"CutFactory()::RegisterEqualityCut(): WARN: Invalid cut name given!"<<endl;					
				return -1;
			}
	
			if(cutValues.empty()){
				cerr<<"CutFactory()::RegisterEqualityCut(): WARN: Empty cut parameters given!"<<endl;
				return -1;
			}
			
			//Check if option with given name already exists
			CutPtr cut= GetCut(name);
			if(cut){
				cerr<<"CutFactory()::RegisterEqualityCut(): WARN: Cut "<<name<<" already registered, skip it!"<<endl;
				return 0;
			}
	
			//Cut does not exist, create one and add to the map
			cut= CreateEqCut<T>(name,cutValues,enabled,reverse,combineInOR);
			m_registeredCuts.insert( std::pair<std::string,CutPtr>(name,cut) );
		
			return 0;

  	}//close RegisterEqualityCut()

		

		/** 
		\brief Create a bound cut
 		*/
		template<typename T>
		static CutPtr CreateBoundCut(std::string name,T cutMinValue,T cutMaxValue,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new BoundCut<T>(name,cutMinValue,cutMaxValue,enabled,reverse,combineInOR));
		}
		/** 
		\brief Create a bound cut
 		*/
		template<typename T>
		static CutPtr CreateBoundCut(std::string name,std::vector<std::pair<T,T>>& cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new BoundCut<T>(name,cutValues,enabled,reverse,combineInOR));
		}
		/** 
		\brief Create an equality cut
 		*/
		template<typename T>
		static CutPtr CreateBoundCut(std::string name,std::initializer_list<std::pair<T,T>> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new BoundCut<T>(name,cutValues,enabled,reverse,combineInOR));
		}

		/** 
		\brief Register bound cut
 		*/
		template <typename T>
  	int RegisterBoundCut(std::string name,T cutMinValue,T cutMaxValue,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{				
			//Check args
			if(name=="") {
				cerr<<"CutFactory()::Register(): WARN: Invalid option name given!"<<endl;					
				return -1;
			}
			if(cutMinValue>=cutMaxValue) {
				cerr<<"CutFactory()::Register(): WARN: Invalid bound cut parameter given!"<<endl;					
				return -1;
			}

			//Check if option with given name already exists
			CutPtr cut= GetCut(name);
			if(cut){
				cerr<<"CutFactory()::RegisterBoundCut(): WARN: Cut "<<name<<" already registered, skip it!"<<endl;
				return 0;
			}
	
			//Cut does not exist, create one and add to the map
			cut= CreateBoundCut<T>(name,cutMinValue,cutMaxValue,enabled,reverse,combineInOR);
			m_registeredCuts.insert( std::pair<std::string,CutPtr>(name,cut) );
		
			return 0;

  	}//close RegisterBoundCut()


		/** 
		\brief Register bound cut
 		*/
		template <typename T>
  	int RegisterBoundCut(std::string name,std::vector<std::pair<T,T>> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{				
			//Check args
			if(name=="") {
				cerr<<"CutFactory()::RegisterBoundCut(): WARN: Invalid option name given!"<<endl;					
				return -1;
			}
			if(cutValues.empty()){
				cerr<<"CutFactory()::RegisterBoundCut(): WARN: Empty cut parameters given!"<<endl;
				return -1;
			}

			//Check if option with given name already exists
			CutPtr cut= GetCut(name);
			if(cut){
				cerr<<"CutFactory()::RegisterBoundCut(): WARN: Cut "<<name<<" already registered, skip it!"<<endl;
				return 0;
			}
	
			//Cut does not exist, create one and add to the map
			cut= CreateBoundCut<T>(name,cutValues,enabled,reverse,combineInOR);
			m_registeredCuts.insert( std::pair<std::string,CutPtr>(name,cut) );
		
			return 0;

  	}//close RegisterBoundCut()
		
		/** 
		\brief Create a single bound cut
 		*/
		template<typename T>
		static CutPtr CreateSingleBoundCut(std::string name,T cutValue,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new SingleBoundCut<T>(name,cutValue,enabled,reverse,combineInOR));
		}
		/** 
		\brief Create a single bound cut
 		*/
		template<typename T>
		static CutPtr CreateSingleBoundCut(std::string name,std::vector<T>& cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{	
			return CutPtr(new SingleBoundCut<T>(name,cutValues,enabled,reverse,combineInOR));
		}

		/** 
		\brief Register single bound cut
 		*/
		template <typename T>
  	int RegisterSingleBoundCut(std::string name,T cutValue,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{				
			//Check args
			if(name=="") {
				cerr<<"CutFactory()::RegisterSingleBoundCut(): WARN: Invalid option name given!"<<endl;					
				return -1;
			}

			//Check if option with given name already exists
			CutPtr cut= GetCut(name);
			if(cut){
				cerr<<"CutFactory()::RegisterSingleBoundCut(): WARN: Cut "<<name<<" already registered, skip it!"<<endl;
				return 0;
			}
	
			//Cut does not exist, create one and add to the map
			cut= CreateSingleBoundCut<T>(name,cutValue,enabled,reverse,combineInOR);
			m_registeredCuts.insert( std::pair<std::string,CutPtr>(name,cut) );
		
			return 0;

  	}//close RegisterSingleBoundCut()

		/** 
		\brief Register single bound cut
 		*/
		template <typename T>
  	int RegisterSingleBoundCut(std::string name,std::vector<T> cutValues,bool enabled=false,bool reverse=false,bool combineInOR=false) 
		{				
			//Check args
			if(name=="") {
				cerr<<"CutFactory()::RegisterSingleBoundCut(): WARN: Invalid option name given!"<<endl;					
				return -1;
			}
			if(cutValues.empty()){
				cerr<<"CutFactory()::RegisterSingleBoundCut(): WARN: Empty cut parameters given!"<<endl;
				return -1;
			}

			//Check if option with given name already exists
			CutPtr cut= GetCut(name);
			if(cut){
				cerr<<"CutFactory()::RegisterSingleBoundCut(): WARN: Cut "<<name<<" already registered, skip it!"<<endl;
				return 0;
			}
	
			//Cut does not exist, create one and add to the map
			cut= CreateSingleBoundCut<T>(name,cutValues,enabled,reverse,combineInOR);
			m_registeredCuts.insert( std::pair<std::string,CutPtr>(name,cut) );
		
			return 0;

  	}//close RegisterSingleBoundCut()

		/** 
		\brief Has cut?
 		*/
		bool HasCut(std::string name)
		{
			CutMap::iterator it= m_registeredCuts.find(name);
			if ( m_registeredCuts.empty() || it==m_registeredCuts.end() ) return false;
			return true;
		}

		/** 
		\brief Get cut pointer
 		*/
		CutPtr GetCut(std::string name)
		{
			CutMap::iterator it= m_registeredCuts.find(name);
			if ( m_registeredCuts.empty() || it==m_registeredCuts.end() ) return nullptr;
			return it->second;
		}

		/** 
		\brief Is cut passed
 		*/
		template<typename T>
		bool IsCutPassed(std::string name,T val)
		{
			CutPtr cut= GetCut(name);
			if(!cut){
				cerr<<"CutFactory()::RegisterBoundCut(): WARN: Cut "<<name<<" not found, returning false!"<<endl;
				return false;
			}
			return cut->isPassed(val);
		}
		
		/** 
		\brief Print registered cuts
 		*/
		void PrintCuts()
		{
			for (CutMap::const_iterator it = m_registeredCuts.begin(); it!=m_registeredCuts.end(); it++){
    		(it->second)->PrintInfo();
				cout<<endl;
			}
		}
		
		/** 
		\brief Returns cut map
 		*/
		CutMap GetCuts(){return m_registeredCuts;}

	protected:

		/** 
		\brief Protected cut factory constructor
 		*/	
		CutFactory(){};
		/** 
		\brief Protected cut factory destructor
 		*/
		virtual ~CutFactory()
		{
			for (CutMap::const_iterator it = m_registeredCuts.begin(); it!=m_registeredCuts.end(); it++){
				if(it->second) delete it->second;
			}
			m_registeredCuts.clear();
		};

	private:

		CutMap m_registeredCuts;		
	
	ClassDef(CutFactory,1)	
};



#ifdef __MAKECINT__
#pragma link C++ class Cut+;
#pragma link C++ class EqualityCut+;
#pragma link C++ class EqualityCut<int>+;
#pragma link C++ class EqualityCut<long int>+;
#pragma link C++ class EqualityCut<float>+;
#pragma link C++ class EqualityCut<double>+;
#pragma link C++ class EqualityCut<bool>+;
#pragma link C++ class BoundCut+;
#pragma link C++ class BoundCut<int>+;
#pragma link C++ class BoundCut<long int>+;
#pragma link C++ class BoundCut<float>+;
#pragma link C++ class BoundCut<double>+;
#pragma link C++ class BoundCut<bool>+;
#pragma link C++ class SingleBoundCut+;
#pragma link C++ class SingleBoundCut<int>+;
#pragma link C++ class SingleBoundCut<long int>+;
#pragma link C++ class SingleBoundCut<float>+;
#pragma link C++ class SingleBoundCut<double>+;
#pragma link C++ class SingleBoundCut<bool>+;
#pragma link C++ class CutFactory+;
#endif

}//close namespace 

#endif

