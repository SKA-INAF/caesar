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
* @file SpectrumUtils.h
* @class SpectrumUtils
* @brief Spectrum utilities
* 
* @author S. Riggi
* @date 03/07/2019
*/

#ifndef _SPECTRUM_UTILS_H_
#define _SPECTRUM_UTILS_H_


#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TF1.h>
#include <TF2.h>

#include <iostream>


//=================================================
//===             RESO MODEL 
//==================================================
enum ResoModel {
	eCONST_RESO= 0,
	eEXP_RESO= 1,
	eEXP_STEP_RESO= 2
};


class ResoPars 
{
	public:
		ResoPars(){};
		virtual ~ResoPars() {};

	public:
		//Pure virtual
		virtual double GetA() const = 0;	
		virtual double GetB() const = 0;	
		virtual double GetC() const = 0;
		virtual std::vector<double> GetPars() const = 0;
		
		//Standard
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return nPars;}
		static int GetParNumber() {return 0;}
		
	protected:
		static int nPars;
		int model;
		double A;
		double B;
		double C;
};

class ConstResoPars : public ResoPars 
{
	public:		
		ConstResoPars(double m_A)	{ 
			A= m_A;
			model= eCONST_RESO;
			nPars= 1;	
		};
		virtual ~ConstResoPars() {};

	public:	
		static int GetParNumber() {return 1;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return 0;}
		virtual double GetC() const {return 0;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			return pars;
		}
		
};//close class ConstResoPars

class ExpResoPars : public ResoPars 
{
	public:		
		ExpResoPars(double m_A,double m_B)	{
			A= m_A;
			B= m_B;
			model= eEXP_RESO;
			nPars= 2;	
		};
		virtual ~ExpResoPars() {};

	public:
		static int GetParNumber() {return 2;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return B;}
		virtual double GetC() const {return 0;}
		
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			pars.push_back(B);
			return pars;
		}	
	
};//close class ExpResoPars


class ExpStepResoPars : public ResoPars 
{
	public:		
		ExpStepResoPars(double m_A,double m_B,double m_C){ 
			A= m_A; 
			B= m_B; 
			C= m_C;
			model= eEXP_STEP_RESO;		
			nPars= 3;
		};
		virtual ~ExpStepResoPars() {};

	public:
		static int GetParNumber() {return 3;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return B;}
		virtual double GetC() const {return C;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			pars.push_back(B);
			pars.push_back(C);
			return pars;
		}

};//close ExpStepBiasPars class
//==================================================



//=================================================
//===             BIAS/RESO MODEL 
//==================================================
enum BiasModel {
	eCONST_BIAS= 0,
	eEXP_BIAS= 1,
	eEXP_STEP_BIAS= 2
};

class BiasPars 
{
	public:
		virtual ~BiasPars() {};

	public:
		//Pure virtual
		virtual double GetA() const = 0;	
		virtual double GetB() const = 0;	
		virtual double GetC() const = 0;
		virtual double GetD() const = 0;	
		virtual std::vector<double> GetPars() const = 0;
		
		//Standard
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return nPars;}
		static int GetParNumber() {return 0;}
		
	protected:
		static int nPars;
		int model;
		double A;
		double B;
		double C;
		double D;
	
};//close BiasPars class

class ConstBiasPars : public BiasPars 
{
	public:		
		ConstBiasPars(double m_A)	{ 
			A= m_A;
			model= eCONST_BIAS;
			nPars= 1;	
		};
		virtual ~ConstBiasPars() {};

	public:	
		static int GetParNumber() {return 1;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return 0;}
		virtual double GetC() const {return 0;}
		virtual double GetD() const {return 0;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			return pars;
		}
	
};//close ConstBiasPars class


class ExpBiasPars : public BiasPars 
{
	public:		
		ExpBiasPars(double m_A,double m_B,double m_C){ 
			A= m_A; 
			B= m_B; 
			C= m_C;
			model= eEXP_BIAS;		
			nPars= 3;
		};
		virtual ~ExpBiasPars() {};

	public:
		static int GetParNumber() {return 3;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return B;}
		virtual double GetC() const {return C;}
		virtual double GetD() const {return 0;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			pars.push_back(B);
			pars.push_back(C);
			return pars;
		}

};//close ExpBias class


class ExpStepBiasPars : public BiasPars 
{
	public:		
		ExpStepBiasPars(double m_A,double m_B,double m_C,double m_D){ 
			A= m_A; 
			B= m_B; 
			C= m_C;
			D= m_D;
			model= eEXP_STEP_BIAS;		
			nPars= 4;
		};
		virtual ~ExpStepBiasPars() {};

	public:
		static int GetParNumber() {return 4;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return B;}
		virtual double GetC() const {return C;}
		virtual double GetD() const {return D;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			pars.push_back(B);
			pars.push_back(C);
			pars.push_back(D);
			return pars;
		}

};//close ExpStepBiasPars class
//==================================================


//=================================================
//===             EFFICIENCY MODEL 
//==================================================
enum EfficiencyModel {
	eCONST_EFF= 0,
	eSIGMOID_EFF= 1
};

class EfficiencyPars 
{
	public:
		virtual ~EfficiencyPars() {};
	
	public:
		//Pure virtual
		virtual double GetNorm() const = 0;
		virtual double GetBreak() const = 0;
		virtual double GetSmooth() const = 0;
		virtual std::vector<double> GetPars() const = 0;

		//Default
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return nPars;}
		static int GetParNumber() {return nPars;}

	protected:
		int model;
		static int nPars;
		double Norm;
		double Break;
		double Smooth;

};//close EfficiencyPars class


class ConstEfficiencyPars : public EfficiencyPars 
{
	public:
		ConstEfficiencyPars(double m_Norm){ 
			Norm= m_Norm;
			model= eCONST_EFF;		
			nPars= 1;
		};
		virtual ~ConstEfficiencyPars() {};

	public:
		static int GetParNumber() {return 1;}
		virtual double GetNorm() const {return Norm;}
		virtual double GetBreak() const {return 0;}	
		virtual double GetSmooth() const {return 0;}	
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(Norm);
			return pars;
		}
	
};//close ConstEfficiencyPars class


class SigmoidEfficiencyPars : public EfficiencyPars 
{
	public:
		SigmoidEfficiencyPars(double m_Norm,double m_Break,double m_Smooth){ 
			Norm= m_Norm; 
			Break= m_Break; 	
			Smooth= m_Smooth;
			model= eSIGMOID_EFF;		
			nPars= 3;
		};
		virtual ~SigmoidEfficiencyPars() {};

	public:
		static int GetParNumber() {return 3;}
		virtual double GetNorm() const {return Norm;}
		virtual double GetBreak() const {return Break;}	
		virtual double GetSmooth() const {return Smooth;}	
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(Norm);
			pars.push_back(Break);
			pars.push_back(Smooth);
			return pars;
		}
	
};

//=================================================
//===      SOURCE COUNTS RESOLUTION BIAS MODEL 
//==================================================
class SourceCountsResoBiasPars 
{
	public:
		SourceCountsResoBiasPars(){
			SetDefaults();
		}
		SourceCountsResoBiasPars(std::vector<double> pars)
		{ 
			if(pars.size()!=12){
				std::cerr<<"WARN: Given par size !=12 in SourceCountsResoBiasPars, using default pars!"<<std::endl;
				SetDefaults();
			}
			else{
				nPars= 12;
				bmaj= pars[0];		
				bmin= pars[1];
				rms= pars[2];
				Sthr= pars[3];
				resolvedSourceThrP0= pars[4];
				resolvedSourceThrP1= pars[5];
				resolvedSourceThrSlope= pars[6];
				nu= pars[7];
				alpha= pars[8];
				phiMedianSlope= pars[9];
				phiMedianSbreak= pars[10];
				phiMedianScaleFactor= pars[11];
			}
		};
		virtual ~SourceCountsResoBiasPars() {};

	public:
		static int GetParNumber() {return nPars;}
		virtual double GetBmaj() const {return bmaj;}
		virtual double GetBmin() const {return bmin;}
		virtual double GetRMS() const {return rms;}
		virtual double GetSthr() const {return Sthr;}
		virtual double GetResolvedSourceP0() const {return resolvedSourceThrP0;}
		virtual double GetResolvedSourceP1() const {return resolvedSourceThrP1;}	
		virtual double GetResolvedSourceSlope() const {return resolvedSourceThrSlope;}	
		virtual double GetNu() const {return nu;}	
		virtual double GetAlpha() const {return alpha;}	
		virtual double GetPhiMedianSlope() const {return phiMedianSlope;}	
		virtual double GetPhiMedianSbreak() const {return phiMedianSbreak;}	
		virtual double GetPhiMedianScaleFactor() const {return phiMedianScaleFactor;}	
		
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(bmaj);
			pars.push_back(bmin);
			pars.push_back(rms);
			pars.push_back(Sthr);
			pars.push_back(resolvedSourceThrP0);
			pars.push_back(resolvedSourceThrP1);
			pars.push_back(resolvedSourceThrSlope);
			pars.push_back(nu);
			pars.push_back(alpha);
			pars.push_back(phiMedianSlope);
			pars.push_back(phiMedianSbreak);
			pars.push_back(phiMedianScaleFactor);
			return pars;
		}

	protected:
		void SetDefaults()
		{
			nPars= 12;
			bmaj= 24.;//arcsec
			bmin= 20.;//arcsec
			rms= 300.e-3;//mJy
			Sthr= 1.5;//mJy
			resolvedSourceThrP0= 1.08;
			resolvedSourceThrP1= 2.03;
			resolvedSourceThrSlope= 1;
			nu= 0.912;//GHz
			alpha= -0.9;
			phiMedianSlope= 0.30;
			phiMedianSbreak= 1;//mJy
			phiMedianScaleFactor= 1;
		}

	protected:
		
		static int nPars;
		double bmaj;
		double bmin;
		double rms;	
		double Sthr;
		double resolvedSourceThrP0;
		double resolvedSourceThrP1;
		double resolvedSourceThrSlope;
		double nu;
		double alpha;
		double phiMedianSlope;
		double phiMedianSbreak;
		double phiMedianScaleFactor;

};//close SourceCountsResoBiasPars class

//======================================================
//==        SPECTRUM FITTING PARS
//======================================================
class FitPar 
{
	public:
		FitPar(){};
		FitPar(std::string name,double val)
			: m_Name(name), m_Value(val)
		{
			m_StepSize= 0.1;//10% of value
			m_IsFixed= false;
			m_IsLimited= false;
			m_MinValue= -1;
			m_MaxValue= -1;
			m_ValueError= 0;
		}
		FitPar(std::string name,double val,double minval,double maxval)
			: m_Name(name), m_Value(val), m_MinValue(minval), m_MaxValue(maxval)
		{
			m_StepSize= 0.1;//10% of value
			m_IsFixed= false;
			m_IsLimited= true;
			m_ValueError= 0;
		}
		virtual ~FitPar() {};

	public:
	
		//Getters
		double GetValue() const {return m_Value;}
		double GetValueError() const {return m_ValueError;}
		std::string GetParName() const {return m_Name;}
		bool IsFixed() const {return m_IsFixed;}
		bool IsLimited() const {return m_IsLimited;}
		double GetMinValue() const {return m_MinValue;}
		double GetMaxValue() const {return m_MaxValue;}
		double GetStepSize() const {return m_StepSize;}

		//Setters
		void SetParName(std::string name){m_Name=name;}
		void SetValue(double val){m_Value=val;}
		void SetValueError(double val){m_ValueError=val;}
		void Fix(){m_IsFixed=true;}
		void Free(){m_IsFixed=false;}
		void SetLimits(double minval,double maxval){
			m_MinValue= minval;
			m_MaxValue= maxval;
			m_IsLimited= true;
		}
		void SetStep(double val){m_StepSize=val;}

	protected:
		std::string m_Name;
		double m_Value;
		double m_ValueError;
		bool m_IsFixed;
		bool m_IsLimited;
		double m_MinValue;
		double m_MaxValue;
		double m_StepSize;

};//close FitPar class


class FitPars 
{
	public:
		FitPars(){
		
		};		
		virtual ~FitPars() {
			Clear();	
		};

	public: 

		struct MatchName {
			MatchName(const std::string& name) : m_name(name) {}
 			bool operator()(const FitPar& obj) const {
   			return obj.GetParName() == m_name;
 			}
 			private:
   			const std::string& m_name;
		};
	
		int GetNPars() const {return (int)(m_FitPars.size());}
		
		void AddPar(FitPar& par){
			m_FitPars.push_back(par);
		}

		FitPar* GetPar(int index) {
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return 0;
			return &m_FitPars[index];
		}

		FitPar* FindPar(int& index,std::string name){
			if(m_FitPars.empty()) return 0;
			std::vector<FitPar>::iterator it = find_if(m_FitPars.begin(), m_FitPars.end(), MatchName(name));
			if(it==m_FitPars.end()) return 0;
			size_t pos= it-m_FitPars.begin();
			index= pos;
			return GetPar(index);
		}

		int SetParValue(int index,double value){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].SetValue(value);
			return 0;
		}
		int SetParValue(std::string name,double value){
			int index= -1;
			FitPar* par= FindPar(index,name);
			if(!par) return -1;
			m_FitPars[index].SetValue(value);
			return 0;
		}
		int FixPar(int index,double value){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].SetValue(value);
			m_FitPars[index].Fix();
			return 0;
		}
		int FixPar(std::string name,double value){
			int index= -1;
			FitPar* par= FindPar(index,name);
			if(!par) return -1;
			m_FitPars[index].SetValue(value);
			m_FitPars[index].Fix();
			return 0;
		}
		int FreePar(int index){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].Free();
			return 0;
		}
		int FreePar(std::string name){
			int index= -1;
			FitPar* par= FindPar(index,name);
			if(!par) return -1;
			m_FitPars[index].Free();
			return 0;
		}
		int LimitPar(int index,double xmin,double xmax){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].SetLimits(xmin,xmax);
			return 0;
		}
		int LimitPar(std::string name,double xmin,double xmax){
			int index= -1;
			FitPar* par= FindPar(index,name);	
			if(!par) return -1;
			m_FitPars[index].SetLimits(xmin,xmax);
			return 0;
		}
	private:
		void Clear(){
			m_FitPars.clear();	
		}
	protected:
		std::vector<FitPar> m_FitPars;
	
};//close FitPars()



//===================================================
//==        SPECTRUM PARS
//===================================================
class SpectrumPars 
{
	public:
		virtual ~SpectrumPars() {};

	public:
		//Pure virtual
		virtual double GetIntegral(double xmin,double xmax) = 0;
		virtual SpectrumPars* Clone() const = 0;

		//Standard fcn
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return pars.GetNPars();}
		static int GetParNumber() {return 0;}
		virtual FitPar* GetPar(int index){return pars.GetPar(index);}
		virtual FitPars GetPars() const {return pars;}
		virtual int SetParValue(int index,double value){return pars.SetParValue(index,value);}
		virtual int FixPar(int index,double val) {
			if(SetParValue(index,val)<0) return -1;
			return pars.FixPar(index,val);
		}
		virtual int FreePar(int index) {return pars.FreePar(index);}
		virtual int LimitPar(int index,double xmin,double xmax) {return pars.LimitPar(index,xmin,xmax);}


		virtual int SetParValue(std::string name,double value){return pars.SetParValue(name,value);}
		virtual int FixPar(std::string name,double val) {
			if(SetParValue(name,val)<0) return -1;
			return pars.FixPar(name,val);
		}
		virtual int FreePar(std::string name) {return pars.FreePar(name);}
		virtual int LimitPar(std::string name,double xmin,double xmax) {return pars.LimitPar(name,xmin,xmax);}

	protected:
		FitPars pars;			
		int model;
		static int nPars;
		
};//close SpectrumPars class



enum SpectrumModel {
	eFlat= 0,
	ePowerLaw= 1,
	ePowerLawWithCutoff= 2,
	eTwoBrokenPowerLaws= 3,
	eBrokenPowerLaws= 4,
	eSmoothBrokenPowerLaws= 5,
	ePol3=6,
	ePol6=7
};

class FlatSpectrumPars : public SpectrumPars 
{
	public:		
		FlatSpectrumPars(FitPar normPar){ 
			model= eFlat;	
			nPars= 1;
			pars.AddPar(normPar);
		};
		virtual ~FlatSpectrumPars() {};

	public:
		//Overridden method
		static int GetParNumber() {return 1;}
		virtual double GetIntegral(double xmin,double xmax);
		SpectrumPars* Clone() const { 
			return new FlatSpectrumPars(*this); 
		}

};//close FlatSpectrumPars class


class PowerLawPars : public SpectrumPars {
	
	public:
		PowerLawPars(FitPar normPar,FitPar gammaPar) 
		{
			model= ePowerLaw;	
			nPars= 2;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
		};
		
		virtual ~PowerLawPars() {};

	public:
		//Overridden method
		static int GetParNumber() {return 2;}
		virtual double GetIntegral(double xmin,double xmax);
		SpectrumPars* Clone() const { 
			return new PowerLawPars(*this); 
		}
		
};//close PowerLawPars class



class BrokenPowerLawsPars : public SpectrumPars {
	
	public:
		BrokenPowerLawsPars(FitPar normPar,FitPar gammaPar,FitPar gamma2Par,FitPar gamma3Par,FitPar breakPar,FitPar cutoffPar) 
		{
			model= eBrokenPowerLaws;	
			nPars= 6;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
			pars.AddPar(gamma2Par);
			pars.AddPar(gamma3Par);
			pars.AddPar(breakPar);
			pars.AddPar(cutoffPar);
		};
		
		virtual ~BrokenPowerLawsPars() {};

	public:
		//Overridden method
		static int GetParNumber() {return 6;}
		virtual double GetIntegral(double xmin,double xmax);
		SpectrumPars* Clone() const { 
			return new BrokenPowerLawsPars(*this); 
		}
		
};//close BrokenPowerLawsPars class


class TwoBrokenPowerLawsPars : public SpectrumPars {
	
	public:
		TwoBrokenPowerLawsPars(FitPar normPar,FitPar gammaPar,FitPar gamma2Par,FitPar breakPar) 
		{
			model= eTwoBrokenPowerLaws;	
			nPars= 4;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
			pars.AddPar(gamma2Par);
			pars.AddPar(breakPar);
		};
		
		virtual ~TwoBrokenPowerLawsPars() {};

	public:
		//Overridden method
		static int GetParNumber() {return 6;}
		virtual double GetIntegral(double xmin,double xmax);
		SpectrumPars* Clone() const { 
			return new TwoBrokenPowerLawsPars(*this); 
		}
		
};//close TwoBrokenPowerLawsPars class



class SmoothCutoffPowerLaws : public SpectrumPars 
{	
	public:
		SmoothCutoffPowerLaws(FitPar normPar,FitPar gammaPar,FitPar gamma2Par,FitPar breakPar,FitPar cutoffPar,FitPar smoothCutoffPar) 
		{
			model= eSmoothBrokenPowerLaws;	
			nPars= 6;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
			pars.AddPar(gamma2Par);
			pars.AddPar(breakPar);
			pars.AddPar(cutoffPar);
			pars.AddPar(smoothCutoffPar);
		};
		
		virtual ~SmoothCutoffPowerLaws() {};

	public:
		//Overridden method
		static int GetParNumber() {return 6;}
		virtual double GetIntegral(double xmin,double xmax) {
			return 1;
		}
		SpectrumPars* Clone() const { 
			return new SmoothCutoffPowerLaws(*this); 
		}
		
};//close SmoothCutoffPowerLaws

class PowerLawWithCutoff : public SpectrumPars 
{	
	public:
		PowerLawWithCutoff(FitPar normPar,FitPar gammaPar,FitPar cutoffPar,FitPar smoothCutoffPar) 
		{
			model= ePowerLawWithCutoff;	
			nPars= 4;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
			pars.AddPar(cutoffPar);
			pars.AddPar(smoothCutoffPar);
		};
		
		virtual ~PowerLawWithCutoff() {};

	public:
		//Overridden method
		static int GetParNumber() {return 4;}
		virtual double GetIntegral(double xmin,double xmax) {
			return 1;
		}
		SpectrumPars* Clone() const { 
			return new PowerLawWithCutoff(*this); 
		}
		
};//close PowerLawWithCutoff


class Pol3SpectrumPars : public SpectrumPars {
	
	public:
		Pol3SpectrumPars(FitPar p0,FitPar p1,FitPar p2,FitPar p3) 
		{
			model= ePol3;	
			nPars= 4;
			pars.AddPar(p0);
			pars.AddPar(p1);
			pars.AddPar(p2);
			pars.AddPar(p3);	
		};
		
		virtual ~Pol3SpectrumPars() {};

	public:
		//Overridden method
		static int GetParNumber() {return 4;}
		virtual double GetIntegral(double xmin,double xmax) {
			return 1;
		}
		SpectrumPars* Clone() const { 
			return new Pol3SpectrumPars(*this); 
		}
		
};//close Pol3SpectrumPars class


class Pol6SpectrumPars : public SpectrumPars {
	
	public:
		Pol6SpectrumPars(FitPar p0,FitPar p1,FitPar p2,FitPar p3,FitPar p4,FitPar p5,FitPar p6) 
		{
			model= ePol6;	
			nPars= 7;
			pars.AddPar(p0);
			pars.AddPar(p1);
			pars.AddPar(p2);
			pars.AddPar(p3);
			pars.AddPar(p4);
			pars.AddPar(p5);
			pars.AddPar(p6);		
		};
		
		virtual ~Pol6SpectrumPars() {};

	public:
		//Overridden method
		static int GetParNumber() {return 4;}
		virtual double GetIntegral(double xmin,double xmax) {
			return 1;
		}
		SpectrumPars* Clone() const { 
			return new Pol6SpectrumPars(*this); 
		}
		
};//close Pol6SpectrumPars class

//===================================================
//==        SPECTRUM UTILS
//===================================================
class SpectrumUtils 
{ 
	public:

		/** 
		\brief Class constructor
 		*/
		SpectrumUtils();
		/** 
		\brief Class destructor
 		*/
		virtual ~SpectrumUtils();
		

	public:
		//==========================================
		//==      HISTO MANIPULATION FUNCTIONS
		//==========================================
		/** 
		\brief Return histogram with counts fluctuated
 		*/
		static TH1D* GetFluctuatedSpectrum(TH1D*);
		/** 
		\brief Return histogram folded with response matrix
 		*/
		static TH1D* GetFoldedSpectrum(TH1D* trueSpectrum,TH2D* responseMatrix);
		/** 
		\brief Copy histogram content
 		*/
		static int CopySpectrumContent(TH1D* h1,TH1D* h2);
		/** 
		\brief Return histogram from model function
 		*/
		static int GetModelSpectrum(TH1D& spectrum,TF1* spectrumModel,bool integrateBins=true);		
		/** 
		\brief Check histogram binnings
 		*/
		static bool CheckSpectrumBinnings(const TAxis* a1, const TAxis* a2);
		/** 
		\brief Compute spectrum model function
 		*/
		static TF1* ComputeSpectrumModel(SpectrumPars& pars,double xmin,double xmax,int npts=1000);
		/** 
		\brief Compute efficiency model function
 		*/
		static TF1* ComputeEfficiencyModel(EfficiencyPars& pars,double xmin,double xmax,int npts=1000);
		/** 
		\brief Compute reso model function
 		*/
		static TF1* ComputeResoModel(ResoPars& resoPars,double xmin,double xmax,int npts=1000);
		/** 
		\brief Compute bias model function
 		*/
		static TF1* ComputeBiasModel(BiasPars& biasPars,double xmin,double xmax,int npts=1000);

		/** 
		\brief Compute response model function
 		*/
		static TF2* ComputeResponseModel(SpectrumPars& spectrumPars,BiasPars& biasPars,ResoPars& resoPars,EfficiencyPars& effPars,double xmin,double xmax,double ymin,double ymax,int npts=1000);
		/** 
		\brief Compute response model matrix
 		*/
		static TH2D* ComputeParametricResponse(SpectrumPars& spectrumPars,BiasPars& biasPars,ResoPars& resoPars,EfficiencyPars& effPars,std::vector<double>& TrueBins, std::vector<double>& RecBins,int npts=1000);
		/** 
		\brief Build response matrix
 		*/
		static TH2D* BuildResponseMatrix(SpectrumPars& spectrumPars,BiasPars& biasPars,ResoPars& resoPars,EfficiencyPars& effPars,std::vector<double>& TrueBins, std::vector<double>& RecBins,int npts=1000);

		//==========================================
		//==      MATH FUNCTIONS
		//==========================================
		// - Spectrum models
		/** 
		\brief Power law spectrum model
 		*/
		static double PowerLawSpectrum(double* x, double* par);	
		/** 
		\brief Power law spectrum with cutoff model
 		*/
		static double PowerLawWithCutoffSpectrum(double* x, double* par);	
		/** 
		\brief Broken power-law spectrum model
 		*/
		static double BrokenPowerLawSpectrum(double* x, double* par);
		/** 
		\brief 2 broken power-law spectrum model
 		*/
		static double TwoBrokenPowerLawSpectrum(double* x, double* par);
		/** 
		\brief Smooth cutoff power-law spectrum model
 		*/
		static double SmoothCutoffPowerLawSpectrum(double* x, double* par);	
		/** 
		\brief 3rd degree polynomial spectrum model
 		*/
		static double Pol3Spectrum(double* x, double* par);
		/** 
		\brief 6th degree polynomial spectrum model
 		*/
		static double Pol6Spectrum(double* x, double* par);	

		/** 
		\brief Compute power law integral
 		*/
		static double GetPowerLawIntegral(double gamma,double lgS_min, double lgS_max);
		/** 
		\brief Compute broken power law integral
 		*/
		static double GetBrokenPowerLawIntegral(double gamma1,double gamma2,double gamma3,double lgS_break,double lgS_cutoff,double lgS_min, double lgS_max);
		/** 
		\brief Compute 2 broken power law integral
 		*/
		static double GetTwoBrokenPowerLawIntegral(double gamma1,double gamma2,double lgS_break,double lgS_min, double lgS_max);
	
		
		// - Response models
		/** 
		\brief Compute flux exp resolution model
 		*/
		static double ExpResolutionModel(double* x, double* par);
		/** 
		\brief Compute flux exp step resolution model
 		*/
		static double ExpStepResolutionModel(double* x, double* par);
		/** 
		\brief Compute flux bias exp model
 		*/
		static double ExpBiasModel(double* x, double* par);
		/** 
		\brief Compute flux bias step+exp model
 		*/
		static double StepExpBiasModel(double* x, double* par);

		/** 
		\brief Compute detection efficiency model
 		*/
		static double EfficiencyModel(double* x, double* par);
		/** 
		\brief Compute 2D response model
 		*/
		static double ResponseModel(double* x, double* par);
		/** 
		\brief Compute response model 1D
 		*/
		static double ResponseModel1D(double* x, double* par);


		// - Source counts resolution bias model
		/** 
		\brief Compute max angular size model fcn
 		*/
		static double PhiMaxModel(double* x,double* pars);
		/** 
		\brief Compute min angular size model fcn
 		*/
		static double PhiMinModel(double* x,double* pars);
		/** 
		\brief Compute median angular size model fcn
 		*/
		static double PhiMedianModel(double* x,double* pars);
		/** 
		\brief Compute source counts reso bias model
 		*/
		static double SourceCountsResoBiasModel(double* x,double* pars);
		/** 
		\brief Compute source counts reso bias corr factor model
 		*/
		static double SourceCountsResoBiasCorrFactor(double* x,double* pars);


};//close class


#endif
