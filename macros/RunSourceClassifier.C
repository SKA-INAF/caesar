#include <TROOT.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <Math/QuantFuncMathCore.h>
//#include <TMultiLayerPerceptron.h>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
//using namespace TMVA;

#include <StatsUtils.h>
//#include <Logger.h>
using namespace Caesar;

#include <vector>
using namespace std;


struct NNOptions {
	
	NNOptions(){
		SetDefaults();
	}

	void SetDefaults(){
		outputFileName= "NNClassifierOutput.root";
		nIters= 5000;
		//nIters= 10;	
		//hiddenLayers= "N,N-1";
		hiddenLayers= "N,N,N";
		neuronType= "sigmoid";
		trainingMethod= "BFGS";
		varTransform= "N";
		useRegulator= false;
		addSourceSNRToNNVars= false;
		NNCut= 0.5;
	}
	
	std::string GetStr(){
		std::stringstream ss;
		ss<<"H:!V:";
		ss<<"NeuronType="<<neuronType<<":";
		ss<<"VarTransform="<<varTransform<<":";
		ss<<"NCycles="<<nIters<<":";
		ss<<"HiddenLayers="<<hiddenLayers<<":";
		ss<<"TrainingMethod="<<trainingMethod<<":";
		if(useRegulator) ss<<"UseRegulator";
		else ss<<"!UseRegulator";
		return ss.str();
	}

	int nIters;
	std::string hiddenLayers;
	std::string neuronType;
	std::string trainingMethod;
	std::string varTransform;
	bool useRegulator;
	bool addSourceSNRToNNVars;
	std::string outputFileName;
	double NNCut;

};//close NNOptions

struct CutOptions {

	CutOptions(){
		SetDefaults();
	}
	
	void SetDefaults(){
		addSourceSNRToCutVars= false;
		outputFileName= "CutsClassifierOutput.root";
	}

	bool addSourceSNRToCutVars;
	std::string outputFileName;

};


struct MacroOptions {

	MacroOptions(){
		SetDefaults();
	}

	void SetDefaults(){
		minSourceSNR= -50000;
		maxSourceSNR= 50000;
		minThetaDiff= -90;
		maxThetaDiff= 90;
		minEccentricityRatio= 0;
		maxEccentricityRatio= 10;
		minSourceToBeamRatio= 0;
		maxSourceToBeamRatio= 50;
		runCutClassifier= false;
		runNNClassifier= true;
	}

	bool runNNClassifier;
	bool runCutClassifier;
	NNOptions nnOptions;
	CutOptions cutOptions;
	float minSourceSNR;
	float maxSourceSNR;
	float minEccentricityRatio;
	float maxEccentricityRatio;
	float minSourceToBeamRatio; 
	float maxSourceToBeamRatio;
	float minThetaDiff;
	float maxThetaDiff;
};
MacroOptions opt;

//Classifier
//TMVA::Factory* factory= 0;
//TMVA::DataLoader* dataloader= 0;
//TMVA::Reader* reader = 0;

//Options
bool normalizeData= false;
double dataNormMin= -1;
double dataNormMax= 1;


//Data
TFile* outputFile_cuts= 0;
TFile* outputFile_nn= 0;
TTree* dataTree= 0;
TFile* inputFile= 0;
std::string inputFileName= "";
TDirectory* currentDir= 0;
TEventList* trainingDataList;
TEventList* testDataList;
TTree* NNDataTree;
TTree* NNDataTree_train;
TTree* NNDataTree_test;
int nVars= 4;
std::vector<TString> varNames {"#delta#theta (deg)","E_{fit}/E_{beam}","A_{source}/A_{beam}","SNR"};
float thetaDiff;
float eccentricityRatio;
float sourceToBeamRatio;
float sourceSNR;
int isTrueSource;
float minSourceSNR= -50000;//5
float maxSourceSNR= 50000;
std::vector<Caesar::BoxStats<double>> dataVars_stats;
std::vector<TH1D*> dataVarHistos_signal;
std::vector<TH1D*> dataVarHistos_bkg;
std::vector<TGraph*> dataVarCorrGraphs_signal;
std::vector<TGraph*> dataVarCorrGraphs_bkg;

//NN vars
double NNCut= 0.5;
double NEventFractionForTraining= 0.5;
TH1D* NNOutputHisto_signal= 0;
TH1D* NNOutputHisto_bkg= 0;

//Functions
int Init();
int ReadData(std::string);
void DrawData();
//int TrainNN(int maxIter=1000);
int RunCutClassifier();
int RunNNClassifier();
int EvaluateNNClassifier();
int RunClassifiers();
void Save();

int RunSourceClassifier(std::string fileName,MacroOptions macroOptions=MacroOptions())
{
	currentDir= gDirectory;
	opt= macroOptions;	
	inputFileName= fileName;

	//## Init
	Init();

	//## Read data
	ReadData(fileName);

	//## Train classifier
	RunClassifiers();

	//## Draw data
	DrawData();
	
	//## Save data
	//Save();

	return 0;

}//close macro

void Save()
{
	//Close file
	if(outputFile_cuts && outputFile_cuts->IsOpen()){
		outputFile_cuts->Close();
	}
	if(outputFile_nn && outputFile_nn->IsOpen()){
		outputFile_nn->Close();
	}

}//close Save()

int Init()
{
	//Init MVA 
	TMVA::Tools::Instance();

	//Initialize reader		
	//reader = new TMVA::Reader( "!Color:!Silent" );

	//Open output file
	outputFile_cuts= new TFile((opt.cutOptions).outputFileName.c_str(),"RECREATE");
	outputFile_nn= new TFile((opt.nnOptions).outputFileName.c_str(),"RECREATE");

	//Initialize event data lists
	trainingDataList= new TEventList;
	testDataList= new TEventList;

	//Initialize data NN tree
	NNDataTree= new TTree("NNDataTree","NNDataTree");
	NNDataTree->Branch("thetaDiff",&thetaDiff,"thetaDiff/D");
	NNDataTree->Branch("eccentricityRatio",&eccentricityRatio,"eccentricityRatio/D");
	NNDataTree->Branch("sourceToBeamRatio",&sourceToBeamRatio,"sourceToBeamRatio/D");
	if((opt.nnOptions).addSourceSNRToNNVars) NNDataTree->Branch("sourceSNR",&sourceSNR,"sourceSNR/D");
	NNDataTree->Branch("isTrueSource",&isTrueSource,"isTrueSource/I");

	NNDataTree_train= new TTree("NNDataTree_train","NNDataTree_train");
	NNDataTree_train->Branch("thetaDiff",&thetaDiff,"thetaDiff/D");
	NNDataTree_train->Branch("eccentricityRatio",&eccentricityRatio,"eccentricityRatio/D");
	NNDataTree_train->Branch("sourceToBeamRatio",&sourceToBeamRatio,"sourceToBeamRatio/D");
	if((opt.nnOptions).addSourceSNRToNNVars) NNDataTree_train->Branch("sourceSNR",&sourceSNR,"sourceSNR/D");
	NNDataTree_train->Branch("isTrueSource",&isTrueSource,"isTrueSource/I");

	NNDataTree_test= new TTree("NNDataTree_test","NNDataTree_test");
	NNDataTree_test->Branch("thetaDiff",&thetaDiff,"thetaDiff/D");
	NNDataTree_test->Branch("eccentricityRatio",&eccentricityRatio,"eccentricityRatio/D");
	NNDataTree_test->Branch("sourceToBeamRatio",&sourceToBeamRatio,"sourceToBeamRatio/D");
	if((opt.nnOptions).addSourceSNRToNNVars) NNDataTree_test->Branch("sourceSNR",&sourceSNR,"sourceSNR/D");
	NNDataTree_test->Branch("isTrueSource",&isTrueSource,"isTrueSource/I");
	
	//Initialize NN out histos
	//NNOutputHisto_signal= new TH1D("NNOutputHisto_signal","NNOutputHisto_signal",100,-3,3);
	//NNOutputHisto_bkg= new TH1D("NNOutputHisto_bkg","NNOutputHisto_bkg",100,-3,3);

	return 0;

}//close Init()

int ReadData(std::string fileName)
{
	//Open file
	inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()){
		//ERROR_LOG("Failed to open file "<<fileName<<"!");
		cerr<<"Failed to open file "<<fileName<<"!"<<endl;
		return -1;	
	}

	//Read tree
	dataTree= (TTree*)inputFile->Get("FitEllipseInfo");
	if(!dataTree){
		//ERROR_LOG("Cannot read tree from file "<<fileName<<"!");
		cerr<<"Cannot read tree from file "<<fileName<<"!"<<endl;
		return -1;
	}
		
	dataTree->SetBranchAddress("isTrueSource",&isTrueSource);
	dataTree->SetBranchAddress("thetaDiff",&thetaDiff);
	dataTree->SetBranchAddress("eccentricityRatio",&eccentricityRatio);
	dataTree->SetBranchAddress("sourceToBeamRatio",&sourceToBeamRatio);
	dataTree->SetBranchAddress("sourceSNR",&sourceSNR);


	//Read and select data
	int nBranches= NNDataTree->GetNbranches();
	nVars= nBranches-1;
	std::vector<std::vector<double>> dataVarList;
	std::vector<std::vector<double>> events;
	for(int j=0;j<nVars;j++){
		dataVarList.push_back( std::vector<double>() );
	}
	int nSelEvents= 0;
	
	//INFO_LOG("Reading #"<<dataTree->GetEntries()<<" data entries present...");
	cout<<"Reading #"<<dataTree->GetEntries()<<" data entries present..."<<endl;
		
	for(int i=0;i<dataTree->GetEntries();i++){
		dataTree->GetEntry(i);
		
		//Check source SNR range
		if((opt.nnOptions).addSourceSNRToNNVars && (sourceSNR<opt.minSourceSNR || sourceSNR>opt.maxSourceSNR || TMath::IsNaN(sourceSNR) || !std::isfinite(sourceSNR)) ){
			cout<<"Skip sourceSNR (SNR="<<sourceSNR<<") as outside selected range ("<<minSourceSNR<<","<<maxSourceSNR<<")"<<endl;
			continue;
		}

		//Check eccentricity ratio
		if(eccentricityRatio<opt.minEccentricityRatio || eccentricityRatio>opt.maxEccentricityRatio || TMath::IsNaN(eccentricityRatio) || !std::isfinite(eccentricityRatio) ){
			continue;
		}

		//Check theta diff
		if(thetaDiff<opt.minThetaDiff || thetaDiff>opt.maxThetaDiff || TMath::IsNaN(thetaDiff) || !std::isfinite(thetaDiff) ){
			continue;
		}

		//Check source to beam area ratio
		if(sourceToBeamRatio<opt.minSourceToBeamRatio || sourceToBeamRatio>opt.maxSourceToBeamRatio || TMath::IsNaN(sourceToBeamRatio) || !std::isfinite(sourceToBeamRatio) ){
			continue;
		}

		//Fill selected events
		dataVarList[0].push_back(thetaDiff);
		dataVarList[1].push_back(eccentricityRatio);
		dataVarList[2].push_back(sourceToBeamRatio);
		if((opt.nnOptions).addSourceSNRToNNVars) dataVarList[3].push_back(sourceSNR);
		
		events.push_back( std::vector<double>() );
		events[nSelEvents].push_back(thetaDiff);
		events[nSelEvents].push_back(eccentricityRatio);
		events[nSelEvents].push_back(sourceToBeamRatio);
		if((opt.nnOptions).addSourceSNRToNNVars) events[nSelEvents].push_back(sourceSNR);
		events[nSelEvents].push_back(isTrueSource);
		nSelEvents++;
		
	}//end loop events

	//Compute data var stats
	std::vector<double> dataVars_min;
	std::vector<double> dataVars_max;

	for(int j=0;j<nVars;j++){	
		BoxStats<double> stats= StatsUtils::ComputeBoxStats(dataVarList[j],false);
		double wmin= stats.minVal;
		double wmax= stats.maxVal;
		dataVars_min.push_back(wmin);
		dataVars_max.push_back(wmax);
		//INFO_LOG("Data var "<<j+1<<": original range ("<<wmin<<","<<wmax<<"), norm range ("<<dataNormMin<<","<<dataNormMax<<")");
		cout<<"Data var "<<j+1<<": original range ("<<wmin<<","<<wmax<<"), norm range ("<<dataNormMin<<","<<dataNormMax<<")"<<endl;

		dataVars_stats.push_back(stats);
	}//end loop vars
	
	
	//Normalize data
	if(normalizeData){
		//INFO_LOG("Normalizing input data (#"<<events.size()<<" data events) ...");
		cout<<"Normalizing input data (#"<<events.size()<<" data events) ..."<<endl;
		for(size_t i=0;i<events.size();i++){//loop on events
			for(size_t j=0;j<events[i].size()-1;j++){//loop on data vars
				double wmin= dataVars_min[j];
				double wmax= dataVars_max[j];
				double w= events[i][j];
				double w_norm= dataNormMin + (dataNormMax-dataNormMin)*(w-wmin)/(wmax-wmin);
				events[i][j]= w_norm;
			}//end loop data vars
		}//end loop events
	}//close if	

	//Switch to original dir	
	currentDir->cd();

	//Initialize data histos
	TH1D* histo= 0;
	int nBins= 200;
	for(int i=0;i<nVars;i++){
		double wmin= dataVars_min[i];
		double wmax= dataVars_max[i];
		if(normalizeData){
			wmin= dataNormMin;
			wmax= dataNormMax;
		}

		TString histoName= Form("dataHisto_var%d_signal",i+1);
		histo= new TH1D(histoName,histoName,nBins,wmin,wmax);
		dataVarHistos_signal.push_back(histo);

		histoName= Form("dataHisto_var%d_bkg",i+1);
		histo= new TH1D(histoName,histoName,nBins,wmin,wmax);
		dataVarHistos_bkg.push_back(histo);
	}//end loop vars


	//Read data and fill event lists
	//INFO_LOG("Setting NN train/test data (#"<<events.size()<<" data entries present) ...");
	cout<<"Setting NN train/test data (#"<<events.size()<<" data entries present) ..."<<endl;
	for(size_t i=0;i<events.size();i++){
		thetaDiff= events[i][0];
		eccentricityRatio= events[i][1];
		sourceToBeamRatio= events[i][2];
		if((opt.nnOptions).addSourceSNRToNNVars){
			sourceSNR= events[i][3];	
			isTrueSource= static_cast<int>(events[i][4]);
		}
		else{
			isTrueSource= static_cast<int>(events[i][3]);
		}

		double rand= gRandom->Uniform(0,1);
		if(rand<=NEventFractionForTraining){//fill training list
			trainingDataList->Enter(i);	
			NNDataTree_train->Fill();
		}
		else{//fill crossval list
			testDataList->Enter(i);	
			NNDataTree_test->Fill();
		}
		NNDataTree->Fill();

		//Fill histos
		if(isTrueSource){
			for(int j=0;j<nVars;j++) dataVarHistos_signal[j]->Fill(events[i][j]);
		}
		else{
			for(int j=0;j<nVars;j++) dataVarHistos_bkg[j]->Fill(events[i][j]);
		}

	}//end loop sel events

	int N_train= trainingDataList->GetN();
	int N_test= testDataList->GetN();

	cout<<"INFO: Train/test events="<<N_train<<"/"<<N_test<<endl;
	
	return 0;

}//close ReadData()

int RunClassifiers()
{
	//Run cuts classifier
	if(opt.runCutClassifier){
		RunCutClassifier();
	}

	//Run NN classifier
	if(opt.runNNClassifier){
		//TrainNN(nnIters);
		RunNNClassifier();	
		EvaluateNNClassifier();
	}

	return 0;

}//close RunClassifiers()

int EvaluateNNClassifier()
{
	//######################################
	//##   EVALUATE CLASSIFIER
	//######################################
	cout<<"INFO: Initializing the MVA reader..."<<endl;
  //if(!reader) reader = new TMVA::Reader( "!Color:!Silent" );
	TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );
	reader->AddVariable("thetaDiff",&thetaDiff);
	reader->AddVariable("eccentricityRatio",&eccentricityRatio);
	reader->AddVariable("sourceToBeamRatio",&sourceToBeamRatio);
	if((opt.nnOptions).addSourceSNRToNNVars) reader->AddVariable("sourceSNR",&sourceSNR);

	TString weightFile= "dataset/weights/SourceNNClassification_MLP.weights.xml";
	cout<<"INFO: Booking MVA method using weight file "<<weightFile.Data()<<" ..."<<endl;
	reader->BookMVA("MLP",weightFile);
	

	NNOutputHisto_signal= new TH1D("NNOutputHisto_signal","NNOutputHisto_signal",100,-3,3);
	NNOutputHisto_bkg= new TH1D("NNOutputHisto_bkg","NNOutputHisto_bkg",100,-3,3);

	int NSig_true= 0;
	int NSig= 0;
	int NSig_right= 0;
	int NSig_wrong= 0;
	int NBkg= 0;
	cout<<"INFO: #"<<dataTree->GetEntries()<<" events will be read to test NN output..."<<endl;
	for (int i=0;i<dataTree->GetEntries();i++) {
		dataTree->GetEntry(i);

		//Check source SNR range
		if((opt.nnOptions).addSourceSNRToNNVars && (sourceSNR<opt.minSourceSNR || sourceSNR>opt.maxSourceSNR || TMath::IsNaN(sourceSNR) || !std::isfinite(sourceSNR)) ){
			cout<<"INFO: Skip sourceSNR (SNR="<<sourceSNR<<") as outside selected range ("<<minSourceSNR<<","<<maxSourceSNR<<")"<<endl;
			continue;
		}

		//Check eccentricity ratio
		if(eccentricityRatio<opt.minEccentricityRatio || eccentricityRatio>opt.maxEccentricityRatio || TMath::IsNaN(eccentricityRatio) || !std::isfinite(eccentricityRatio) ){
			continue;
		}

		//Check theta diff
		if(thetaDiff<opt.minThetaDiff || thetaDiff>opt.maxThetaDiff || TMath::IsNaN(thetaDiff) || !std::isfinite(thetaDiff) ){
			continue;
		}

		//Check source to beam area ratio
		if(sourceToBeamRatio<opt.minSourceToBeamRatio || sourceToBeamRatio>opt.maxSourceToBeamRatio || TMath::IsNaN(sourceToBeamRatio) || !std::isfinite(sourceToBeamRatio) ){
			continue;
		}

		if(i%10000==0) cout<<"INFO: Reading "<<i+1<<"/"<<dataTree->GetEntries()<<" event..."<<endl;

		double NNOut= 0;
		try{
			NNOut = reader->EvaluateMVA("MLP");
		}
		catch(...){
			cerr<<"WARN: Invalid NN evaluation on event "<<i+1<<"!"<<endl;
		}

		//Fill NN out histo
		if(isTrueSource==1){
			NSig_true++;
			NNOutputHisto_signal->Fill(NNOut);
		}
		else{
			NNOutputHisto_bkg->Fill(NNOut);
		}

		//Update efficiency
		if(NNOut>=(opt.nnOptions).NNCut){//identified as signal
			NSig++;
			if(isTrueSource==1) NSig_right++;//correct identification
			else NSig_wrong++;//wrong identification
		}
		else{
			NBkg++;
		}

	}//end loop events

	//Compute efficiency & purity
	double Eff= (double)NSig_right/(double)NSig_true;
	double Purity= (double)NSig_right/(double)NSig;
	cout<<"INFO: NSig_true="<<NSig_true<<", NSig_right="<<NSig_right<<", NSig="<<NSig<<", Efficiency="<<Eff<<", Purity="<<Purity<<endl;

	//Draw NN output histos
	TCanvas* NNOutPlot= new TCanvas("NNOutPlot","NNOutPlot");
	NNOutPlot->cd();

	TH2D* NNOutPlotBkg= new TH2D("NNOutPlotBkg","",100,-2,2,100,0,1);
	NNOutPlotBkg->SetXTitle("NNOut");
	NNOutPlotBkg->SetYTitle("entries");
	NNOutPlotBkg->SetStats(0);
	NNOutPlotBkg->Draw();

	NNOutputHisto_signal->SetFillStyle(3001);
	NNOutputHisto_signal->SetLineColor(kRed);
	NNOutputHisto_signal->DrawNormalized("hist same");

	NNOutputHisto_bkg->SetLineColor(kBlack);
	NNOutputHisto_bkg->DrawNormalized("hist same");

	
	//Clear stuff
	if(reader) {
		cout<<"DEBUG: Deleting reader..."<<endl;
		delete reader;
		reader= 0;
		cout<<"DEBUG: done!"<<endl;
	}
	

	return 0;

}//close EvaluateNNClassifier()


int RunNNClassifier()
{
	//######################################
	//##   SET DATA
	//######################################
	//Initialize factory
	TMVA::Factory* factory = new TMVA::Factory( 
		"SourceNNClassification", 
		outputFile_nn,
   	"!V:!Silent:Color:DrawProgressBar:Transformations=N:AnalysisType=Classification" 
	);

	//Set data
	TMVA::DataLoader* dataloader= new TMVA::DataLoader("dataset");
	dataloader->AddVariable("thetaDiff", 'F');
	dataloader->AddVariable("eccentricityRatio", 'F');
	dataloader->AddVariable("sourceToBeamRatio", 'F');
	if( (opt.nnOptions).addSourceSNRToNNVars) dataloader->AddVariable("sourceSNR", 'F');

	TCut signalCut = "isTrueSource==1"; // how to identify signal events
	TCut bkgCut = "isTrueSource==0"; // how to identify background events
	dataloader->SetInputTrees(dataTree, signalCut, bkgCut);

	// Apply additional cuts on data
	TCut sourceSNRCut= Form("sourceSNR>=%f && sourceSNR<=%f",opt.minSourceSNR,opt.maxSourceSNR);
	TCut thetaDiffCut= Form("thetaDiff>=%f && thetaDiff<=%f",opt.minThetaDiff,opt.maxThetaDiff);
	TCut eccentricityRatioCut= Form("eccentricityRatio>=%f && eccentricityRatio<=%f",opt.minEccentricityRatio,opt.maxEccentricityRatio);
	TCut sourceToBeamRatioCut= Form("sourceToBeamRatio>=%f && sourceToBeamRatio<=%f",opt.minSourceToBeamRatio,opt.maxSourceToBeamRatio);

  TCut preselectionCut= thetaDiffCut && eccentricityRatioCut && sourceToBeamRatioCut;
	if((opt.nnOptions).addSourceSNRToNNVars) preselectionCut= thetaDiffCut && eccentricityRatioCut && sourceToBeamRatioCut && sourceSNRCut;
	dataloader->PrepareTrainingAndTestTree( preselectionCut, "SplitMode=Random:NormMode=NumEvents:!V");


	//######################################
	//##   ADD CLASSIFIERS
	//######################################
	//Add NN classifier
	std::string nnOptions= (opt.nnOptions).GetStr();
	factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", nnOptions.c_str() );
	
	//######################################
	//##   TRAIN CLASSIFIERS
	//######################################
	//Train classifier
	factory->TrainAllMethods();
   
	//######################################
	//##   TEST CLASSIFIERS
	//######################################
	// Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

	// Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

	//Close file
	if(outputFile_nn && outputFile_nn->IsOpen()) outputFile_nn->Close();

	
	//Clear stuff
	if(factory) {
		delete factory;
		factory= 0;
	}
	if(dataloader){
 		delete dataloader;
		dataloader= 0;
	}
		
	return 0;

}//close RunNNClassifier()

int RunCutClassifier()
{
	//######################################
	//##   SET DATA
	//######################################
	//Initialize factory
	TMVA::Factory* factory = new TMVA::Factory( 
		"SourceCutsClassification", 
		outputFile_cuts,
   	"!V:!Silent:Color:DrawProgressBar:Transformations=N;I;D;P;G,D:AnalysisType=Classification" 
	);

	//Set data
	TMVA::DataLoader* dataloader= new TMVA::DataLoader("dataset");
	dataloader->AddVariable("thetaDiff", 'F');
	dataloader->AddVariable("eccentricityRatio", 'F');
	dataloader->AddVariable("sourceToBeamRatio", 'F');
	if( (opt.cutOptions).addSourceSNRToCutVars ) dataloader->AddVariable("sourceSNR", 'F');

	TCut signalCut = "isTrueSource==1"; // how to identify signal events
	TCut bkgCut = "isTrueSource==0"; // how to identify background events
	dataloader->SetInputTrees(dataTree, signalCut, bkgCut);

	// Apply additional cuts on data
	TCut sourceSNRCut= Form("sourceSNR>=%f && sourceSNR<=%f",opt.minSourceSNR,opt.maxSourceSNR);
	TCut thetaDiffCut= Form("thetaDiff>=%f && thetaDiff<=%f",opt.minThetaDiff,opt.maxThetaDiff);
	TCut eccentricityRatioCut= Form("eccentricityRatio>=%f && eccentricityRatio<=%f",opt.minEccentricityRatio,opt.maxEccentricityRatio);
	TCut sourceToBeamRatioCut= Form("sourceToBeamRatio>=%f && sourceToBeamRatio<=%f",opt.minSourceToBeamRatio,opt.maxSourceToBeamRatio);

  TCut preselectionCut= thetaDiffCut && eccentricityRatioCut && sourceToBeamRatioCut;
	if((opt.nnOptions).addSourceSNRToNNVars) preselectionCut= thetaDiffCut && eccentricityRatioCut && sourceToBeamRatioCut && sourceSNRCut;

	dataloader->PrepareTrainingAndTestTree( preselectionCut, "SplitMode=Random:NormMode=NumEvents:!V");
	

	//######################################
	//##   BOOK CLASSIFIER
	//######################################
	//Select classifier
	//factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
  //                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
	factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=GA:EffSel");
	//factory->BookMethod(dataloader,TMVA::Types::kCuts, "CutsGA",
  //                        "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

	//######################################
	//##   TRAIN CLASSIFIERS
	//######################################
	//Train classifier
	factory->TrainAllMethods();
   
	//######################################
	//##   TEST CLASSIFIERS
	//######################################
	// Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

	// Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

	//Close file
	if(outputFile_cuts && outputFile_cuts->IsOpen()) outputFile_cuts->Close();

	//Clear stuff
	delete factory;
  delete dataloader;

	return 0;

}//close RunCutClassifier()

/*
int TrainNN(int maxIter)
{
	//## Create NN
	TMultiLayerPerceptron* NN= new TMultiLayerPerceptron("thetaDiff,eccentricityRatio,sourceToBeamRatio,sourceSNR:4:4:4:isTrueSource",NNDataTree, trainingDataList, testDataList, TNeuron::kTanh);
	NN->SetLearningMethod(TMultiLayerPerceptron::kBFGS);
	//NN->LoadWeights("weights.dat");
	//NN->SetTau(0.001);
 	//NN->SetEta(0.001);
	//NN->SetReset(250);
	
	//## Train network
	//NN->Train(5,"text,graph,update=5");
	//NN->Train(maxIter,"text,graph,update=5");
	NN->Train(maxIter,"text,graph,update=5, +");

	//## Dump final weights to ascii file
	NN->DumpWeights("NNWeights.dat");


	//## Analyze output
	double NNInput[nVars];
	double NSig_right= 0;
	double NSig_wrong= 0;
	double NSig= 0;
	double NBkg= 0;
	double NSig_true= 0;
	
	for(int i=0;i<NNDataTree_test->GetEntries();i++){
		NNDataTree_test->GetEntry(i);

		NNInput[0]= thetaDiff;
		NNInput[1]= eccentricityRatio;
		NNInput[2]= sourceToBeamRatio;
		if((opt.nnOptions).addSourceSNRToNNVars) NNInput[3]= sourceSNR;
		double NNOut = NN->Evaluate(0,NNInput);

		//Fill NN out histo
		if(isTrueSource==1){
			NSig_true++;
			NNOutputHisto_signal->Fill(NNOut);
		}
		else{
			NNOutputHisto_bkg->Fill(NNOut);
		}

		//Update efficiency
		if(NNOut>=NNCut){//identified as signal
			NSig++;
			if(isTrueSource==1) NSig_right++;//correct identification
			else NSig_wrong++;//wrong identification
		}
		else{
			NBkg++;
		}
	
	}//end loop events

	//Compute efficiency & purity
	double Eff= NSig_right/NSig_true;
	double Purity= NSig_right/NSig;
	INFO_LOG("NSig_true="<<NSig_true<<", NSig_right="<<NSig_right<<", NSig="<<NSig<<", Efficiency="<<Eff<<", Purity="<<Purity);

	//Draw NN output histos
	TCanvas* NNOutPlot= new TCanvas("NNOutPlot","NNOutPlot");
	NNOutPlot->cd();

	TH2D* NNOutPlotBkg= new TH2D("NNOutPlotBkg","",100,-3,3,100,0,1);
	NNOutPlotBkg->SetXTitle("NNOut");
	NNOutPlotBkg->SetYTitle("entries");
	NNOutPlotBkg->Draw();

	NNOutputHisto_signal->SetFillStyle(3001);
	NNOutputHisto_signal->SetLineColor(kRed);
	NNOutputHisto_signal->DrawNormalized("hist same");

	NNOutputHisto_bkg->SetLineColor(kBlack);
	NNOutputHisto_bkg->DrawNormalized("hist same");
	
	return 0;

}//close TrainNN()
*/

void DrawData()
{
	//Define canvas
	std::vector<TCanvas*> histoPlots;
	TCanvas* histoPlot= 0;	
	TLegend* histoLegend= 0;
	
	//Draw var histos
	for(size_t j=0;j<dataVarHistos_signal.size();j++){
		TString canvasName= Form("Var%d",(int)(j+1));
		histoPlot= new TCanvas(canvasName,canvasName);

		histoPlot->cd();
		dataVarHistos_signal[j]->GetXaxis()->SetTitle(varNames[j]);
		dataVarHistos_signal[j]->SetLineColor(kRed);
		dataVarHistos_signal[j]->DrawNormalized("histo");
		dataVarHistos_bkg[j]->DrawNormalized("histo same");

		histoLegend= new TLegend(0.6,0.6,0.7,0.7);
		histoLegend->SetFillColor(0);
		histoLegend->SetTextSize(0.04);
		histoLegend->SetTextFont(52);
		histoLegend->AddEntry(dataVarHistos_signal[j],"real sources","L");
		histoLegend->AddEntry(dataVarHistos_bkg[j],"false sources","L");
		histoLegend->Draw("same");


		histoPlots.push_back(histoPlot);
	}//end loop vars


}//close DrawData()



