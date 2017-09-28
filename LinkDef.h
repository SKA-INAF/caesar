
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedef;
#pragma link C++ namespace Caesar;
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;


//== IMG CLASSES ==
#pragma link C++ class Caesar::MetaData+;
#pragma link C++ class Caesar::ImgMetaData+;
#pragma link C++ class Caesar::ImgStats+;
#pragma link C++ class Caesar::Image+;
#pragma link C++ class vector<Caesar::Image>+;
#pragma link C++ class vector<Caesar::Image*>+;
#pragma link C++ class Caesar::ImgRange+;
//================

//== UTILS CLASSES ==
#pragma link C++ class Caesar::FileInfo+;
#pragma link C++ class Caesar::SysUtils+;
#pragma link C++ class Caesar::AstroUtils+;
#pragma link C++ class Caesar::CodeUtils+;
#pragma link C++ class Caesar::MathUtils+;
#pragma link C++ class Caesar::GraphicsUtils+;
#pragma link C++ class Caesar::StatsUtils+;

#pragma link C++ class Caesar::ClippedStats<int>+;
#pragma link C++ class Caesar::ClippedStats<long int>+;
#pragma link C++ class Caesar::ClippedStats<float>+;
#pragma link C++ class Caesar::ClippedStats<double>+;

#pragma link C++ class Caesar::StatMoments<int>+;
#pragma link C++ class Caesar::StatMoments<long int>+;
#pragma link C++ class Caesar::StatMoments<float>+;
#pragma link C++ class Caesar::StatMoments<double>+;

#pragma link C++ function Caesar::StatsUtils::GetMedianFast<float>;
#pragma link C++ function Caesar::StatsUtils::GetMedianFast<double>;
#pragma link C++ function Caesar::StatsUtils::GetMedianFast<int>;
#pragma link C++ function Caesar::StatsUtils::GetMedianFast<long int>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<float>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<double>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<int>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<long int>;
#pragma link C++ function Caesar::StatsUtils::GetMAD<float>;
#pragma link C++ function Caesar::StatsUtils::GetMAD<double>;
#pragma link C++ function Caesar::StatsUtils::GetMAD<int>;
#pragma link C++ function Caesar::StatsUtils::GetMAD<long int>;
#pragma link C++ function Caesar::StatsUtils::GetMADFast<float>;
#pragma link C++ function Caesar::StatsUtils::GetMADFast<double>;
#pragma link C++ function Caesar::StatsUtils::GetMADFast<int>;
#pragma link C++ function Caesar::StatsUtils::GetMADFast<long int>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<float>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<double>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<int>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<long int>;
#pragma link C++ function Caesar::StatsUtils::GetBiWeightEstimators<float>;
#pragma link C++ function Caesar::StatsUtils::GetBiWeightEstimators<double>;
#pragma link C++ function Caesar::StatsUtils::GetBiWeightEstimators<int>;
#pragma link C++ function Caesar::StatsUtils::GetBiWeightEstimators<long int>;

#pragma link C++ class Caesar::ClippedStats+;
#pragma link C++ class Caesar::ClippedStats*+;
#pragma link C++ class Caesar::ClippedStats<int>+;
#pragma link C++ class Caesar::ClippedStats<long int>+;
#pragma link C++ class Caesar::ClippedStats<float>+;
#pragma link C++ class Caesar::ClippedStats<double>+;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<int>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<long int>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<float>;
#pragma link C++ function Caesar::StatsUtils::GetClippedEstimators<double>;

#pragma link C++ class Caesar::Graph+;
#pragma link C++ class Caesar::Graph*;
//==============

//== LOGGER ==
#pragma link C++ class Caesar::Logger+;
#pragma link C++ class Caesar::SysLogger+;
#pragma link C++ class Caesar::FileLogger+;
#pragma link C++ class Caesar::ConsoleLogger+;
#pragma link C++ class Caesar::LoggerManager+;

//== CONFIG PARSER ==
#pragma link C++ class Caesar::OptionBase+;
#pragma link C++ class Caesar::OptionBase*+;
#pragma link C++ class std::map<std::string,Caesar::OptionBase*>+;

#pragma link C++ class Caesar::Option<int>+;
#pragma link C++ class Caesar::Option<long int>+;
#pragma link C++ class Caesar::Option<float>+;
#pragma link C++ class Caesar::Option<double>+;
#pragma link C++ class Caesar::Option<std::string>+;
#pragma link C++ class Caesar::Option<bool>+;
#pragma link C++ class Caesar::Option<char*>+;

#pragma link C++ class Caesar::Option<std::vector<int>>+;
#pragma link C++ class Caesar::Option<std::vector<long int>>+;
#pragma link C++ class Caesar::Option<std::vector<float>>+;
#pragma link C++ class Caesar::Option<std::vector<double>>+;
#pragma link C++ class Caesar::Option<std::vector<std::string>>+;
#pragma link C++ class Caesar::Option<std::vector<bool>>+;
#pragma link C++ class Caesar::Option<std::vector<char*>>+;

#pragma link C++ class Caesar::Option<int>*+;
#pragma link C++ class Caesar::Option<long int>*+;
#pragma link C++ class Caesar::Option<float>*+;
#pragma link C++ class Caesar::Option<double>*+;
#pragma link C++ class Caesar::Option<std::string>*+;
#pragma link C++ class Caesar::Option<bool>*+;
#pragma link C++ class Caesar::Option<char*>*+;

#pragma link C++ class Caesar::Option<std::vector<int>>*+;
#pragma link C++ class Caesar::Option<std::vector<long int>>*+;
#pragma link C++ class Caesar::Option<std::vector<float>>*+;
#pragma link C++ class Caesar::Option<std::vector<double>>*+;
#pragma link C++ class Caesar::Option<std::vector<std::string>>*+;
#pragma link C++ class Caesar::Option<std::vector<bool>>*+;
#pragma link C++ class Caesar::Option<std::vector<char*>>*+;


#pragma link C++ class Caesar::OptionHelper<int>+;
#pragma link C++ class Caesar::OptionHelper<long int>+;
#pragma link C++ class Caesar::OptionHelper<float>+;
#pragma link C++ class Caesar::OptionHelper<double>+;
#pragma link C++ class Caesar::OptionHelper<std::string>+;
#pragma link C++ class Caesar::OptionHelper<bool>+;
#pragma link C++ class Caesar::OptionHelper<char*>+;

#pragma link C++ class Caesar::OptionHelper<std::vector<int>>+;
#pragma link C++ class Caesar::OptionHelper<std::vector<long int>>+;
#pragma link C++ class Caesar::OptionHelper<std::vector<float>>+;
#pragma link C++ class Caesar::OptionHelper<std::vector<double>>+;
#pragma link C++ class Caesar::OptionHelper<std::vector<std::string>>+;
#pragma link C++ class Caesar::OptionHelper<std::vector<bool>>+;
#pragma link C++ class Caesar::OptionHelper<std::vector<char*>>+;

#pragma link C++ class Caesar::OptionHelper<int>*+;
#pragma link C++ class Caesar::OptionHelper<long int>*+;
#pragma link C++ class Caesar::OptionHelper<float>*+;
#pragma link C++ class Caesar::OptionHelper<double>*+;
#pragma link C++ class Caesar::OptionHelper<std::string>*+;
#pragma link C++ class Caesar::OptionHelper<bool>*+;
#pragma link C++ class Caesar::OptionHelper<char*>*+;

#pragma link C++ class Caesar::OptionHelper<std::vector<int>>*+;
#pragma link C++ class Caesar::OptionHelper<std::vector<long int>>*+;
#pragma link C++ class Caesar::OptionHelper<std::vector<float>>*+;
#pragma link C++ class Caesar::OptionHelper<std::vector<double>>*+;
#pragma link C++ class Caesar::OptionHelper<std::vector<std::string>>*+;
#pragma link C++ class Caesar::OptionHelper<std::vector<bool>>*+;
#pragma link C++ class Caesar::OptionHelper<std::vector<char*>*>+;

#pragma link C++ class Caesar::OptionFactory+;
#pragma link C++ class Caesar::ConfigParser+;
#pragma link C++ MACRO REGISTER_OPTION defined_in Caesar;

#pragma link C++ class Caesar::SBuffer+;
#pragma link C++ class Caesar::TaskData+;
//==============



//== IMG IO CLASSES ==
#pragma link C++ class Caesar::FITSHeader+;
#pragma link C++ class Caesar::FITSFileInfo+;
#pragma link C++ class Caesar::FITSReader+;
#pragma link C++ class Caesar::FITSWriter+;
//==============

//== FILTER CLASSES ==
#pragma link C++ class Caesar::WTFilter+;
#pragma link C++ class Caesar::KirschFilter+;
#pragma link C++ class Caesar::LoGFilter+;
#pragma link C++ class Caesar::GradientFilter+;
#pragma link C++ class Caesar::MorphFilter+;
#pragma link C++ class Caesar::SaliencyFilter+;
//====================

//== IMG PROC CLASSES ==
//bkg finder
#pragma link C++ class Caesar::BkgSampleData+;
#pragma link C++ class vector<Caesar::BkgSampleData>+;
//#pragma link C++ class Caesar::BkgData+;
#pragma link C++ class Caesar::ImgBkgData+;
#pragma link C++ class Caesar::BkgFinder+;
//#pragma link C++ function Caesar::BlobFinder::FindBlob<Blob>(Caesar::Img*,std::vector<Blob*>&,Caesar::Img*,Caesar::BkgData*,double,double,int,bool,bool);
//#pragma link C++ function Caesar::BlobFinder::FindBlob<Source>(Caesar::Img*,std::vector<Source*>&,Caesar::Img*,Caesar::BkgData*,double,double,int,bool,bool);
#pragma link C++ function Caesar::BlobFinder::FindBlob<Blob>(Caesar::Img*,std::vector<Blob*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool);
#pragma link C++ function Caesar::BlobFinder::FindBlob<Source>(Caesar::Img*,std::vector<Source*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool);

//data class
#pragma link C++ class Caesar::Pixel+;
#pragma link C++ class vector<Caesar::Pixel*>+;
#pragma link C++ class map<int,Caesar::Pixel*>+;

#pragma link C++ class Caesar::Region+;
#pragma link C++ class Caesar::vector<Region>+;
#pragma link C++ class Caesar::vector<Region*>+;
#pragma link C++ class Caesar::RegionCollection+;

#pragma link C++ class Caesar::Blob+;
#pragma link C++ class Caesar::vector<Blob>+;
#pragma link C++ class Caesar::vector<Blob*>+;
#pragma link C++ class Caesar::Source+;
#pragma link C++ class Caesar::vector<Source>+;
#pragma link C++ class Caesar::vector<Source*>+;
#pragma link C++ enum Caesar::Source::SourceType+;
#pragma link C++ enum Caesar::Source::SourceFlag+;
#pragma link C++ enum Caesar::Source::SimSourceType+;
#pragma link C++ class Caesar::Contour+;

//slic generator
#pragma link C++ class Caesar::SLICData+;
#pragma link C++ class Caesar::SLICContourData+;
#pragma link C++ class Caesar::SLICNeighborData+;
#pragma link C++ class Caesar::SLICNeighborCollection+;
#pragma link C++ class Caesar::SLICSimilarityData+;
#pragma link C++ enum Caesar::SLICEdgeModel+;
#pragma link C++ class Caesar::SLIC+;
#pragma link C++ class Caesar::SLICSegmenter+;

//blob/source finder
#pragma link C++ class Caesar::BlobFinder+;
#pragma link C++ class Caesar::SFinder+;
//======================

//ChanVese segmenter
#pragma link C++ class Caesar::ChanVeseSegmenter+;
#pragma link C++ class Caesar::LRACSegmenter+;
//======================

//Source fitter
#pragma link C++ class Caesar::SourceComponentPars+;
#pragma link C++ class Caesar::SourceFitPars+;
#pragma link C++ class Caesar::SourceFitter+;
#pragma link C++ enum Caesar::FitStatusFlag;
//======================


#endif
