
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
#pragma link C++ class Caesar::ImgPeak+;
#pragma link C++ class Caesar::ImgBkgPars+;
//================

//== UTILS CLASSES ==
#pragma link C++ class Caesar::FileInfo+;
#pragma link C++ class Caesar::ProcMemInfo+;
#pragma link C++ class Caesar::SysUtils+;
#pragma link C++ class Caesar::AstroUtils+;
#pragma link C++ class Caesar::CodeUtils+;
#pragma link C++ class Caesar::MathUtils+;
#pragma link C++ class Caesar::GraphicsUtils+;
#pragma link C++ class Caesar::StatsUtils+;
#pragma link C++ class Caesar::EllipseUtils+;
#pragma link C++ class Caesar::ImgUtils+;

#pragma link C++ class Caesar::WCSUtils+;
#pragma link C++ class Caesar::WCS+;


#pragma link C++ class Caesar::ClippedStats<int>+;
#pragma link C++ class Caesar::ClippedStats<long int>+;
#pragma link C++ class Caesar::ClippedStats<float>+;
#pragma link C++ class Caesar::ClippedStats<double>+;

#pragma link C++ class Caesar::StatMoments<int>+;
#pragma link C++ class Caesar::StatMoments<long int>+;
#pragma link C++ class Caesar::StatMoments<float>+;
#pragma link C++ class Caesar::StatMoments<double>+;

#pragma link C++ class Caesar::BoxStats<int>+;
#pragma link C++ class Caesar::BoxStats<long int>+;
#pragma link C++ class Caesar::BoxStats<float>+;
#pragma link C++ class Caesar::BoxStats<double>+;


#pragma link C++ function Caesar::StatsUtils::GetMedianFast<float>;
#pragma link C++ function Caesar::StatsUtils::GetMedianFast<double>;
#pragma link C++ function Caesar::StatsUtils::GetMedianFast<int>;
#pragma link C++ function Caesar::StatsUtils::GetMedianFast<long int>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<float>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<double>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<int>;
#pragma link C++ function Caesar::StatsUtils::GetMedian<long int>;

#pragma link C++ function Caesar::StatsUtils::GetRangeMedian<float>;
#pragma link C++ function Caesar::StatsUtils::GetRangeMedian<double>;
#pragma link C++ function Caesar::StatsUtils::GetRangeMedian<int>;
#pragma link C++ function Caesar::StatsUtils::GetRangeMedian<long int>;

#pragma link C++ function Caesar::StatsUtils::ComputeBoxStats<float>;
#pragma link C++ function Caesar::StatsUtils::ComputeBoxStats<double>;
#pragma link C++ function Caesar::StatsUtils::ComputeBoxStats<int>;
#pragma link C++ function Caesar::StatsUtils::ComputeBoxStats<long int>;

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
//#pragma link C++ MACRO REGISTER_OPTION defined_in Caesar;
#pragma link C++ global REGISTER_OPTION defined_in Caesar;
#pragma link C++ global GET_OPTION defined_in Caesar;
#pragma link C++ global SET_OPTION defined_in Caesar;
#pragma link C++ global PRINT_OPTIONS defined_in Caesar;
#pragma link C++ global GET_OPTION_VALUE defined_in Caesar;

#pragma link C++ class Caesar::SBuffer+;
#pragma link C++ class Caesar::TaskData+;
//==============



//== IMG IO CLASSES ==
#pragma link C++ class Caesar::FITSHeader+;
#pragma link C++ class Caesar::FITSFileInfo+;
#pragma link C++ class Caesar::FITSReader+;
#pragma link C++ class Caesar::FITSWriter+;

#ifdef CASACORE_ENABLED
#pragma link C++ class Caesar::CasaReader+;
#endif
//==============

//== FILTER CLASSES ==
#pragma link C++ class Caesar::WTFilter+;
#pragma link C++ class Caesar::KirschFilter+;
#pragma link C++ class Caesar::LoGFilter+;
#pragma link C++ class Caesar::GradientFilter+;
#pragma link C++ class Caesar::MorphFilter+;
#pragma link C++ class Caesar::SaliencyFilter+;
#pragma link C++ class Caesar::GausFilter+;
//====================

//== IMG PROC CLASSES ==
//bkg finder
#pragma link C++ class Caesar::BkgSampleData+;
#pragma link C++ class vector<Caesar::BkgSampleData>+;
#pragma link C++ class Caesar::ImgBkgData+;
#pragma link C++ class Caesar::BkgFinder+;


//data class
#pragma link C++ class Caesar::Pixel+;
#pragma link C++ class vector<Caesar::Pixel*>+;
#pragma link C++ class map<int,Caesar::Pixel*>+;

#pragma link C++ class Caesar::Region+;
#pragma link C++ class vector<Caesar::Region>+;
#pragma link C++ class vector<Caesar::Region*>+;
#pragma link C++ class Caesar::RegionCollection+;

#pragma link C++ class Caesar::Blob+;
#pragma link C++ class vector<Caesar::Blob>+;
#pragma link C++ class vector<Caesar::Blob*>+;
#pragma link C++ class Caesar::Source+;
#pragma link C++ class vector<Caesar::Source>+;
#pragma link C++ class vector<Caesar::Source*>+;
#pragma link C++ enum Caesar::Source::SourceType+;
#pragma link C++ enum Caesar::Source::SourceFlag+;
#pragma link C++ enum Caesar::Source::SimSourceType+;
#pragma link C++ class Caesar::Contour+;
#pragma link C++ class vector<Caesar::Contour>+;
#pragma link C++ class vector<Caesar::Contour*>+;
#pragma link C++ struct ContourPointMatcher+;

#pragma link C++ class Caesar::SourceCube+;
#pragma link C++ class vector<Caesar::SourceCube>+;
#pragma link C++ class vector<Caesar::SourceCube*>+;

// - Source match classes
#pragma link C++ class Caesar::SourceMatchData+;
#pragma link C++ class vector<Caesar::SourceMatchData>+;
#pragma link C++ class vector<Caesar::SourceMatchData*>+;
#pragma link C++ class Caesar::SourceGroup+;
#pragma link C++ class vector<Caesar::SourceGroup>+;
#pragma link C++ class vector<Caesar::SourceGroup*>+;
#pragma link C++ class Caesar::ComponentMatchIndex+;
#pragma link C++ class vector<Caesar::ComponentMatchIndex>+;
#pragma link C++ class vector<Caesar::ComponentMatchIndex*>+;
#pragma link C++ class Caesar::ComponentMatchIndexGroup+;
#pragma link C++ class vector<Caesar::ComponentMatchIndexGroup>+;
#pragma link C++ class vector<Caesar::ComponentMatchIndexGroup*>+;
#pragma link C++ class Caesar::SpectralIndexData+;
#pragma link C++ class vector<Caesar::SpectralIndexData>+;
#pragma link C++ class vector<Caesar::SpectralIndexData*>+;
#pragma link C++ class Caesar::PLSpectralIndexData+;
#pragma link C++ class vector<Caesar::PLSpectralIndexData>+;
#pragma link C++ class vector<Caesar::PLSpectralIndexData*>+;
#pragma link C++ class Caesar::PolSpectralIndexData+;
#pragma link C++ class vector<Caesar::PolSpectralIndexData>+;
#pragma link C++ class vector<Caesar::PolSpectralIndexData*>+;
#pragma link C++ class Caesar::SASpectralIndexData+;
#pragma link C++ class vector<Caesar::SASpectralIndexData>+;
#pragma link C++ class vector<Caesar::SASpectralIndexData*>+;

#pragma link C++ class Caesar::SourceMatchPars+;
#pragma link C++ class vector<Caesar::SourceMatchPars>+;
#pragma link C++ class vector<Caesar::SourceMatchPars*>+;
#pragma link C++ class Caesar::SourceComponentMatchPars+;
#pragma link C++ class vector<Caesar::SourceComponentMatchPars>+;
#pragma link C++ class vector<Caesar::SourceComponentMatchPars*>+;
#pragma link C++ class Caesar::SourceMatchParsGroup+;
#pragma link C++ class vector<Caesar::SourceMatchParsGroup>+;
#pragma link C++ class vector<Caesar::SourceMatchParsGroup*>+;
#pragma link C++ class Caesar::SourceComponentMatchParsGroup+;
#pragma link C++ class vector<Caesar::SourceComponentMatchParsGroup>+;
#pragma link C++ class vector<Caesar::SourceComponentMatchParsGroup*>+;

// - Source export classes
#pragma link C++ class Caesar::SourceExporter+;
#pragma link C++ class Caesar::SourceImporter+;

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
#pragma link C++ function Caesar::BlobFinder::FindBlobsST<Blob>(Caesar::Img*,std::vector<Blob*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool,Caesar::Image*);
#pragma link C++ function Caesar::BlobFinder::FindBlobsST<Source>(Caesar::Img*,std::vector<Source*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool,Caesar::Image*);
#pragma link C++ function Caesar::BlobFinder::FindBlobs<Blob>(Caesar::Img*,std::vector<Blob*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool,Caesar::Image*);
#pragma link C++ function Caesar::BlobFinder::FindBlobs<Source>(Caesar::Img*,std::vector<Source*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool,Caesar::Image*);

#ifdef OPENMP_ENABLED
#pragma link C++ function Caesar::BlobFinder::FindBlobsMT<Blob>(Caesar::Img*,std::vector<Blob*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool,Caesar::Image*);
#pragma link C++ function Caesar::BlobFinder::FindBlobsMT<Source>(Caesar::Img*,std::vector<Source*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool,Caesar::Image*);
#endif

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
#pragma link C++ struct Caesar::SourceFitter::SourceFitData+;
#pragma link C++ enum Caesar::FitStatusFlag;
#pragma link C++ struct SourceFitOptions+;
//======================

//DS9Region
#pragma link C++ class Caesar::DS9RegionParser+;
#pragma link C++ class Caesar::DS9RegionMetaData+;
#pragma link C++ class Caesar::DS9Region+;
#pragma link C++ class vector<Caesar::DS9Region>+;
#pragma link C++ class vector<Caesar::DS9Region*>+;
#pragma link C++ class Caesar::DS9PolygonRegion+;
#pragma link C++ class vector<Caesar::DS9PolygonRegion>+;
#pragma link C++ class vector<Caesar::DS9PolygonRegion*>+;
#pragma link C++ class Caesar::DS9BoxRegion+;
#pragma link C++ class vector<DS9BoxRegion>+;
#pragma link C++ class vector<DS9BoxRegion*>+;
#pragma link C++ class Caesar::DS9CircleRegion+;
#pragma link C++ class vector<Caesar::DS9CircleRegion>+;
#pragma link C++ class vector<Caesar::DS9CircleRegion*>+;
#pragma link C++ class Caesar::DS9EllipseRegion+;
#pragma link C++ class vector<Caesar::DS9EllipseRegion>+;
#pragma link C++ class vector<Caesar::DS9EllipseRegion*>+;

//AstroObject
#pragma link C++ class Caesar::AstroObject+;
#pragma link C++ class Caesar::AstroObject*+;
#pragma link C++ class vector<Caesar::AstroObject>+;
#pragma link C++ class vector<Caesar::AstroObject*>+;
#pragma link C++ class Caesar::AstroObjectParser+;
#pragma link C++ class vector<Caesar::AstroObjectParser>+;
#pragma link C++ class vector<Caesar::AstroObjectParser*>+;
#pragma link C++ class vector<vector<Caesar::AstroObjectParser>>+;
#pragma link C++ class vector<vector<Caesar::AstroObjectParser*>>+;

//== CUT PARSER ==
#pragma link C++ class Caesar::Cut+;
#pragma link C++ class Caesar::Cut*+;
#pragma link C++ class std::map<std::string,Caesar::Cut*>+;

#pragma link C++ class Caesar::EqualityCut<int>+;
#pragma link C++ class Caesar::EqualityCut<long int>+;
#pragma link C++ class Caesar::EqualityCut<float>+;
#pragma link C++ class Caesar::EqualityCut<double>+;
#pragma link C++ class Caesar::EqualityCut<bool>+;
#pragma link C++ class Caesar::EqualityCut<int>*+;
#pragma link C++ class Caesar::EqualityCut<long int>*+;
#pragma link C++ class Caesar::EqualityCut<float>*+;
#pragma link C++ class Caesar::EqualityCut<double>*+;
#pragma link C++ class Caesar::EqualityCut<bool>*+;

#pragma link C++ class Caesar::BoundCut<int>+;
#pragma link C++ class Caesar::BoundCut<long int>+;
#pragma link C++ class Caesar::BoundCut<float>+;
#pragma link C++ class Caesar::BoundCut<double>+;
#pragma link C++ class Caesar::BoundCut<bool>+;
#pragma link C++ class Caesar::BoundCut<int>*+;
#pragma link C++ class Caesar::BoundCut<long int>*+;
#pragma link C++ class Caesar::BoundCut<float>*+;
#pragma link C++ class Caesar::BoundCut<double>*+;
#pragma link C++ class Caesar::BoundCut<bool>*+;
#pragma link C++ class Caesar::SingleBoundCut<int>+;
#pragma link C++ class Caesar::SingleBoundCut<long int>+;
#pragma link C++ class Caesar::SingleBoundCut<float>+;
#pragma link C++ class Caesar::SingleBoundCut<double>+;
#pragma link C++ class Caesar::SingleBoundCut<bool>+;
#pragma link C++ class Caesar::SingleBoundCut<int>*+;
#pragma link C++ class Caesar::SingleBoundCut<long int>*+;
#pragma link C++ class Caesar::SingleBoundCut<float>*+;
#pragma link C++ class Caesar::SingleBoundCut<double>*+;
#pragma link C++ class Caesar::SingleBoundCut<bool>*+;

#pragma link C++ class Caesar::CutFactory+;
#pragma link C++ class Caesar::CutParser+;
#pragma link C++ class Caesar::SourceSelector+;



#endif
