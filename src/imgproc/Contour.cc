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
* @file Contour.cc
* @class Contour
* @brief Contour
*
* Class representing image contour with methods for morphological parameter extraction 
* @author S. Riggi
* @date 11/07/2015
*/

#include <Contour.h>
#include <MathUtils.h>
#include <Logger.h>

#include <TMath.h>
#include <TEllipse.h>
#include <TPolyLine.h>
#include <TPaveText.h>
#include <TVectorD.h>
#include <TVector2.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;


ClassImp(Caesar::Contour)

namespace Caesar {

Contour::Contour() : TObject() {

	Init();

}//close costructor


Contour::~Contour(){

}//close destructor


Contour::Contour(const Contour& contour) : TObject(contour) {
  // Contour copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((Contour&)contour).Copy(*this);
}


void Contour::Copy(TObject &obj) const {
   
	// Copy this contour to contour
  TObject::Copy((Contour&)obj);
  ((Contour&)obj).HasParameters = HasParameters;
	((Contour&)obj).Area = Area;
  ((Contour&)obj).Perymeter = Perymeter;
  ((Contour&)obj).IsConvexContour = IsConvexContour;
	((Contour&)obj).CircularityRatio = CircularityRatio;
	((Contour&)obj).BoundingBoxCenter = BoundingBoxCenter;
	((Contour&)obj).BoundingBoxMaj = BoundingBoxMaj;
	((Contour&)obj).BoundingBoxMin = BoundingBoxMin;
	((Contour&)obj).BoundingBoxAngle = BoundingBoxAngle;
	((Contour&)obj).Elongation = Elongation;
	((Contour&)obj).Rectangularity = Rectangularity;
	((Contour&)obj).Roundness = Roundness;
	((Contour&)obj).Eccentricity = Eccentricity;
	((Contour&)obj).TiltAngle = TiltAngle;
	((Contour&)obj).HasEllipseFit = HasEllipseFit;
	((Contour&)obj).EllipseCenter = EllipseCenter;
	((Contour&)obj).EllipseMajAxis = EllipseMajAxis;
	((Contour&)obj).EllipseMinAxis = EllipseMinAxis;
	((Contour&)obj).EllipseRotAngle = EllipseRotAngle;
	((Contour&)obj).EllipseFitRedChi2 = EllipseFitRedChi2;
	((Contour&)obj).EllipseAreaRatio = EllipseAreaRatio;

	((Contour&)obj).m00 = m00;	
	((Contour&)obj).m10 = m10;
	((Contour&)obj).m01 = m01;
	((Contour&)obj).m20 = m20;
	((Contour&)obj).m11 = m11;
	((Contour&)obj).m02 = m02;
	((Contour&)obj).m30 = m30;
	((Contour&)obj).m21 = m21;
	((Contour&)obj).m12 = m12;
	((Contour&)obj).m03 = m03;
	
	((Contour&)obj).mu20 = mu20;
	((Contour&)obj).mu11 = mu11;
	((Contour&)obj).mu02 = mu02;
	((Contour&)obj).mu30 = mu30;
	((Contour&)obj).mu21 = mu21;
	((Contour&)obj).mu12 = mu12;
	((Contour&)obj).mu03 = mu03;
	
	((Contour&)obj).nu20 = nu20;
	((Contour&)obj).nu11 = nu11;
	((Contour&)obj).nu02 = nu02;
	((Contour&)obj).nu30 = nu30;
	((Contour&)obj).nu21 = nu21;
	((Contour&)obj).nu12 = nu12;
	((Contour&)obj).nu03 = nu03;

	((Contour&)obj).HuMoments = HuMoments;
	((Contour&)obj).BoundingBoxVertex = BoundingBoxVertex;
	((Contour&)obj).Centroid = Centroid;
	((Contour&)obj).RealFDs = RealFDs;
	((Contour&)obj).ImagFDs = ImagFDs;
	((Contour&)obj).ModFDs = ModFDs;
	((Contour&)obj).BendingEnergies = BendingEnergies;
	((Contour&)obj).CentroidDistanceModFDs = CentroidDistanceModFDs;

	((Contour&)obj).HasBEPars= HasBEPars;	
	((Contour&)obj).HasCentroidDistanceFDPars= HasCentroidDistanceFDPars;		
	((Contour&)obj).HasFDPars= HasFDPars;	

	((Contour&)obj).m_Points = m_Points;

}//close Copy()

Contour& Contour::operator=(const Contour& contour) { 
	// Operator =
  if (this != &contour)  ((Contour&)contour).Copy(*this);
  return *this;
}


void Contour::Init(){

	m_Points.clear();	
	HasParameters= false;	
	HasFDPars= false;
	HasEllipseFit= false;
	EllipseMajAxis= 0;
	EllipseMinAxis= 0;
	EllipseRotAngle= 0;
	EllipseFitRedChi2= 1.e+99;
	EllipseAreaRatio= 0;
	Area= 0.;
	Perymeter= 0.;
	IsConvexContour= false;
	CircularityRatio= -999;
	BoundingBoxMaj= -999;
	BoundingBoxMin= -999;
	BoundingBoxAngle= -999;
	Elongation= -999;
	Rectangularity= -999;
	Roundness= -999;
	Eccentricity= -999;
	TiltAngle= -999;
	//Moments= ...
	for(int k=0;k<7;k++) HuMoments.push_back(0.);//HuMoments[k]= 0;
	for(int k=0;k<4;k++) BoundingBoxVertex.push_back(TVector2(0,0));//BoundingBoxVertex[k]= TVector2(0,0);	
	BoundingBoxCenter= TVector2(0,0);
	Centroid= TVector2(0,0);
	//FDs.clear();
	RealFDs.clear();
	ImagFDs.clear();
	ModFDs.clear();
		
	HasFDPars= false;
	HasBEPars= false;
	HasCentroidDistanceFDPars= false;

}//close Init()

TGraph* Contour::GetGraph(){

	//Check number of contour pts
	int nContourPts= m_Points.size();
	if(nContourPts<=0) {
		WARN_LOG("No contour points available (did you fill the contour?), returning null ptr graph!");
		return 0;
	}
	//Fill contour graph
	TGraph* ContourGraph= new TGraph(nContourPts+1);
	for(int i=0;i<nContourPts;i++){
		ContourGraph->SetPoint(i,m_Points[i].X(),m_Points[i].Y());
	}
	ContourGraph->SetPoint(nContourPts,m_Points[0].X(),m_Points[0].Y());//Add another point as OpenCV does not close the contour

	return ContourGraph;

}//close Contour::GetGraph()


TPolyLine* Contour::GetBoundingBoxLine(){

	if(!HasPoints() || !HasParameters){
		WARN_LOG("No contour points and/or parameters available (did you fill the contour and compute its parameters?)!");
		return 0;
	}

	TPolyLine* boundingRectangle= new TPolyLine(5);	
	for(int k=0;k<4;k++) boundingRectangle->SetPoint(k,BoundingBoxVertex[k].X(),BoundingBoxVertex[k].Y());
	boundingRectangle->SetPoint(4,BoundingBoxVertex[0].X(),BoundingBoxVertex[0].Y());	
	return boundingRectangle;

}//close Contour::GetBoundingBoxLine()

TPaveText* Contour::GetParamInfoBox(){

	TPaveText* infoBox = new TPaveText(0.15,0.7,0.75,0.85,"NDC");
	infoBox->AddText(Form("Eccentricity: %1.2f",Eccentricity));
	infoBox->AddText(Form("Elongation: %1.2f",Elongation));
	infoBox->AddText(Form("TiltAngle(deg): %1.2f",TiltAngle));
	infoBox->AddText(Form("Rectangularity: %1.2f",Rectangularity));
	infoBox->AddText(Form("Roundness: %1.2f, CircRatio: %1.2f",Roundness,CircularityRatio));
	infoBox->AddText(Form("BoundingBox: (%1.2f,%1.2f,%1.2f)",BoundingBoxMaj,BoundingBoxMin,BoundingBoxAngle));
	infoBox->AddText(Form("HuMoments: (%1.2g,%1.2g,%1.2g,%1.2g,%1.2g,%1.2g,%1.2g)",HuMoments[0],HuMoments[1],HuMoments[2],HuMoments[3],HuMoments[4],HuMoments[5],HuMoments[6]));
	if(HasEllipseFit) infoBox->AddText(Form("(x0,y0,a,b,Theta): (%1.2f,%1.2f,%1.2f,%1.2f,%1.2f), EllipseAreaRatio: %1.2f",EllipseCenter.X(),EllipseCenter.Y(),EllipseMajAxis/2,EllipseMinAxis/2,EllipseRotAngle,EllipseAreaRatio));
	infoBox->SetTextAlign(12);
	infoBox->SetTextSize(0.02);
	infoBox->SetTextFont(52);
	infoBox->SetFillColor(0);
	infoBox->SetBorderSize(1);

	return infoBox;

}//close GetParamInfoBox()

int Contour::ComputeParameters(){

	int nContourPts= m_Points.size();
	if(nContourPts<=0) {
		WARN_LOG("No contour points available (did you fill the contour?)!");
		return -1;
	}
	if(nContourPts<4) {
		WARN_LOG("Too few contour points available (n="<<nContourPts<<") to get any reliable parameter estimate!");
		return -1;
	}

	//## Copy points to cv::Point list
	DEBUG_LOG("Copying contour points in cv::Point list...");
	std::vector<cv::Point2f> points;
	for(unsigned int i=0;i<m_Points.size();i++){
		double x= m_Points[i].X();
		double y= m_Points[i].Y();
		points.push_back( cv::Point2f(x,y) );	
	}

	int status= 0;


	try{
		//## Compute Shape Parameters
		ComputeShapeParams(points);
	
		//## Compute fitted ellipse
		//if(ComputeFittedEllipse()<0){
		//	WARN_LOG("Failed to fit an ellipse to the contour!");
		//	status= -1;
		//}

		//## Compute moment params (moments, HuMoment, eccentricity)
		ComputeMomentParams(points);

		//## Compute Fourier descriptor
		//ComputeFourierDescriptors();
		//ComputeCentroidDistanceFD();

		//## Compute average bending energy
		//ComputeBendingEnergy();		

	}//close try block
	catch(cv::Exception ex){//something goes wrong!
  	ERROR_LOG("Computing contour parameters failed (err: "<<ex.msg <<")");
		HasParameters= false;
		return -1;
  }		
	catch(...){//something goes wrong!
  	ERROR_LOG("Unknown exception caught while computing contour parameters!");
		HasParameters= false;
		return -1;
  }

	HasParameters= true;

	return 0;

}//close ComputeParameters()

//void Contour::ComputeShapeParams(){
void Contour::ComputeShapeParams(std::vector<cv::Point2f>const & points){

	/*
	//Copy points to cv::Point list
	std::vector<cv::Point2f> points;
	for(unsigned int i=0;i<m_Points.size();i++){
		double x= m_Points[i].X();
		double y= m_Points[i].Y();
		points.push_back( cv::Point2f(x,y) );	
	}
	*/

	//Compute Area
	Area= cv::contourArea(points,false);

	//Compute perimeter
	bool isContourClosed= true;//assuming contours are always closed
	Perymeter= cv::arcLength(points, isContourClosed);

	//Compute circularity ratio
	CircularityRatio= 4*TMath::Pi()*Area/pow(Perymeter,2);

	//Compute Bounding Box
	cv::RotatedRect MinBoundingRect= cv::minAreaRect(points);//rotated bounding box

	cv::Point2f bb[4]; 
	MinBoundingRect.points(bb);//bounding box vertexes
	for(int k=0;k<4;k++) BoundingBoxVertex[k]= TVector2(bb[k].x,bb[k].y);

	double MinBoundingRect_height= MinBoundingRect.size.height;
	double MinBoundingRect_width= MinBoundingRect.size.width;
	BoundingBoxCenter= TVector2(MinBoundingRect.center.x,MinBoundingRect.center.y);
	BoundingBoxMaj= std::max(MinBoundingRect_height,MinBoundingRect_width);
	BoundingBoxMin= std::min(MinBoundingRect_height,MinBoundingRect_width);
	BoundingBoxAngle= MinBoundingRect.angle;//counterclockwise in degree 
	if(MinBoundingRect.size.width < MinBoundingRect.size.height){
    BoundingBoxAngle+= 180;
  }
	else{
		BoundingBoxAngle+= 90;  
	}
	
	//Compute elongation
	Elongation= 1.-BoundingBoxMin/BoundingBoxMaj;

	//Compute rectangularity
	Rectangularity= Area/(BoundingBoxMaj*BoundingBoxMin);

	//Compute roundness
	Roundness= 4.*Area/(TMath::Pi()*pow(BoundingBoxMaj,2));

}//close ComputeShapeParams()


void Contour::ComputeMomentParams(std::vector<cv::Point2f>const & points){

	//Compute moments & HuMoments
	ComputeMoments(points);

	//Compute eccentricity
	ComputeEccentricity();

}//close ComputeMomentParams()


//void Contour::ComputeMoments(){
void Contour::ComputeMoments(std::vector<cv::Point2f>const & points){
	
	//Copy points to cv::Point list
	/*
	std::vector<cv::Point2f> points;
	for(size_t i=0;i<m_Points.size();i++){
		double x= m_Points[i].X();
		double y= m_Points[i].Y();
		points.push_back( cv::Point2f(x,y) );	
	}
	*/

	//====================================
	//==   COMPUTE MOMENTS
	//====================================
	// - spatial moments: m00, m10, m01, m20, m11, m02, m30, m21, m12, m03
  // - central moments: mu20, mu11, mu02, mu30, mu21, mu12, mu03
  // - central normalized moments: nu20, nu11, nu02, nu30, nu21, nu12, nu03
	cv::Moments moments = cv::moments(points);
	
	
	//====================================
	//==   COMPUTE HU MOMENTS
	//====================================
	double humoments_array[7];
	cv::HuMoments(moments, humoments_array);
	for(int i=0;i<7;i++) HuMoments[i]= humoments_array[i]; 
	
	m00= moments.m00;
	m10= moments.m10;
	m01= moments.m01;
	m20= moments.m20;
	m11= moments.m11;
	m02= moments.m02;
	m30= moments.m30;
	m21= moments.m21;
	m12= moments.m12;
	m03= moments.m03;
	Centroid= TVector2(m10/m00,m01/m00);
		
	mu20= moments.mu20;
	mu11= moments.mu11;
	mu02= moments.mu02;
	mu30= moments.mu30;
	mu21= moments.mu21;
	mu12= moments.mu12;
  mu03= moments.mu03;
	
	nu20= moments.nu20;
	nu11= moments.nu11;
	nu02= moments.nu02;
	nu30= moments.nu30;
	nu21= moments.nu21;
	nu12= moments.nu12;
	nu03= moments.nu03;

}//close ComputeMoments()


void Contour::ComputeEccentricity(){

	//Compute covariance matrix and its eigenvectors
	double Cxx= mu20/m00;
	double Cyy= mu02/m00;
	double Cxy= mu11/m00;
	double delta= sqrt(4*Cxy*Cxy+(Cxx-Cyy)*(Cxx-Cyy));
	double lambda1= ((Cxx+Cyy) + delta)/2.; 
	double lambda2= ((Cxx+Cyy) - delta)/2.;	
	Eccentricity= sqrt(1-lambda2/lambda1);
	TiltAngle= 0.5*atan(2.*Cxy/(Cxx-Cyy))*TMath::RadToDeg();

}//close ComputeEccentricity()


TEllipse* Contour::GetFittedEllipse(){

	//Check if ellipse fit succeeded
	if(!HasEllipseFit) return nullptr;

	//Return the ellipse
	TEllipse* ellipse= new TEllipse;
	ellipse->SetX1(EllipseCenter.X());
	ellipse->SetY1(EllipseCenter.Y());
	ellipse->SetTheta(EllipseRotAngle);
	ellipse->SetR1(EllipseMajAxis/2.);
	ellipse->SetR2(EllipseMinAxis/2.);	

	return ellipse;

}//close GetFittedEllipse()

int Contour::ComputeFittedEllipse(){
	
	//## Fit an ellipse to contour points
	HasEllipseFit= true;//checked in internal routines and swithed to false if the fit fails
	
	// Get the contour graph
	TGraph* contourGraph= GetGraph();
	if(!contourGraph){
		ERROR_LOG("Cannot get contour graph!");
		return -1;
	}

	//Get ellipse conic
	TVectorD conic;
	if(EllipseFitter(conic,contourGraph)<0){
		ERROR_LOG("Failed to compute the ellipse conic!");
		if(contourGraph) contourGraph->Delete();
		return -1;
	}

	//Convert to parametric representation
  TVectorD ellipseParams;
	if(ConicToParametric(ellipseParams,conic)<0){
		WARN_LOG("Failed to convert ellipse from conic to parametric representation!");
		if(contourGraph) contourGraph->Delete();	
		return -1;
	}

	//Check for errors
	if(!HasEllipseFit){
		WARN_LOG("Ellipse fit failed!");
		if(contourGraph) contourGraph->Delete();	
		return -1;
	}

	//ellipse[0] = x0; // ellipse's "x" center
  //ellipse[1] = y0; // ellipse's "y" center
  //ellipse[2] = a; // ellipse's "semimajor" axis along "x"
  //ellipse[3] = b; // ellipse's "semiminor" axis along "y"
  //ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)
	//EllipseCenter= cv::Point2f(ellipseParams[0],ellipseParams[1]);

	EllipseCenter= TVector2(ellipseParams[0],ellipseParams[1]);

	//EllipseMajAxis= 2*std::max(ellipseParams[2],ellipseParams[3]);
	//EllipseMinAxis= 2*std::min(ellipseParams[2],ellipseParams[3]);
	EllipseMajAxis= 2*ellipseParams[2];
	EllipseMinAxis= 2*ellipseParams[3];
	EllipseRotAngle= ellipseParams[4];

	//Compute the area ratio
	double EllipseArea= TMath::Pi()*ellipseParams[2]*ellipseParams[3];
	EllipseAreaRatio= Area/EllipseArea;
	
	//Compute the fit residual
	double chi2= EllipseFitChi2(contourGraph,ellipseParams);
	double ndf= contourGraph->GetN()-ellipseParams.GetNoElements();
	EllipseFitRedChi2= chi2/ndf;

	std::stringstream sstream;
	sstream<<"ELLIPSE FIT: ";
	sstream<<"(x0,y0)=("<<ellipseParams[0]<<","<<ellipseParams[1]<<"), ";
	sstream<<"a="<<ellipseParams[2]<<" b="<<ellipseParams[3]<<", ";
	sstream<<"theta(deg)="<<ellipseParams[4]<<", ";
	sstream<<"chi2="<<chi2<<" redchi2="<<EllipseFitRedChi2;
	DEBUG_LOG(sstream.str());
	
	/*
	//OpenCV method
	std::vector<cv::Point2f> hull;
	cv::convexHull(m_Points, hull);
	//cv::RotatedRect fittedEllipse = cv::fitEllipse(cv::Mat(m_Points));
	cv::RotatedRect fittedEllipse = cv::fitEllipse(hull);
	EllipseCenter= fittedEllipse.center;
	EllipseMajAxis= std::max(fittedEllipse.size.width,fittedEllipse.size.height);
	EllipseMinAxis= std::min(fittedEllipse.size.width,fittedEllipse.size.height);
	EllipseRotAngle= fittedEllipse.angle;//in degrees
	EllipseArea= TMath::Pi()*fittedEllipse.size.width/2*fittedEllipse.size.height/2;
	EllipseAreaRatio= Area/EllipseArea;
	cout<<"*** OPENCV ELLIPSE FIT ***"<<endl;
	cout<<"(x0,y0)=("<<EllipseCenter.x<<","<<EllipseCenter.y<<")"<<endl;
	cout<<"a="<<EllipseMajAxis/2<<" b="<<EllipseMinAxis/2<<endl;
	cout<<"theta(deg)="<<EllipseRotAngle<<endl;
	cout<<"*******************"<<endl;	
	*/
	if(contourGraph) contourGraph->Delete();
	
	return 0;

}//close ComputeFittedEllipse()

 
int Contour::EllipseFitter(TVectorD& ellipse,TGraph* contourGraph){

  //TVectorD ellipse;
  if (!contourGraph) {
		ERROR_LOG("Null ptr to input contour graph, no ellipse fit will be performed!");
		HasEllipseFit= false;
		return -1;
	}

	//Check number of points
	int N = contourGraph->GetN();
  if(N<6) {
		WARN_LOG("Contour has less then 6 points, cannot fit ellipse!");
		HasEllipseFit= false;
		return -1;
	}
  
	int i= 0;
  double tmp= 0;
  double xmin, xmax, ymin, ymax, X0, Y0;
  contourGraph->ComputeRange(xmin, ymin, xmax, ymax);

	#if 1 // 0 or 1 
	  X0 = (xmax + xmin) / 2.0;
  	Y0 = (ymax + ymin) / 2.0;
	#else // 0 or 1 
  	X0 = Y0 = 0.0;
	#endif // 0 or 1
  
  TMatrixD D1(N, 3); // quadratic part of the design matrix
  TMatrixD D2(N, 3); // linear part of the design matrix
  
  for(i=0;i<N;i++) {
    double x = (contourGraph->GetX())[i] - X0;
    double y = (contourGraph->GetY())[i] - Y0;
    D1[i][0] = x*x;
    D1[i][1] = x*y;
    D1[i][2] = y*y;
    D2[i][0] = x;
    D2[i][1] = y;
    D2[i][2] = 1.0;
  }
  
  // Quadratic part of the scatter matrix
  TMatrixD S1(TMatrixD::kAtA, D1);

  // Combined part of the scatter matrix
  TMatrixD S2(D1, TMatrixD::kTransposeMult, D2);

  // Linear part of the scatter matrix
  TMatrixD S3(TMatrixD::kAtA, D2);
  S3.Invert(&tmp); 
	S3*= -1.0;
  if (tmp == 0.0) {
    WARN_LOG("Linear part of the scatter matrix is singular!");
		HasEllipseFit= false;
    return -1;
  }
  // For getting a2 from a1
  TMatrixD T(S3, TMatrixD::kMultTranspose, S2);

  // Reduced scatter matrix
  TMatrixD M(S2, TMatrixD::kMult, T); M += S1;

  // Premultiply by inv(C1)
  for (i=0;i<3;i++) {
    tmp = M[0][i]/2.0;
    M[0][i] = M[2][i]/2.0;
    M[2][i] = tmp;
    M[1][i]*= -1.0;
  }

  // Solve eigensystem
  TMatrixDEigen eig(M); // note: eigenvectors are not normalized
  const TMatrixD &evec = eig.GetEigenVectors();
  if ((eig.GetEigenValuesIm()).Norm2Sqr() != 0.0) {
    WARN_LOG("Eigenvalues have nonzero imaginary parts!");
		HasEllipseFit= false;
    return -1;
  }

  // Evaluate a’Ca (in order to find the eigenvector for min. pos. eigenvalue)
  for (i=0; i<3; i++) {
    tmp = 4.0 * evec[0][i] * evec[2][i] - evec[1][i] * evec[1][i];
    if (tmp > 0.0) break;
  }
  if (i > 2) {
    WARN_LOG("No min. pos. eigenvalue found!");
    // i = 2;
		HasEllipseFit= false;
    return -1;
  }

  // Eigenvector for min. pos. eigenvalue
  TVectorD a1(TMatrixDColumn_const(evec, i));
  tmp = a1.Norm2Sqr();
  if (tmp > 0.0) {
    a1*= 1.0/sqrt(tmp); // normalize this eigenvector
  } 
	else {
    WARN_LOG("Eigenvector for min. pos. eigenvalue is NULL!");
		HasEllipseFit= false;
    return -1;
  }
  TVectorD a2(T*a1);
  
  // Ellipse coefficients
  ellipse.ResizeTo(8);
  ellipse[0] = X0; // "X0"
  ellipse[1] = Y0; // "Y0"
  ellipse[2] = a1[0]; // "A"
  ellipse[3] = a1[1]; // "B"
  ellipse[4] = a1[2]; // "C"
  ellipse[5] = a2[0]; // "D"
  ellipse[6] = a2[1]; // "E"
  ellipse[7] = a2[2]; // "F"
  
  //return ellipse;
	return 0;

}//close EllipseFitter()


int Contour::ConicToParametric(TVectorD& ellipse,const TVectorD &conic) {
  
	//TVectorD ellipse;  
  if (conic.GetNrows() != 8) {
    ERROR_LOG("Improper input vector length!");
		HasEllipseFit= false;
    //return ellipse;
		return -1;
  }
  
  double a, b, theta;
  double x0 = conic[0]; // = X0
  double y0 = conic[1]; // = Y0
  
  // http://mathworld.wolfram.com/Ellipse.html
  double A = conic[2];
  double B = conic[3] / 2.0;
  double C = conic[4];
  double D = conic[5] / 2.0;
  double F = conic[6] / 2.0;
  double G = conic[7];
  
  double J = B * B - A * C;
  double Delta = A * F * F + C * D * D + J * G - 2.0 * B * D * F;
  double I = - (A + C);
  
  // http://mathworld.wolfram.com/QuadraticCurve.html
  if (!( (Delta != 0.0) && (J < 0.0) && (I != 0.0) && (Delta / I < 0.0) )) {
    ERROR_LOG("Ellipse (real) specific constraints not met!");
		HasEllipseFit= false;
    //return ellipse;
		return -1;
  }
  
  x0 += (C * D - B * F) / J;
  y0 += (A * F - B * D) / J;
  
  double tmp = sqrt((A - C) * (A - C) + 4.0 * B * B);
  a = sqrt(2.0 * Delta / J / (I + tmp));
  b = sqrt(2.0 * Delta / J / (I - tmp));
  
  theta = 0.0;
  if (B != 0.0) {
    tmp = (A - C) / 2.0 / B;
    theta = -45.0 * (std::atan(tmp) / TMath::PiOver2());
    if (tmp < 0.0) { theta -= 45.0; } else { theta += 45.0; }
    if (A > C) theta += 90.0;
  } 
	else if (A > C) theta = 90.0;
  
  // try to keep "a" > "b"
  if (a < b) { tmp = a; a = b; b = tmp; theta -= 90.0; }
  // try to keep "theta" = -45 ... 135 degrees
  if (theta < -45.0) theta += 180.0;
  if (theta > 135.0) theta -= 180.0;
  
  // ellipse coefficients
  ellipse.ResizeTo(5);
  ellipse[0] = x0; // ellipse's "x" center
  ellipse[1] = y0; // ellipse's "y" center
  ellipse[2] = a; // ellipse's "semimajor" axis along "x"
  ellipse[3] = b; // ellipse's "semiminor" axis along "y"
  ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)
  
  //return ellipse;
	return 0;

}//close ConicToParametric()



void Contour::ComputeCentroidDistanceFD(){

	//Reset lists
	CentroidDistanceModFDs.clear();
	int N= (int)m_Points.size();
	int n= N;//computing all Fourier descriptors (truncate to 10 in case)

	//Compute centroid distance function r(t)= sqrt( (x-xc)^2 + (y-yc)^2 )
	std::vector<double> r;
	double Xc= Centroid.X();
	double Yc= Centroid.Y();
	for(int i=0;i<N;i++){
		double X= m_Points[i].X();
		double Y= m_Points[i].Y();
		double dist= sqrt( pow(X-Xc,2) + pow(Y-Yc,2) );
		r.push_back(dist);
	}//end loop contour points
	
	//Compute the Discrete Fourier Transform of r	
	std::vector< std::complex<double> > Fn;
	for(int i=0; i<n; i++) {//loop over n
		Fn.push_back( std::complex<double>(0.,0.) );
		//int s= -floor(N/2.) + i;
		int s= i;

		for(int j=0;j<N;j++) {//loop over data size
			int k= j;
			double arg= 2.*TMath::Pi()*s*k/N;
			std::complex<double> prod= std::polar(1.,-arg);
			Fn[i]+= r[j]*prod;
		}//end loop data size
		Fn[i]/= N;
	}//end loop n

	for(unsigned int i=0; i<Fn.size(); i++) {
		//Normalize FD dividing by first FD
		Fn[i]/= Fn[1];

		//Compute modulus
		double FDMod= std::abs(Fn[i]); 
		CentroidDistanceModFDs.push_back(FDMod);
	}//end loop n
		
	HasCentroidDistanceFDPars= true;

}//close Contour::ComputeCentroidDistanceFD()



void Contour::ComputeFourierDescriptors(){

	//Reset lists
	//FDs.clear();
	RealFDs.clear();
	ImagFDs.clear();
	ModFDs.clear();
	int N= (int)m_Points.size();
	int n= N;//computing all Fourier descriptors (truncate to 10 in case)

	//Put contour point in a vector of complex numbers U
	std::vector< std::complex<double> > U= GetComplexPointRepresentation(true);//shift in centroid for convenience
		
	//Compute the Discrete Fourier Transform of contour	complex points	
	std::vector< std::complex<double> > Fn= MathUtils::DFTShifted(U,n);
		
	//Compute the Fourier descriptors 
	//- translational invariance ==> set F(0)=0
	//- scale invariance ==> set F(k)= F(k)/|F(1)|
	//- rotational invariance ==>?
	//- invariance against shifting of the initial point ==> set F(k)= F(k) exp(-i phi(1) k) 
	double phi_1= std::arg(Fn[1+floor(N/2.)]);//in radians
	double phi_minus1= std::arg(Fn[-1+floor(N/2.)]);//in radians
	double mod_1= std::abs(Fn[1+floor(N/2.)]);
	//double scaleInvariance= mod_1;
	
	for(unsigned int k=0;k<Fn.size();k++){
		int index= -floor(N/2.) + k;
		//std::complex<double> shiftInvariance= std::polar(1.,-phi_1*index);
		//if(index==0) Fn[k]= std::complex<double>(0.,0.);
		//else Fn[k]*= shiftInvariance/scaleInvariance;

		double Fn_mod= std::abs(Fn[k]);
		double Fn_phase= std::arg(Fn[k]);
		
		std::complex<double> Fn_norm= std::complex<double>(0.,0.);
		if(index!=0) Fn_norm= std::polar(Fn_mod/mod_1,Fn_phase-0.5*(phi_minus1+phi_1)+index*0.5*(phi_minus1-phi_1) );
		 
		//if(index==0) Fn[k]= std::complex<double>(0.,0.);
		//else Fn[k]/= scaleInvariance;
		//double FDMod= std::abs(Fn[k]);
		double FDMod= std::abs(Fn_norm);
	
		//FDs.push_back(Fn_norm);
		RealFDs.push_back(std::real(Fn_norm));
		ImagFDs.push_back(std::imag(Fn_norm));

		ModFDs.push_back(FDMod);
		DEBUG_LOG("FD no. "<<k<<" scale="<<index<<" FD="<<std::real(Fn[k])<<" + i "<<std::imag(Fn[k])<<" |FD|="<<FDMod);
	}//end loop fourier coeff

	HasFDPars= true;

}//close ComputeFourierDescriptors()


void Contour::ComputeBendingEnergy(){

	int N= (int)m_Points.size();
	
	//Put contour point in a vector of complex numbers ut
	std::vector< std::complex<double> > ut= GetComplexPointRepresentation(true);//translate points in centroid coordinate system
			
	//Compute ut'= IDST(i x 2pi x s x Us)	
	std::vector< std::complex<double> > Us= MathUtils::DFT(ut,N);
	std::vector< std::complex<double> > Us_firstDeriv; 
	for(int i=0;i<N;i++){
		int s= -floor(N/2.) + i;
		double arg= 2*TMath::Pi()*s;
		std::complex<double> z(0,arg); 
		std::complex<double> thisUs_firstDeriv= z*Us[i];	
		Us_firstDeriv.push_back(thisUs_firstDeriv);
	}//end loop points

	std::vector< std::complex<double> > ut_firstDeriv= MathUtils::IDFT(Us_firstDeriv,N);
	double L= 0;
	for(unsigned int i=0;i<ut_firstDeriv.size();i++) L+= std::abs(ut_firstDeriv[i]);	
	L*= 2.*TMath::Pi()/N;
		

	//Compute average bending energy
	const int NSCALES= 18;
	double SigmaScales[]= {
		1.e-7,5.e-7,1.e-6,5.e-6,1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2,5.e-2,0.1,0.5,1,2,5,10
	};
	BendingEnergies.clear();

	for(int k=0;k<NSCALES;k++){
		double smoothPar= SigmaScales[k];

		//Compute curvature at each contour point		
		std::vector<double> Curvature= MathUtils::GetContourCurvature(ut,smoothPar);

		//Compute average bending energy = L^2/N x sum(curv^2)
		double BendingEnergy= 0;
		for(unsigned int i=0;i<Curvature.size();i++){
			BendingEnergy+= pow(Curvature[i],2);
		}//end loop points
		BendingEnergy*= pow(L,2)/N;
		//BendingEnergy/= (double)N;
		BendingEnergies.push_back(BendingEnergy);
		DEBUG_LOG("Scale no. "<<k<<"="<<smoothPar<<" BE="<<BendingEnergy);
	}//end loop scales

	HasBEPars= true;

}//close ComputeBendingEnergy()

}//close namespace
