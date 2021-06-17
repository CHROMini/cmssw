#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TArrayD.h"



class SimplePlane {


  //to define the equation of the plane ax+by+cz+d = 0
  double a, b, c, d;
  
  //point taken as reference for the plane :
  TVector3 Origin;


  //angles of the plane:
  double theta_plane, phi_plane;
  bool isRotated;
  
  int layerNum_;
  int detID1_;
//  int detID2_;
  
  TMatrixD *rotationX_, *rotationY_;
  TMatrixD *rotationInvX_, *rotationInvY_;
  
public:
    //Default Constructor 
    SimplePlane(){
      
      isRotated = false;
      Origin = TVector3 (0, 0, 0);
      a=0;
      b=0;
      c=0;
      d=0;
      theta_plane=0;
      phi_plane=0;
    }
    
    

    SimplePlane(double x, double y, double z, double theta, double phi,  int layerNum,  int detID1){
      
      isRotated = false;
      Origin = TVector3 (x, y, z);
      theta_plane=theta;
      phi_plane=phi;
      a=0;
      b=0;
      c=1;
      d=-1.*a*x-1.*b*y-1.*c*z;
      layerNum_ = layerNum;
      detID1_   = detID1;
      
      rotationX_ = new TMatrixD(3, 3);
      rotationY_ = new TMatrixD(3, 3);

      TArrayD dataRx(9);
      TArrayD dataRy(9);

      double pi=acos(-1);
      double costheta = cos(theta*pi/180.);
      double sintheta = sin(theta*pi/180.);
      double cosphi = cos(phi*pi/180.);
      double sinphi = sin(phi*pi/180.);

      //https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

      dataRx[0] = 1; dataRx[1] = 0;	   dataRx[2] = 0;
      dataRx[3] = 0; dataRx[4] = costheta; dataRx[5] = -1*sintheta;
      dataRx[6] = 0; dataRx[7] = sintheta; dataRx[8] = costheta;

      rotationX_->SetMatrixArray(dataRx.GetArray());

      dataRy[0] = cosphi; dataRy[1] = 0; dataRy[2] = sinphi;
      dataRy[3] = 0;	  dataRy[4] = 1; dataRy[5] = 0;
      dataRy[6] = -1.*sinphi; dataRy[7] = 0; dataRy[8] = cosphi;

      rotationY_->SetMatrixArray(dataRy.GetArray());

      rotationInvX_ = new TMatrixD(3, 3);
      rotationInvY_ = new TMatrixD(3, 3);
      TArrayD dataIRx(9);
      TArrayD dataIRy(9);

      double costheta2 = costheta;
      double sintheta2 = -1* sintheta;


      double cosphi2 = cosphi;
      double sinphi2 = -1*sinphi;


      dataIRx[0] = 1; dataIRx[1] = 0;	      dataIRx[2] = 0;
      dataIRx[3] = 0; dataIRx[4] = costheta2; dataIRx[5] = -1*sintheta2;
      dataIRx[6] = 0; dataIRx[7] = sintheta2; dataIRx[8] = costheta2;

      rotationInvX_->SetMatrixArray(dataIRx.GetArray());

      dataIRy[0] = cosphi2; dataIRy[1] = 0; dataIRy[2] = sinphi2;
      dataIRy[3] = 0;	    dataIRy[4] = 1; dataIRy[5] = 0;
      dataIRy[6] = -1.*sinphi2; dataIRy[7] = 0; dataIRy[8] = cosphi2;

      rotationInvY_->SetMatrixArray(dataIRy.GetArray());
      
    }
    
    int getLayerNum(){return layerNum_;};
    int getDetID1(){return detID1_;};
    double getParamA() {return a;}
    double getParamB() {return b;}
    double getParamC() {return c;}
    double getParamD() {return d;}
    double getOriginX() {return Origin.X();}
    double getOriginY() {return Origin.Y();}
    double getOriginZ() {return Origin.Z();}
    double getTheta() {return theta_plane;}
    double getPhi()   {return phi_plane;}
    
    void doRotation();
    void applyAlignment();
    
    TVector3 getLocalPointPosition(TVector3 globalPos);
    TVector3 getLocalTrackPosition(double *par);
    
} ;
  
