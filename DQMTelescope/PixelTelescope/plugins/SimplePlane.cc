#include "DQMTelescope/PixelTelescope/plugins/SimplePlane.h"



void SimplePlane::doRotation(){

  if (isRotated == false)  {
   
   TArrayD dataPlaneVec(3);
   TMatrixD planeVec(3, 1);
   dataPlaneVec[0] = a;
   dataPlaneVec[1] = b;
   dataPlaneVec[2] = c;
   
   planeVec.SetMatrixArray(dataPlaneVec.GetArray());
   
   planeVec = (*rotationY_)*(*rotationX_)*planeVec;
   
   a=TMatrixDColumn(planeVec, 0)[0];
   b=TMatrixDColumn(planeVec, 0)[1];
   c=TMatrixDColumn(planeVec, 0)[2];
   d=-1.*a*Origin.X()-1.*b*Origin.Y()-1.*c*Origin.Z();
  
   isRotated = true;
  }
}
void SimplePlane::applyAlignment(){
          int detId= detID1_;
          double neworiginx=Origin.X();
          double neworiginy=Origin.Y();         
          double neworiginz=Origin.Z();         
   
          if (detId==353114116) {
            neworiginx-=-2.56917e-02;
            neworiginy-=-2.38370e-02;
          }
          else if (detId==353113092) {
            neworiginx-=-2.02377e-02;
            neworiginy-=-2.70933e-02;
          }
          if (detId==352851972) {
           neworiginx-=2.47805e-02;
           neworiginy-=-6.97738e-03;
          }
          else if (detId==352850948) {
           neworiginx-=1.85113e-02;
           neworiginy-=-9.17572e-03;
          }
          else if (detId==352589828) {
           neworiginx-=5.50173e-02;
           neworiginy-=-3.05621e-02;
          }
          else if (detId==352588804) {
           neworiginx-=6.73395e-02;
           neworiginy-=-3.06608e-02;
          }
          else if (detId==344200196) {
           neworiginx-=-2.77220e-02;
           neworiginy-=-9.31173e-02;
          }
          else if (detId==344201220) {
           neworiginx-=-2.24462e-02;
           neworiginy-=-8.72844e-02;
          }
          else if (detId==344462340) {
           neworiginx-=7.14737e-03;
           neworiginy-=-1.09290e-01;
          }
          else if (detId==344724484) {
           neworiginx-=-7.60673e-02;
           neworiginy-=-1.20672e-01;
          }
          else if (detId==344986628) {
           neworiginx-=7.71767e-02;
           neworiginy-=-1.20759e-01;
          }
          else if (detId==344987652) {
           neworiginx-=9.50454e-02;
           neworiginy-=-1.33751e-01;
          }
          Origin = TVector3 (neworiginx, neworiginy, neworiginz);

}

TVector3 SimplePlane::getLocalPointPosition(TVector3 globalPos){
  
  TArrayD dataPlaneVec(3);
  
  dataPlaneVec[0] = globalPos.X()-Origin.X();
  dataPlaneVec[1] = globalPos.Y()-Origin.Y();
  dataPlaneVec[2] = globalPos.Z()-Origin.Z();
  
  TMatrixD planeVec(3, 1);
  planeVec.SetMatrixArray(dataPlaneVec.GetArray());
  planeVec = (*rotationInvX_)*(*rotationInvY_)*planeVec;
   
  TVector3 localpos(TMatrixDColumn(planeVec, 0)[0], TMatrixDColumn(planeVec, 0)[1], TMatrixDColumn(planeVec, 0)[2]);
   
/*
  std::cout << "debut getLocalPointPosition() : rotationInvY_ " << std::endl;
  //rotationInvY_->Print();
  std::cout << "globalPos.X() " << globalPos.X() << " localpos.X() " << localpos.X() << std::endl;
  std::cout << "globalPos.Y() " << globalPos.Y() << " localpos.Y() " << localpos.Y() << std::endl;
  std::cout << "globalPos.Z() " << globalPos.Z() << " localpos.Z() " << localpos.Z() << std::endl;
*/

  return localpos;
  
}


TVector3 SimplePlane::getLocalTrackPosition(double *par){ 
  
  double t = -1* (a*par[0] + b*par[2] + c*par[4] + d)/(a*par[1] + b*par[3] + c*par[5]); 
  
  TVector3 globalPos(par[0] + t*par[1], par[2] + t*par[3], par[4] + t*par[5] );

  TArrayD dataPlaneVec(3);
  
  dataPlaneVec[0] = globalPos.X()-Origin.X();
  dataPlaneVec[1] = globalPos.Y()-Origin.Y();
  dataPlaneVec[2] = globalPos.Z()-Origin.Z();
  
  TMatrixD planeVec(3, 1);
  planeVec.SetMatrixArray(dataPlaneVec.GetArray());
  planeVec = (*rotationInvX_)*(*rotationInvY_)*planeVec;
   
  TVector3 localpos(TMatrixDColumn(planeVec, 0)[0], TMatrixDColumn(planeVec, 0)[1], TMatrixDColumn(planeVec, 0)[2]);
  
  return localpos;
  
  
}
