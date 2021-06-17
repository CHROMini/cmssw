#include "DQMTelescope/PixelTelescope/plugins/TelescopeTracks.h"
#include "TMultiGraph.h"
#include "TList.h"
#include "TGraph2DErrors.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "Math/Vector3D.h"
#include "TMath.h"

///////////////////////////////////////////////////////
/////////// external function definition //////////////
/////////// needed for the fit           //////////////
///////////////////////////////////////////////////////

using namespace ROOT::Math;

// calculate distance line-point 
// not used anymore...
double distance2(double x,double y,double z, double *p) { 
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
 
   XYZVector xp(x,y,z); 
   XYZVector x0(p[0], p[2], p[4]); 
   XYZVector x1(p[0] + p[1], p[2] + p[3], p[4] + p[5]); 
   XYZVector u = (x1-x0).Unit(); 
   double d2 = ((xp-x0).Cross(u)) .Mag2();
   
   //std::cout << "d2 ******* " << d2 << std::endl;
   return d2; 
}
	
// function to be minimized 
void SumDistance2(int &, double *, double & sum, double * par, int ) { 
//     TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
     TMultiGraph * gr = dynamic_cast<TMultiGraph*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
     TList* tlist_gr = gr->GetListOfGraphs();
     TObjOptLink *lnk = (TObjOptLink*)tlist_gr->FirstLink();
     TGraphErrors *obj = (TGraphErrors*) lnk->GetObject();;
     double * x = obj->GetX();
     double * y = obj->GetY();
     double * xerr = obj->GetEX();
     double * yerr = obj->GetEY();
     int npoints = obj->GetN();
     lnk = (TObjOptLink*)lnk->Next();
     TGraphErrors *obj2 = (TGraphErrors*) lnk->GetObject();;
     double * aparam= obj2->GetX();
     double * bparam= obj2->GetY();
     double * cparam= obj2->GetEX();
     double * dparam= obj2->GetEY();
     lnk = (TObjOptLink*)lnk->Next();
     TGraphErrors *obj3 = (TGraphErrors*) lnk->GetObject();;
     double * x0param= obj3->GetX();
     double * y0param= obj3->GetY();
     double * z0param= obj3->GetEX();
     lnk = (TObjOptLink*)lnk->Next();
     TGraphErrors *obj4 = (TGraphErrors*) lnk->GetObject();;
     double * thetaparam= obj4->GetX();
     double * phiparam= obj4->GetY();

//     std::cout << " in SumDistance " << x[0] << "  "  << y[0]  << "  " << aparam[0] << "  " << bparam[0] <<"  " << cparam[0] <<"  " << dparam[0] <<std::endl; 
    
     sum = 0;
     double pi=acos(-1);
     for (int i  = 0; i < npoints; ++i) { 
        // reverse matrices
        TMatrixD *rotationInvX_0 = new TMatrixD(3, 3);
        TMatrixD *rotationInvY_0 = new TMatrixD(3, 3);
        TArrayD dataIRx0(9);
        TArrayD dataIRy0(9);
        /*
        double costheta2 = cos(30*pi/180.);
        double sintheta2 = sin(30*pi/180.);
        double cosphi2 = cos(-20*pi/180.);
        double sinphi2 = sin(-20*pi/180.);
        */
        double costheta2 = cos(-1.*thetaparam[i]*pi/180.);
        double sintheta2 = sin(-1.*thetaparam[i]*pi/180.);
        double cosphi2 = cos(-1.*phiparam[i]*pi/180.);
        double sinphi2 = sin(-1.*phiparam[i]*pi/180.);
        dataIRx0[0] = 1; dataIRx0[1] = 0;         dataIRx0[2] = 0;
        dataIRx0[3] = 0; dataIRx0[4] = costheta2; dataIRx0[5] = -1*sintheta2;
        dataIRx0[6] = 0; dataIRx0[7] = sintheta2; dataIRx0[8] = costheta2;
        rotationInvX_0->SetMatrixArray(dataIRx0.GetArray());
        dataIRy0[0] = cosphi2; dataIRy0[1] = 0; dataIRy0[2] = sinphi2;
        dataIRy0[3] = 0;       dataIRy0[4] = 1; dataIRy0[5] = 0;
        dataIRy0[6] = -1.*sinphi2; dataIRy0[7] = 0; dataIRy0[8] = cosphi2;
        rotationInvY_0->SetMatrixArray(dataIRy0.GetArray());

//     	double d = distance2(x[i],y[i],z[i],par); 

        // intersection of the track with the telescope plane
        double tval = -1* (aparam[i]*par[0] + bparam[i]*par[2] + cparam[i]*par[4] + dparam[i])/(aparam[i]*par[1] + bparam[i]*par[3] + cparam[i]*par[5]);
        double globx = par[0] + par[1]*tval;
        double globy = par[2] + par[3]*tval;
        double globz = par[4] + par[5]*tval;
        // rotation of the intersection position to consider coordonates (x,y) within the telescope plane
        // after a translation to z=0

        TArrayD dataPlaneVec(3);
        dataPlaneVec[0] = globx-x0param[i];
        dataPlaneVec[1] = globy-y0param[i];
        dataPlaneVec[2] = globz-z0param[i];
        TMatrixD planeVec0(3, 1);
        planeVec0.SetMatrixArray(dataPlaneVec.GetArray());
        planeVec0 = (*rotationInvX_0)*(*rotationInvY_0)*planeVec0;
        //TVector3 localpos(TMatrixDColumn(planeVec, 0)[0], TMatrixDColumn(planeVec, 0)[1], TMatrixDColumn(planeVec, 0)[2]);
        double xtemp=TMatrixDColumn(planeVec0, 0)[0];
        double ytemp=TMatrixDColumn(planeVec0, 0)[1];

        // distance in the telescope plane
//      double d = pow((x[i]-xtemp)*(x[i]-xtemp) + (y[i]-ytemp)*(y[i]-ytemp),0.5);
        double d = (x[i]-xtemp)*(x[i]-xtemp)/(xerr[i]*xerr[i]) + (y[i]-ytemp)*(y[i]-ytemp)/(yerr[i]*yerr[i]);
	sum += d;
	//if (first) std::cout << "point " << i << "\t" << x[i] << "\t"  << y[i] << "\t" << z[i] << "\t"  << std::sqrt(d) << std::endl; 
        delete rotationInvX_0;
        delete rotationInvY_0;
     }
     // how to determine correctly the number of degrees of freedom ?
     sum /=npoints;

     //first = false;
     //std::cout << "sum " << sum << std::endl;
}

///////////////////////////////////////////////////////////
/////////// def of members of Telescopetrack //////////////
///////////////////////////////////////////////////////////

void TelescopeTracks::setTrackParam(int ipara, double value){

  if(ipara == 0)  p0_ = value;
  else if(ipara == 1)  p1_= value;
  else if(ipara == 2)  p2_= value;
  else if(ipara == 3)  p3_= value;
  else if(ipara == 4)  p4_= value;
  else if(ipara == 5)  p5_= value;

}

void TelescopeTracks::addCluster(SiPixelCluster theCluster){ clusterList_.push_back(theCluster);}

void TelescopeTracks::cleanClusterList(){ clusterList_.clear(); }

std::vector<SiPixelCluster> TelescopeTracks::getclusterList(){ return clusterList_;}

void     	      TelescopeTracks::addGlobalPoint(TVector3 gp ){globalPoints.push_back(gp); }
void 		      TelescopeTracks::cleanGlobalPoints(){ globalPoints.clear(); }
std::vector<TVector3> TelescopeTracks::getGlobalPoints(){ return globalPoints;} 
    
void        	      TelescopeTracks::addGlobalPointErr(TVector3 gp ){globalPoints_err.push_back(gp); }
void                  TelescopeTracks::cleanGlobalPointsErr(){ globalPoints_err.clear(); }
std::vector<TVector3> TelescopeTracks::getGlobalPointsErr(){ return globalPoints_err;}
    
void     	      TelescopeTracks::addPseudoLocalPoint(TVector3 gp ){pseudolocalPoints.push_back(gp); }
void 		      TelescopeTracks::cleanPseudoLocalPoints(){ pseudolocalPoints.clear(); }
std::vector<TVector3> TelescopeTracks::getPseudoLocalPoints(){ return pseudolocalPoints;}
     
void        	      TelescopeTracks::addPseudoLocalPointErr(TVector3 gp ){pseudolocalPoints_err.push_back(gp); }
void                  TelescopeTracks::cleanPseudoLocalPointsErr(){ pseudolocalPoints_err.clear(); }
std::vector<TVector3> TelescopeTracks::getPseudoLocalPointsErr(){ return pseudolocalPoints_err;}
    

double TelescopeTracks::getParameter(int ipara){

  if(ipara == 0) return p0_;
  else if(ipara == 1) return p1_;
  else if(ipara == 2) return p2_;
  else if(ipara == 3) return p3_;
  else if(ipara == 4) return p4_;
  else if(ipara == 5) return p5_;
  else return -10000;
  
}

void TelescopeTracks::intersection(double *planeq, double *p, double &x, double &y, double &z) { 
  // a parameteric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
  //
  //

  double a=planeq[0];
  double b=planeq[1];
  double c=planeq[2];
  double d=planeq[3];

  double t = -1* (a*p[0] + b*p[2] + c*p[4] + d)/(a*p[1] + b*p[3] + c*p[5]);
  x = p[0] + p[1]*t;
  y = p[2] + p[3]*t;
  z = p[4] + p[5]*t;  
}

void TelescopeTracks::line(double t, double *p, double &x, double &y, double &z) { 
  // a parameteric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
  x = p[0] + p[1]*t;
  y = p[2] + p[3]*t;
  z = p[4] + p[5]*t;  
}

void TelescopeTracks::fitTrack(){
  
  std::vector<double> track;
  
  const int nCluster = pseudolocalPoints.size();

  TMultiGraph* gr = new TMultiGraph();
  TGraphErrors* gr1pos = new TGraphErrors();
  TGraphErrors* grab= new TGraphErrors();
  TGraphErrors* grcd= new TGraphErrors();
  TGraphErrors* gref= new TGraphErrors();
  
//  std::cout << "debug fitTrack " << nCluster << std::endl;   
  for(int icls=0; icls<nCluster; icls++){
    gr1pos->SetPoint(icls, pseudolocalPoints[icls].X(),     pseudolocalPoints[icls].Y());
    gr1pos->SetPointError(icls, pseudolocalPoints_err[icls].X(), pseudolocalPoints_err[icls].Y());
    grab->SetPoint(icls,param_ass_a_[icls],param_ass_b_[icls]);
    grab->SetPointError(icls,param_ass_c_[icls],param_ass_d_[icls]);
    grcd->SetPoint(icls,param_ass_x0_[icls],param_ass_y0_[icls]);
    grcd->SetPointError(icls,param_ass_z0_[icls],param_ass_z0_[icls]);
    gref->SetPoint(icls,param_ass_theta_[icls],param_ass_phi_[icls]);
//    std::cout << " x " << globalPoints[icls].X() << "   "  << pseudolocalPoints[icls].X() << std::endl; 
//    std::cout << " y " << globalPoints[icls].Y() << "   "  << pseudolocalPoints[icls].Y() << std::endl; 
//    std::cout << " abcd " << param_ass_a_[icls] << "   " << param_ass_b_[icls] << "   "  << param_ass_c_[icls] << "   "  << param_ass_d_[icls] << std::endl;
  }
  
  gr->Add(gr1pos);
  gr->Add(grab);
  gr->Add(grcd);
  gr->Add(gref);
  
  TVirtualFitter *min = TVirtualFitter::Fitter(0,4);

  min->SetObjectFit(gr);
  min->SetFCN(SumDistance2);
  
  Double_t arglist[10];
//  arglist[0] = 3;
//  min->ExecuteCommand("SET PRINT",arglist,-1);
  arglist[0] = -1;
  min->ExecuteCommand("SET PRINT",arglist,1);
  
  double pStart[6] = {1,1,1,1,1,1};
  min->SetParameter(0,"p0",pStart[0],0.01,0,0);
  min->SetParameter(1,"p1",pStart[1],0.01,0,0);
  min->SetParameter(2,"p2",pStart[2],0.01,0,0);
  min->SetParameter(3,"p3",pStart[3],0.01,0,0);
  min->SetParameter(4,"p4",pStart[4],0.01,0,0);
  min->SetParameter(5,"p5",pStart[5],0.01,0,0);
  arglist[0] = 10000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  min->ExecuteCommand("MIGRAD",arglist,2);
  
  int nvpar,nparx; 
  double amin,edm, errdef;
  min->GetStats(amin,edm,errdef,nvpar,nparx);
  int ndf = nCluster-nvpar;
  
  normChi2_ =amin/ndf;
  chi2_=amin;
  //min->PrintResults(1,amin);
  
  // get fit parameters
  //double parFit[4];
  //for (int i = 0; i <4; ++i) parFit[i] = min->GetParameter(i);
 
  p0_ = min->GetParameter(0);
  p1_ = min->GetParameter(1);
  p2_ = min->GetParameter(2);
  p3_ = min->GetParameter(3);
  p4_ = min->GetParameter(4);
  p5_ = min->GetParameter(5);

/*
  std::cout << "tout va bien? " << std::endl;
  std::cout << "grab "<< grab << std::endl;
  delete grab;
  std::cout << "grab "<< grab << std::endl;
  std::cout << "gr1pos "<< gr1pos << std::endl;
  delete gr1pos;
  std::cout << "gr1pos "<< gr1pos << std::endl;
  std::cout << "grcd "<< grcd << std::endl;
  delete grcd;
  std::cout << "grcd "<< grcd << std::endl;
  std::cout << "TgraphError delete ok ? " << std::endl;
  gr1pos= 0;
  std::cout << "gr1pos "<< gr1pos << std::endl;
  grab=0;
  grcd=0;
  std::cout << "Link TgraphError a 0 ok ? " << std::endl;
*/
  delete gr;
  gr=0;
}



