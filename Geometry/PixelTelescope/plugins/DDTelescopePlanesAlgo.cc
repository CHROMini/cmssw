//////////////////////////////////////////////////////////////////////////////////////////
// DDTelescopePlanesAlgo
// Description:  Places pixel telescope planes inside a given telescope arm.
// The planes can be tilted (rotation around CMS_X) then skewed (rotation around CMS_Y).
// The planes are all centered in (CMS_X,CMS_Y) = (0,0) and shifted along CMS_Z by deltaZ.
// Author: Gabrielle Hugo
//////////////////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DetectorDescription/Core/interface/DDCurrentNamespace.h"
#include "DetectorDescription/Core/interface/DDSplit.h"
#include "Geometry/PixelTelescope/plugins/DDTelescopePlanesAlgo.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"
#include "DetectorDescription/Core/interface/DDTransform.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"


DDTelescopePlanesAlgo::DDTelescopePlanesAlgo() {
  LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo info: Creating an instance";
}


DDTelescopePlanesAlgo::~DDTelescopePlanesAlgo() {}


void DDTelescopePlanesAlgo::initialize(const DDNumericArguments & nArgs,
				  const DDVectorArguments & vArgs,
				  const DDMapArguments & ,
				  const DDStringArguments & sArgs,
				  const DDStringVectorArguments & ) {

  n             = int(nArgs["N"]);
  tiltAngle     = nArgs["tiltAngle"];
  skewAngle     = nArgs["skewAngle"];
  deltaZ        = nArgs["deltaZ"];
  
  LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo debug: Parameters for position"
			  << "ing:: n " << n << " telescope planes with deltaZ "
			  << deltaZ << ", tiltAngle "
			  << tiltAngle/CLHEP::deg << ", skew angle " 
			  << skewAngle/CLHEP::deg;

  idNameSpace = DDCurrentNamespace::ns();
  childName   = sArgs["ChildName"];

  DDName parentName = parent().name();
  LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo debug: Parent " << parentName
			  << "\tChild " << childName << " NameSpace "
			  << idNameSpace;
}


void DDTelescopePlanesAlgo::execute(DDCompactView& cpv) {

  DDRotation prepaRot, tiltRot, skewRot, globalRot;                      // Default-constructed to identity in SO(3).
  DDRotationMatrix prepaMatrix, tiltMatrix, skewMatrix, globalRotMatrix; // Default-constructed to identity matrix.
  std::string rotstr = "RTelescopePlanesAlgo";

  // prepaMatrix calculus
  // Rotation around CMS_X of 90°, in counter-trigonometric sense.
  // This is because the Phase 1 modules frame of reference is the following: width along X, length along Z, thickness along Y.
  // TBM towards Y+.
  // After rotation, the Phase 1 module is placed with: width along X, length along Y, thickness along Z.
  // TBM towards Z-.
  std::string prepaRotstr = rotstr + "Prepa";
  prepaRot = DDRotation(DDName(prepaRotstr, idNameSpace));
  if (!prepaRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new rotation: " << prepaRotstr
			    << "\t90., 0., "
			    << "180., 0., "
			    << "90., 90.";
    prepaRot = DDrot(DDName(prepaRotstr, idNameSpace), 
		     90.*CLHEP::deg, 0., 180.*CLHEP::deg, 0.*CLHEP::deg, 90.*CLHEP::deg, 90.*CLHEP::deg);
  }
  prepaMatrix = *prepaRot.matrix();  // matrix of prepaRot

  // tiltMatrix calculus
  // Rotation around CMS_X. 
  // tiltAngle is counted in the counter-trigonometric sense. tiltAngle = 0 on (XY) plane. 
  // titlAngle ‎∈ [0° 90°].
  std::string tiltRotstr = rotstr + "Tilt" + std::to_string(tiltAngle/CLHEP::deg);
  tiltRot = DDRotation(DDName(tiltRotstr, idNameSpace));
  if (!tiltRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new rotation: " << tiltRotstr
			    << "\t90., 0., "
			    << 90. + tiltAngle/CLHEP::deg << ", 90., "
			    << tiltAngle/CLHEP::deg << ", 90.";
    tiltRot = DDrot(DDName(tiltRotstr, idNameSpace), 
		    90.*CLHEP::deg, 0., 90.*CLHEP::deg + tiltAngle, 90.*CLHEP::deg, tiltAngle, 90.*CLHEP::deg);
    // trigonometric sense:
    // tiltRot = DDrot(DDName(tiltRotstr, idNameSpace), 
    // 90.*CLHEP::deg, 0., 90.*CLHEP::deg - tiltAngle, 90.*CLHEP::deg, tiltAngle, 270.*CLHEP::deg);
  }
  tiltMatrix = *tiltRot.matrix();   // matrix of tiltRot
  tiltMatrix *= prepaMatrix;        // matrix of (tiltRot ◦ prepaRot)

  // skewMatrix calculus
  // Rotation around CMS_Y. 
  // skewAngle is counted in the trigonometric sense. skewAngle = 0 on (XY) plane. 
  // skewAngle ‎∈ [0° 90°].
  std::string skewRotstr = rotstr + "Skew" + std::to_string(skewAngle/CLHEP::deg);
  skewRot = DDRotation(DDName(skewRotstr, idNameSpace));
  if (!skewRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new rotation: " << skewRotstr
			    << "\t" << 90. + skewAngle/CLHEP::deg << ", 0., "
			    << "90., 90., "
			    << skewAngle/CLHEP::deg << ", 0.";
    skewRot = DDrot(DDName(skewRotstr, idNameSpace), 
		    90.*CLHEP::deg + skewAngle, 0., 90.*CLHEP::deg, 90.*CLHEP::deg, skewAngle, 0.);
    // counter-trigonometric sense:
    // skewRot = DDrot(DDName(skewRotstr, idNameSpace), 
    // 90.*CLHEP::deg - skewAngle, 0., 90.*CLHEP::deg, 90.*CLHEP::deg, skewAngle, 180.*CLHEP::deg);
  }
  skewMatrix = *skewRot.matrix();   // matrix of skewRot
  skewMatrix *= tiltMatrix;         // matrix of (skewRot ◦ tiltRot ◦ prepaRot)

  // globalRot def
  std::string globalRotstr = rotstr + "Global";
  globalRot = DDRotation(DDName(globalRotstr, idNameSpace));
  if (!globalRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new "
			    << "rotation: " << globalRotstr;
    globalRotMatrix = skewMatrix;   
    // Can finally create globalRot. globalRot = skewRot ◦ tiltRot ◦ prepaRot
    globalRot = DDrot(DDName(globalRotstr, idNameSpace), new DDRotationMatrix(globalRotMatrix)); 
  }

  // Loops on all n telescope planes
  DDName mother = parent().name();  // telescope arm
  DDName child(DDSplit(childName).first, DDSplit(childName).second);
  int    copy   = 1;

  for (int i=0; i<n; i++) {
  
    // translation def
    double xpos = 0.;
    double ypos = 0.;
    double zpos = (-(n-1.)/2. + i) * deltaZ;
    DDTranslation tran(xpos, ypos, zpos);
  
    // Positions child with respect to parent
    cpv.position(child, mother, copy, tran, globalRot); // Rotate child with globalRot, then translate it with tran
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test " << child << " number "
			    << copy << " positioned in " << mother << " at "
			    << tran  << " with " << globalRot;

    copy += 1;
  }
}
