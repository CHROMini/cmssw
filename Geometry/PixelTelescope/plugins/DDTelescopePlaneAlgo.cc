//////////////////////////////////////////////////////////////////////////////////////////
// DDTelescopePlaneAlgo
// Description:  Places telescope plane frame of reference inside telescope arm frame of reference.
// This is the composition of: * a preparatory rotation (transition between Phase 1 module frame of reference and what we need),
//                             * then a tilt (rotation around CMS_X), 
//                             * then a skew (rotation around CMS_Y),
//                             * then a translation.
// Author: Gabrielle Hugo
//////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DetectorDescription/Core/interface/DDCurrentNamespace.h"
#include "DetectorDescription/Core/interface/DDSplit.h"
#include "Geometry/PixelTelescope/plugins/DDTelescopePlaneAlgo.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"
#include "DetectorDescription/Core/interface/DDTransform.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"


DDTelescopePlaneAlgo::DDTelescopePlaneAlgo() {
  LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo info: Creating an instance";
}


DDTelescopePlaneAlgo::~DDTelescopePlaneAlgo() {}


void DDTelescopePlaneAlgo::initialize(const DDNumericArguments & nArgs,
				  const DDVectorArguments & vArgs,
				  const DDMapArguments & ,
				  const DDStringArguments & sArgs,
				  const DDStringVectorArguments & ) {

  childIndex        = int(nArgs["ChildIndex"]);
  tiltAngle         = nArgs["planeTiltAngle"];
  skewAngle         = nArgs["planeSkewAngle"];
  planeTranslation  = vArgs["planeTranslation"];
  
  LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo debug: Parameters for positioning"
			  << " telescope plane index = " << childIndex << " in arm."
			  << " tiltAngle = " << tiltAngle/CLHEP::deg
			  << ", skew angle = " << skewAngle/CLHEP::deg;
  
  if (planeTranslation.size() != 3) {
    edm::LogError("DDTelescopePlaneAlgo") << " ERROR - Parameter planeTranslation should be of size 3.";
  }
  else {
    LogDebug("DDTelescopePlaneAlgo") << ", translation in X = " << planeTranslation.at(0)
				      << ", translation in Y = " << planeTranslation.at(1)
				      << ", translation in Z = " << planeTranslation.at(2);
  }

  idNameSpace = DDCurrentNamespace::ns();
  childName   = sArgs["ChildName"];

  DDName parentName = parent().name();
  LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo debug: Parent " << parentName
			  << "\tChild " << childName << " NameSpace "
			  << idNameSpace;
}


void DDTelescopePlaneAlgo::execute(DDCompactView& cpv) {

  DDRotation prepaRot, tiltRot, skewRot, globalRot;                      // Default-constructed to identity in SO(3).
  DDRotationMatrix prepaMatrix, tiltMatrix, skewMatrix, globalRotMatrix; // Default-constructed to identity matrix.
  std::string rotstr = "RTelescopePlaneAlgo";

  // prepaMatrix calculus
  // Rotation around CMS_X of 90°, in counter-trigonometric sense.
  // This is because the Phase 1 modules frame of reference is the following: width along X, length along Z, thickness along Y.
  // TBM towards Y+.
  // After rotation, the Phase 1 module is placed with: width along X, length along Y, thickness along Z.
  // TBM towards Z-.
  std::string prepaRotstr = rotstr + "Prepa";
  prepaRot = DDRotation(DDName(prepaRotstr, idNameSpace));
  if (!prepaRot) {
    LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo test: Creating a new rotation: " << prepaRotstr
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
    LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo test: Creating a new rotation: " << tiltRotstr
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
    LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo test: Creating a new rotation: " << skewRotstr
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
    LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo test: Creating a new "
			    << "rotation: " << globalRotstr;
    globalRotMatrix = skewMatrix;   
    // Can finally create globalRot. globalRot = skewRot ◦ tiltRot ◦ prepaRot
    globalRot = DDrot(DDName(globalRotstr, idNameSpace), new DDRotationMatrix(globalRotMatrix)); 
  }

  // Places plane within telescope
  DDName mother = parent().name();  // telescope arm
  DDName child(DDSplit(childName).first, DDSplit(childName).second);  // telescope plane

  // translation def
  double xpos = 0;
  double ypos = 0;
  double zpos = 0;
  if (planeTranslation.size() != 3) {
    edm::LogError("DDTelescopePlaneAlgo") << " ERROR - Parameter planeTranslation should be of size 3.";
  }
  else {
    xpos = planeTranslation.at(0);
    ypos = planeTranslation.at(1);
    zpos = planeTranslation.at(2);
  }
  DDTranslation tran(xpos, ypos, zpos);
  
  // Positions child with respect to parent
  cpv.position(child, mother, childIndex, tran, globalRot); // Rotate child with globalRot, then translate it with tran
  LogDebug("DDTelescopePlaneAlgo") << "DDTelescopePlaneAlgo test " << child << " number "
			  << childIndex << " positioned in " << mother << " at "
			  << tran  << " with " << globalRot;
}
