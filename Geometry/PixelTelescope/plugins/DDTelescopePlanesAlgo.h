#ifndef DD_TelescopePlanesAlgo_h
#define DD_TelescopePlanesAlgo_h

//////////////////////////////////////////////////////////////////////////////////////////
// DDTelescopePlanesAlgo
// Description:  Places telescope plane frame of reference inside telescope arm frame of reference.
// This is the composition of: * a preparatory rotation (transition between Phase 1 module frame of reference and what we need),
//                             * then a tilt (rotation around CMS_X), 
//                             * then a skew (rotation around CMS_Y),
//                             * then a translation.
// Author: Gabrielle Hugo
//////////////////////////////////////////////////////////////////////////////////////////

#include <map>
#include <string>
#include <vector>
#include "DetectorDescription/Core/interface/DDTypes.h"
#include "DetectorDescription/Core/interface/DDAlgorithm.h"

class DDTelescopePlanesAlgo : public DDAlgorithm {
 
public:
  DDTelescopePlanesAlgo(); 
  ~DDTelescopePlanesAlgo() override;
  
  void initialize(const DDNumericArguments & nArgs,
		  const DDVectorArguments & vArgs,
		  const DDMapArguments & mArgs,
		  const DDStringArguments & sArgs,
		  const DDStringVectorArguments & vsArgs) override;

  void execute(DDCompactView& cpv) override;

private:

  int           childIndex;              // Plane index within arm.
  double        tiltAngle;               // Plane rotation around CMS_X. Angle is counted in the counter-trigonometric sense.
                                         // Angle = 0 on (XY) plane. Must be in [0° 90°].
  double        skewAngle;               // Plane rotation around CMS_Y. Angle is counted in the trigonometric sense. 
                                         // Angle = 0 on (XY) plane. Must be in [0° 90°].  
  std::vector<double> planeTranslation;  // Plane translation within arm.

  std::string   idNameSpace;    // Namespace of this and ALL sub-parts.
  std::string   childName;      // Child name (ie, telescope plane name).
};

#endif
