#ifndef DD_TelescopePlanesAlgo_h
#define DD_TelescopePlanesAlgo_h

//////////////////////////////////////////////////////////////////////////////////////////
// DDTelescopePlanesAlgo
// Description:  Places pixel telescope planes inside a given telescope arm.
// The planes can be tilted (rotation around CMS_X) then skewed (rotation around CMS_Y).
// The planes are all centered in (CMS_X,CMS_Y) = (0,0) and shifted along CMS_Z by deltaZ.
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

  int           n;              // Number of telescope planes.
  double        tiltAngle;      // Rotation around CMS_X. Angle is counted in the counter-trigonometric sense. Angle = 0 on (XY) plane. Must be in [0° 90°].
  double        skewAngle;      // Rotation around CMS_Y. Angle is counted in the trigonometric sense. Angle = 0 on (XY) plane. Must be in [0° 90°].
  double        deltaZ;         // Distance in Z between the centers of 2 consecutive planes.

  std::string   idNameSpace;    // Namespace of this and ALL sub-parts.
  std::string   childName;      // Child name (ie, telescope plane name).
};

#endif
