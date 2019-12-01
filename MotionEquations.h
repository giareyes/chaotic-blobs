#include "SETTINGS.h"
#include "TRIANGLE.h"
#include "MATERIAL.h"
#include <vector>

// regular equation of motion
MATRIX massMatrix();

// D(u, u') = (alpha*M + beta*K(u))u'
MATRIX dampingForce();

// K 
MATRIX internalForce();

// external force is calculated by STVK

// reduced equation of motion
MATRIX reducedDampingForce();

MATRIX reducedInternalForce();
