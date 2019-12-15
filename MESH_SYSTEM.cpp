#include "MESH_SYSTEM.h"
#include <iostream>

#include <float.h>
#include <random>


MESH_SYSTEM::MESH_SYSTEM(const Real poissonsRatio, const Real youngsModulus, int numBlobs) : _DOFs(0)
{
  for(int i = 0; i < numBlobs; i++)
  {
    blobs.push_back(TRIANGLE_MESH(poissonsRatio, youngsModulus));
  }
}
