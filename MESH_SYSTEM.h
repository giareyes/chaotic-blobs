#ifndef MESH_SYSTEM_H
#define MESH_SYSTEM_H

#include "SETTINGS.h"
#include "TRIANGLE.h"
#include "MATERIAL.h"
#include "TRIANGLE_MESH.h"
#include "WALL.h"
#include <vector>
#include <map>

using namespace std;

class MESH_SYSTEM {
public:
  MESH_SYSTEM(const Real poissonsRatio = 0.3, const Real youngsModulus = 1e6, int numBlobs = 2);
  ~MESH_SYSTEM();

  void allBodyForce(const VEC2& bodyForce);
  void allStep(int scene, const Real deform);
  void allMove(float dt, const VEC2& outerForce);
  void allBuild(int sceneNum);
  void addWalls(const WALL& wall);

  

  std::vector<TRIANGLE_MESH> blobs;
};

#endif
