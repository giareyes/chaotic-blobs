#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "SETTINGS.h"
#include "MATERIAL.h"
#include "EXTRAFUNCTIONS.h"
#include <vector>

// TRIANGLE vertex ordering is COUNTER CLOCKWISE
//
//          v1
//          o
//         / \
//        /   \
//    e1 /     \ e0
//      /       \
//     /         \
//    /           \
//   o-------------o
// v2      e2      v0
//
class TRIANGLE
{
public:
  TRIANGLE(MATERIAL* material, const std::vector<VEC2*>& vertices);

  MATRIX2 computeF() const;

  VEC2 vertexAverage();
  VEC2* vertex(int i) { return _vertices[i]; };
  const VEC2& vertex(int i) const { return *_vertices[i]; };

  MATRIX computeForceJacobian();

  VECTOR computeForceVector();

  VECTOR precomputedLinearCoef();

  // compute rest area of this triangle
  Real restArea() const;
  Real area() const;

private:
  MATRIX pFpuVectorized();

  std::vector<MATRIX> pFpu();

  VEC2* _vertices[3];
  VEC2 _restPose[3];

  MATRIX2 _Dm;
  MATRIX2 _DmInverse;
  MATRIX6 _linearCoef;

  // material model
  MATERIAL* _material;

  // vector interface for getting and setting positions
  // (only used by unit test)
  VECTOR getDisplacements() const;
  void setDisplacements(const VECTOR& u);
};

#endif
