#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include "SETTINGS.h"
#include "TRIANGLE.h"
#include "MATERIAL.h"
#include <vector>
#include <map>

// TRIANGLE vertex ordering is CLOCKWISE
class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH(const Real poissonsRatio = 0.3, const Real youngsModulus = 1e6);
  ~TRIANGLE_MESH();

  // build the different kinds of tests
  void buildBlob(const Real xPos);

  bool stepQuasistatic();

  // Euler's equation of motion
  void stepMotion(float dt, const VEC2& outerForce);

  // regular equation of motion
  void setMassMatrix();

  // void setVelocity();

  // D(u, u') = (alpha*M + beta*K(u))u'
  MATRIX dampingForce();

  // advance the constrained nodes for the stretch test
  void stepStretchTest(const Real stretch);
  void stretch2(const Real stretch);
  void stepSquashTest(const Real squash);
  void stepShearTest(const Real stretch);
  void addBodyForce(const VEC2& bodyForce);

  MATERIAL* material() { return _material; };
  VECTOR& fExternal() { return _fExternal; };

  // get a specific triangle
  TRIANGLE& getTriangle(const int index) { return _triangles[index]; };

  const int DOFs() { return _DOFs; };
  const std::vector<TRIANGLE>& triangles() { return _triangles; };
  const std::vector<VEC2>& vertices() { return _vertices; };
  const std::vector<int>& constrainedVertices() { return _constrainedVertices; };

private:
  // scatter displacement u to the vertices
  void uScatter();

  // gather displacement u from the vertices
  void uGather();

  void computeMaterialForces();

  // rebuild the vertex-to-index lookup
  void computeVertexToIndexTable();

  void computeStiffnessMatrix(MATRIX& K);

  // how many degrees of freedom are there?
  int _DOFs;

  // the displacement vector
  VECTOR _u;

  // the force vector
  VECTOR _f;
  VECTOR _fExternal;

  // mass matrix
  MATRIX _mass;
  VECTOR _velocity;

  // the geometry
  std::vector<VEC2>   _vertices;
  std::vector<VEC2>   _restVertices;
  std::vector<TRIANGLE> _triangles;

  // these vertices are constrained
  std::vector<int> _constrainedVertices;

  // these vertices are unconstrained
  std::vector<int> _unconstrainedVertices;

  // back-index from a vertex to its position in _u
  std::map<VEC2*, int> _vertexToIndex;

  // the material
  MATERIAL* _material;
};

#endif
