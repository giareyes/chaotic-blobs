#include "TRIANGLE.h"

#include <iostream>
using namespace std;

TRIANGLE::TRIANGLE(MATERIAL* material, const vector<VEC2*>& vertices) :
  _material(material)
{
  assert(vertices.size() == 3);
  _vertices[0] = vertices[0];
  _vertices[1] = vertices[1];
  _vertices[2] = vertices[2];

  // store these as the rest pose, for now
  for (unsigned int x = 0; x < 3; x++)
    _restPose[x] = *_vertices[x];

    VEC2 e0 = _restPose[1] - _restPose[0];
    VEC2 e1 = _restPose[2] - _restPose[0];

    // put e0, e1 into Dm
    _Dm.col(0) = e0;
    _Dm.col(1) = e1;

    // calculate Dm inverse
    _DmInverse = _Dm.inverse().eval();

    Real lambda = _material->getLambda();
    Real mu = _material->getMu();
    //Real det = _Dm.determinant();
    //MATRIX2 _DmInverse = det*_DmInverse;

    // create the linear coefficient matrix below:
    _linearCoef.setZero();

    _linearCoef(0,0) = _DmInverse(0,0)*_DmInverse(0,0) + _DmInverse(1,0)*_DmInverse(1,0) + _DmInverse(0,1)*_DmInverse(0,1)
                        + _DmInverse(1,1)*_DmInverse(1,1) + 2*_DmInverse(0,0)*_DmInverse(1,0) + 2*_DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(1,1) = _linearCoef(0,0);

    _linearCoef(2,0) = -1*_DmInverse(0,0)*_DmInverse(0,0) - _DmInverse(0,1)*_DmInverse(0,1) - _DmInverse(0,0)*_DmInverse(1,0) - _DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(3,1) = _linearCoef(2,0);

    _linearCoef(4,0) = -1*_DmInverse(1,0)*_DmInverse(1,0) - _DmInverse(1,1)*_DmInverse(1,1) - _DmInverse(0,0)*_DmInverse(1,0) - _DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(5,1) = _linearCoef(4,0);

    _linearCoef(0,2) = -1*_DmInverse(0,0)*_DmInverse(0,0) - _DmInverse(0,1)*_DmInverse(0,1) - _DmInverse(0,0)*_DmInverse(1,0) - _DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(1,3) = _linearCoef(0,2);

    _linearCoef(2,2) = _DmInverse(0,0)*_DmInverse(0,0) + _DmInverse(0,1)*_DmInverse(0,1);
    _linearCoef(3,3) = _linearCoef(2,2);

    _linearCoef(4,2) = _DmInverse(0,0)*_DmInverse(1,0) + _DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(5,3) = _linearCoef(4,2);

    _linearCoef(0,4) = -1*_DmInverse(1,0)*_DmInverse(1,0) - _DmInverse(1,1)*_DmInverse(1,1) - _DmInverse(0,0)*_DmInverse(1,0) - _DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(1,5) = _linearCoef(0,4);

    _linearCoef(2,4) = _DmInverse(0,0)*_DmInverse(1,0) + _DmInverse(0,1)*_DmInverse(1,1);
    _linearCoef(3,5) = _linearCoef(2,4);

    _linearCoef(4,4) = _DmInverse(1,0)*_DmInverse(1,0) + _DmInverse(1,1)*_DmInverse(1,1);
    _linearCoef(5,5) = _linearCoef(4,4);


    // get restArea
    // get triangle normal in R^3
    VEC3 triangleNormal;
    VEC3 restPose3[3];
    for (int x = 0; x < 3; x++)
    {
      restPose3[x].setZero();
      for (int y = 0; y < 2; y++)
        restPose3[x][y] = _restPose[x][y];
    }
    triangleNormal = (restPose3[2] - restPose3[0]).cross(restPose3[1] - restPose3[0]);

    Real restArea = triangleNormal.norm() * 0.5;

    // multiply by (-4*mu - 4*lambda)*-1*restArea
    _linearCoef = restArea*(4*mu + 4*lambda)*_linearCoef;

}

///////////////////////////////////////////////////////////////////////
// get the deformation gradient at Xi
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::computeF() const
{
  MATRIX2 F; //force gradient; will be returned
  MATRIX2 Ds; // spatial matrix

  //deformed e0 and e1
  VEC2 de0 = *_vertices[1] - *_vertices[0];
  VEC2 de1 = *_vertices[2] - *_vertices[0];

  //put de0 and de1 into our spacial matrix
  Ds.col(0) = de0;
  Ds.col(1) = de1;

  //multiply Ds and Dm inverse to calcuate F
  F = Ds*_DmInverse;

  return F;
}

///////////////////////////////////////////////////////////////////////
// take the average of the vertices
///////////////////////////////////////////////////////////////////////
VEC2 TRIANGLE::vertexAverage()
{
  VEC2 final = (*_vertices[0]);
  for (int x = 1; x < 3; x++)
    final += (*_vertices[x]);

  return final * 1.0 / 3.0;
}

///////////////////////////////////////////////////////////////////////
// use a precomputed Linear Coefficient
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::precomputedLinearCoef()
{
  // figure out how to get lamda and mu
  VECTOR linearForce(6);
  VECTOR pos(6);
  for(int i = 0; i < 3; i++)
  {
    pos[2*i] = (*_vertices[i])[0];
    pos[2*i + 1] = (*_vertices[i])[1];
  }

  // linearForce = _linearCoef*getDisplacements();
  linearForce = _linearCoef*pos;
  return linearForce;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVector()
{
  VECTOR forceVector(6);
  MATRIX2 F = computeF();
  MATRIX2 pk1 = _material->PK1(F);
  MATRIX dfdu(4,6);

  //multiply pk1 by -1 and area as discussed in OH
  pk1 = -1*restArea()*pk1;

  //get vectorized df/du
  dfdu = pFpuVectorized();

  forceVector = dfdu.transpose()*vectorize(pk1) + precomputedLinearCoef();

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobian()
{
  MATRIX2 F = computeF();
  MATRIX jacobian(6,6);
  MATRIX dfdu(4,6);

  //get 4x4 matrix of derivative of pk1. This is matrix B from class
  MATRIX dpdf = _material->DPDF(F);

  //multiply dpdf by -1 and area as discussed in OH
  dpdf = -1*restArea()*dpdf;

  //get vectorized df/du. This is matrix A from class
  dfdu = pFpuVectorized();

  // d^2 psi/ dx^2 = A transpose * B * A
  jacobian = dfdu.transpose() * dpdf * dfdu;

  return jacobian;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
vector<MATRIX> TRIANGLE::pFpu()
{
  vector<MATRIX> result; //the end 3rd order tensor

  //compute dF/x0 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x0;
  x0.setZero();
  x0(0,0) = -1;
  x0(0,1) = -1;
  x0 = x0*_DmInverse;
  result.push_back(x0);

  //compute dF/x1 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x1;
  x1.setZero();
  x1(1,0) = -1;
  x1(1,1) = -1;
  x1 = x1*_DmInverse;
  result.push_back(x1);

  //compute dF/x2 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x2;
  x2.setZero();
  x2(0,0) = 1;
  x2 = x2*_DmInverse;
  result.push_back(x2);

  //compute dF/x3 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x3;
  x3.setZero();
  x3(1,0) = 1;
  x3 = x3*_DmInverse;
  result.push_back(x3);

  //compute dF/x4 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x4;
  x4.setZero();
  x4(0,1) = 1;
  x4 = x4*_DmInverse;
  result.push_back(x4);

  //compute dF/x5 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x5;
  x5.setZero();
  x5(1,1) = 1;
  x5 = x5*_DmInverse;
  result.push_back(x5);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpuVectorized()
{
  //take result from pFpu and vectorize it using function from EXTRAFUNCTIONS.cpp
  vector<MATRIX> pfpu = pFpu();
  MATRIX vectorized(4,6);
  for(int x = 0; x < 6; x++)
  {
    vectorized.col(x) = vectorize(pfpu[x]);
  }
  return vectorized;
}

///////////////////////////////////////////////////////////////////////
// compute rest area of this triangle
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::restArea() const
{
  // get triangle normal in R^3
  VEC3 triangleNormal;
  VEC3 restPose3[3];
  for (int x = 0; x < 3; x++)
  {
    restPose3[x].setZero();
    for (int y = 0; y < 2; y++)
      restPose3[x][y] = _restPose[x][y];
  }
  triangleNormal = (restPose3[2] - restPose3[0]).cross(restPose3[1] - restPose3[0]);
  return triangleNormal.norm() * 0.5;
}

///////////////////////////////////////////////////////////////////////
// compute deformed area of this triangle
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::area() const
{
  // get triangle normal in R^3
  VEC3 triangleNormal;
  VEC3 pose3[3];
  for (int x = 0; x < 3; x++)
  {
    pose3[x].setZero();
    for (int y = 0; y < 2; y++)
      pose3[x][y] = (*_vertices[x])[y];
  }
  triangleNormal = (pose3[2] - pose3[0]).cross(pose3[1] - pose3[0]);
  return triangleNormal.norm() * 0.5;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::getDisplacements() const
{
  VECTOR result(6);
  result[0] = (*_vertices[0])[0] - _restPose[0][0];
  result[1] = (*_vertices[0])[1] - _restPose[0][1];
  result[2] = (*_vertices[1])[0] - _restPose[1][0];
  result[3] = (*_vertices[1])[1] - _restPose[1][1];
  result[4] = (*_vertices[2])[0] - _restPose[2][0];
  result[5] = (*_vertices[2])[1] - _restPose[2][1];

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::setDisplacements(const VECTOR& u)
{
  (*_vertices[0])[0] = _restPose[0][0] + u[0];
  (*_vertices[0])[1] = _restPose[0][1] + u[1];
  (*_vertices[1])[0] = _restPose[1][0] + u[2];
  (*_vertices[1])[1] = _restPose[1][1] + u[3];
  (*_vertices[2])[0] = _restPose[2][0] + u[4];
  (*_vertices[2])[1] = _restPose[2][1] + u[5];
}
