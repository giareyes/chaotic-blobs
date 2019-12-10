#include "TRIANGLE_MESH.h"
#include "STVK.h"
#include <iostream>

#include <float.h>
#include <random>

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const Real poissonsRatio, const Real youngsModulus) : _DOFs(0)
{
  const Real E = youngsModulus;
  const Real nu = poissonsRatio;

  cout << " Young's Modulus: " << youngsModulus << endl;
  cout << " Poisson's Ratio: " << poissonsRatio << endl;

  const Real lambda = (0.25)*( E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)));
  const Real mu = (0.25)*(E / (2.0 * (1 + nu)));
  _material = new STVK(lambda, mu);
}

TRIANGLE_MESH::~TRIANGLE_MESH()
{
  delete _material;
}

void TRIANGLE_MESH::buildBlob(const Real xPos)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  // build bottom (constrained) vertices
  vector<VEC2> blobButt;
  for(int i = 0; i < 5; i++)
  {
    VEC2 v0(xPos + 0.1*i, 0);
    blobButt.push_back(v0);
    _vertices.push_back(v0);
    _restVertices.push_back(v0);
    _constrainedVertices.push_back(i);
  }

  // build tops to bottom triangles
  int vCount = 5;
  float scale = 5;
  float base = xPos + 0.05;
  for(int j = 4; j > 0; j--)
  {
    for(int i = 0; i < j; i++)
    {
      VEC2 v0(base + i*0.1, (sqrt(3.0) / 20.0)*(scale-j));
      _vertices.push_back(v0);
      _restVertices.push_back(v0);
      (j == 1)? _constrainedVertices.push_back(vCount + i) : _unconstrainedVertices.push_back(vCount + i);
    }
    vCount += j;
    base += 0.05;
  }

  int endofBaseVerts = vCount;

  // build left curved side vertices
  for(int i = 0; i < 3; i++)
  {
    float y;
    y = (i == 2)? (sqrt(3.0) / 20.0) + (sqrt(3.0) / 40.0) : (sqrt(3.0) / 40.0) + (i*(sqrt(3.0) / 20.0)/(i+1));
    VEC2 vo(xPos - 0.03 + (i+i)*0.01, y);
    _vertices.push_back(vo);
    _restVertices.push_back(vo);
    _unconstrainedVertices.push_back(vCount);
    vCount++;
  }

  VEC2 vo(xPos + 0.05, (sqrt(3.0) / 20.0)*2);
  _vertices.push_back(vo);
  _restVertices.push_back(vo);
  _unconstrainedVertices.push_back(vCount);
  vCount++;

  // build right curved side vertices
  for(int i = 0; i < 3; i++)
  {
    float y;
    y = (i == 2)? (sqrt(3.0) / 20.0) + (sqrt(3.0) / 40.0) : (sqrt(3.0) / 40.0) + (i*(sqrt(3.0) / 20.0)/(i+1));
    VEC2 vo2(xPos + 0.4 + 0.03 - (i+i)*0.01, y);
    _vertices.push_back(vo2);
    _restVertices.push_back(vo2);
    _unconstrainedVertices.push_back(vCount);
    vCount++;
  }

  VEC2 vo2(xPos + 0.35, (sqrt(3.0) / 20.0)*2);
  _vertices.push_back(vo2);
  _restVertices.push_back(vo2);
  _unconstrainedVertices.push_back(vCount);
  vCount++;

  // build the base triangles
  int oldCount = 0;
  int newCount = 5;
  for(int max = 4; max > 0; max--)
  {
  for(int i = 0; i < max; i++)
  {
      vector<VEC2*> v(4);
      v[0] = &_vertices[oldCount + i];
      v[1] = &_vertices[oldCount + i + 1];
      v[2] = &_vertices[newCount + i];
      if(i < max - 1)
        v[3] = &_vertices[newCount + 1 + i];

      //bottom
      vector<VEC2*> triangle;
      triangle.push_back(v[0]);
      triangle.push_back(v[1]);
      triangle.push_back(v[2]);
      _triangles.push_back(TRIANGLE(_material, triangle));

      //top
      if(i < max - 1)
      {
        vector<VEC2*> top;
        top.push_back(v[2]);
        top.push_back(v[3]);
        top.push_back(v[1]);
        _triangles.push_back(TRIANGLE(_material, top));
      }
    }
    oldCount = newCount;
    newCount += max;
  }

  // build curved vertices
  for(int i = 0; i < 2; i++)
  {
    vector<VEC2*> v(12);
    v[0] = &_vertices[5 + i*3];
    v[1] = &_vertices[0 + i*4];
    v[2] = &_vertices[endofBaseVerts + i*4];
    v[3] = &_vertices[endofBaseVerts + 1 + i*4];
    v[4] = &_vertices[endofBaseVerts + 2 + i*4];
    v[5] = &_vertices[endofBaseVerts + 3 + i*4];
    v[6] = &_vertices[9 + i*2];
    v[7] = &_vertices[endofBaseVerts - (3 - i)];

    vector<VEC2*> triangle;
    triangle.push_back(v[0]);
    triangle.push_back(v[1]);
    triangle.push_back(v[2]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[0]);
    triangle.push_back(v[2]);
    triangle.push_back(v[3]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[0]);
    triangle.push_back(v[3]);
    triangle.push_back(v[4]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[0]);
    triangle.push_back(v[6]);
    triangle.push_back(v[4]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[4]);
    triangle.push_back(v[6]);
    triangle.push_back(v[5]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[5]);
    triangle.push_back(v[6]);
    triangle.push_back(v[7]);
    _triangles.push_back(TRIANGLE(_material, triangle));
  }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  setBasisReduction();
  setMassMatrix();

  VECTOR zeros(_DOFs);
  VECTOR z2(_U.cols());

  z2.setZero();

  zeros.setZero();
  _u          = zeros;
  _q          = z2;
  _f          = zeros;
  _fExternal  = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}

///////////////////////////////////////////////////////////////////////
// rebuild the vertex-to-index lookup
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeVertexToIndexTable()
{
  _vertexToIndex.clear();
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    VEC2* vertex = &_vertices[_unconstrainedVertices[x]];
    _vertexToIndex[vertex] = 2 * x;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::uScatter()
{
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    _vertices[index][0] = _restVertices[index][0] + _u[2 * x];
    _vertices[index][1] = _restVertices[index][1] + _u[2 * x + 1];
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::uGather()
{
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    _u[2 * x]     = _vertices[index][0] - _restVertices[index][0];
    _u[2 * x + 1] = _vertices[index][1] - _restVertices[index][1];
  }
}

void TRIANGLE_MESH::setMassMatrix()
{
  // matrix of size 2Nx2N
  MATRIX M(_unconstrainedVertices.size()*2,_unconstrainedVertices.size()*2);

  //for now we will make every vertex has mass 1
  M.setIdentity();
  _mass = M;
}

void TRIANGLE_MESH::setBasisReduction()
{
  MATRIX U(34,9);
  MATRIX svddiag(34, 9);
  svddiag.setIdentity();
  // matrix of deformations
  U << 0.011446,  -0.002633,  -0.016704,   0.032098,  -0.009118,  -0.054865,  -0.054865,    0.088559,     0.067614,
     0.005271,  -0.001795,  -0.006527,   0.005844,  -0.015733,  -0.019406,  -0.019406,    0.006542,    -0.021152,
     0.004064,  -0.000470,  -0.017430,   0.022485,  -0.002303,  -0.085461,  -0.055144,    0.085461,     0.059172,
     0.009012,  -0.001448,  -0.000137,   0.003981,  -0.011465,  -0.000642,  -0.014926,    0.000642,    -0.017341,
    -0.004064,   0.000470,  -0.022485,   0.017430,   0.002303,  -0.085461,  -0.059079,    0.085461,     0.055184,
     0.009012,  -0.001448,   0.003981,  -0.000137,  -0.011465,   0.000642,  -0.017255,   -0.000642,    -0.014943,
    -0.011446,   0.002633,  -0.032098,   0.016704,   0.009118,  -0.088559,  -0.067503,    0.088559,     0.054882,
     0.005271,  -0.001795,   0.005844,  -0.006527,  -0.015733,   0.006542,  -0.020597,   -0.006542,    -0.019407,
     0.012217,  -0.003448,  -0.042263,   0.064138,  -0.005182,  -0.134988,  -0.079721,    0.134988,     0.101109,
     0.019461,  -0.005735,  -0.006665,   0.013395,  -0.022996,  -0.001420,  -0.013228,    0.001420,    -0.019055,
     0.000000,   0.000000,  -0.048845,   0.048845,   0.000000,  -0.138178,  -0.088842,    0.138178,     0.088817,
     0.022649,  -0.002830,   0.004947,   0.004947,  -0.019589,   0.000000,  -0.014737,   -0.000000,    -0.014741,
    -0.012217,   0.003448,  -0.064138,   0.042263,   0.005182,  -0.134988,  -0.101219,    0.134988,     0.079705,
     0.019461,  -0.005735,   0.013395,  -0.006665,  -0.022996,   0.001420,  -0.018808,   -0.001420,    -0.013217,
     0.012745,  -0.008094,  -0.085459,   0.113785,   0.002792,  -0.134270,  -0.063082,    0.134270,     0.086800,
     0.039573,  -0.011360,  -0.004205,   0.013985,  -0.020846,   0.013096,   0.000750,   -0.013096,    -0.009037,
    -0.012745,   0.008094,  -0.113785,   0.085459,  -0.002792,  -0.134270,  -0.086941,    0.134270,     0.063019,
     0.039573,  -0.011360,   0.013985,  -0.004205,  -0.020846,  -0.013096,  -0.008899,    0.013096,     0.000746,
     0.008654,  -0.002118,  -0.012365,   0.025727,  -0.019616,  -0.066554,  -0.030756,    0.066554,     0.069918,
     0.003497,  -0.001276,  -0.014071,   0.008843,  -0.032046,  -0.041551,  -0.032057,    0.041551,    -0.002250,
     0.011032,  -0.002679,  -0.019256,   0.036220,  -0.017551,  -0.098034,  -0.055586,    0.098034,     0.091227,
     0.004048,  -0.001544,  -0.015372,   0.010157,  -0.034182,  -0.037453,  -0.036060,    0.037453,    -0.005837,
     0.015469,  -0.003856,  -0.029304,   0.052425,  -0.017910,  -0.124694,  -0.070971,    0.124694,     0.114063,
     0.003343,  -0.001598,  -0.017875,   0.008649,  -0.034542,  -0.031420,  -0.033724,    0.031420,    -0.013374,
     0.014198,  -0.003726,  -0.039106,   0.064256,  -0.010398,  -0.140364,  -0.081399,    0.140364,     0.113799,
     0.010418,  -0.003566,  -0.013914,   0.010963,  -0.032107,  -0.019283,  -0.025475,    0.019283,    -0.013915,
    -0.008654,   0.002118,  -0.025727,   0.012365,   0.019616,  -0.066554,  -0.072897,    0.066554,     0.030760,
     0.003498,  -0.001276,   0.008843,  -0.014071,  -0.032046,   0.041551,  -0.005429,   -0.041551,    -0.032057,
    -0.011032,   0.002679,  -0.036220,   0.019256,   0.017551,  -0.098034,  -0.091957,    0.098034,     0.055597,
     0.004048,  -0.001544,   0.010157,  -0.015372,  -0.034182,   0.037453,  -0.006969,   -0.037453,    -0.036057,
    -0.015469,   0.003856,  -0.052425,   0.029304,   0.017910,  -0.124694,  -0.114055,    0.124694,     0.070974,
     0.003343,  -0.001598,   0.008649,  -0.017875,  -0.034542,   0.031420,  -0.012945,   -0.031420,    -0.033711,
    -0.014198,   0.003726,  -0.064256,   0.039106,   0.010398,  -0.140364,  -0.113929,    0.140364,     0.081385,
     0.010418,  -0.003566,   0.010963,  -0.013914,  -0.032107,   0.019283,  -0.013573,   -0.019283,    -0.025460;

  JacobiSVD<MatrixXd> svd( U, ComputeFullV | ComputeFullU );
  svddiag = svd.matrixU();
  svddiag.conservativeResize(svddiag.rows(),9);
  _U = svddiag;

  U.resize(0,0);
}

void TRIANGLE_MESH::qTou()
{
  _u = _U * _q;
}
// void TRIANGLE_MESH::setVelocity(const VEC2& v)
// {
//   _velocity[0] = v[0];
//   _velocity[1] = v[1];
// }

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addBodyForce(const VEC2& bodyForce)
{
  for (unsigned int x = 0; x < _fExternal.size() / 2; x++)
  {
    _fExternal[2 * x]     += bodyForce[0];
    _fExternal[2 * x + 1] += bodyForce[1];
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the shear test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepShearTest(const Real shear)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][1] > 0)
      _vertices[right][0] += shear;
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepStretchTest(const Real stretch)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.5)
      _vertices[right][0] += stretch;
  }
}

void TRIANGLE_MESH::stretch2(const Real stretch)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][1] > 0)
      _vertices[right][1] += stretch;
  }
}

///////////////////////////////////////////////////////////////////////
//this will find the global force vector of forces on each unrestrained
//vertex.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMaterialForces()
{
  //the global vector will be v= [f_0, f_1 ...] where each f_i = [x,y] (column vectors)
  //and forces are only for unconstrained vertices. This means the size of this vector is
  // size(unconstrained vertices)*2 since each has an x,y component
  VECTOR global_vector(_unconstrainedVertices.size()*2);
  global_vector.setZero();

  //loop through triangles, calculate force on each local point
  //convert local vertices to global vertex
  //add force computed into global vertex spot of global_vector
  int n_of_triangles = _triangles.size();

  for(int x = 0; x < n_of_triangles; x++)
  {
    //get our current triangle
    TRIANGLE current = getTriangle(x);

    //find the force vector that's being applied to this triangle
    VECTOR current_force = current.computeForceVector();

    //add the forces into the right place in the global force vector
    for(int y = 0; y < 3; y++)
    {
      //get current triangle vertex y
      VEC2* current_v = _triangles[x].vertex(y);

      //if the vertex we're on is unconstrained, add it to the global vector
      if(_vertexToIndex.find(current_v) != _vertexToIndex.end())
      {
        //find the global index using the vertexToIndex map
        int global_index = _vertexToIndex[current_v];
        global_vector[global_index] += current_force[2*y];
        global_vector[global_index + 1] += current_force[2*y + 1];
      }
    }
  }

  //store global vector into _f
  _f = global_vector;
}

// my attempt at writing a motion thing.... ignore for now.
void TRIANGLE_MESH::stepMotion(float dt, const VEC2& outerForce)
{
  //first find displacements of constrained vertices using velocity... figure out later
  // vector<VEC2> oldVerts(_u.size());
  // for(int i = 0; i < _u.size(); i++)
  //   oldVerts[i] = _u[i];
  //
  // // then find u
  // for(int i = 0; i < 5; i++)
  // {
  //   vector<VEC2> du(_unconstrainedVertices.size());
  //
  //   vector<VEC2> ddu(_unconstrainedVertices.size());

    // tiny time step
    // for(int i = 0; i < _unconstrainedVertices.size(); i++)
    // {
    //   _u[i][0] = _vertices[_unconstrainedVertices[i]][0] + dt*(_velocity[0]);
    //   _u[i][1] = _vertices[_unconstrainedVertices[i]][1] + dt*(_velocity[1]);
    // }
    //
    // for(int i = 0; i < _unconstrainedVertices.size(); i++)
    // {
    //   du[i][0] = (_u[i][0] - oldVerts[i][0])/dt;
    //   du[i][1] = (_u[i][1] - oldVerts[i][1])/dt;
    // }
    //
    // for(int i = 0; i < _unconstrainedVertices.size(); i++)
    // {
    //   ddu[i][0] = (du[i][0] - oldVerts[i][0])/dt;
    //   ddu[i][1] = (du[i][1] - oldVerts[i][1])/dt;
    // }

  // }
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::stepQuasistatic()
{
  // printMatrix(_U); <---- _U is the same
  //make stiffness Matrix K. size is 2*unrestrained vertices x  2*unrestrained vertices
  MATRIX K(_unconstrainedVertices.size()*2,_unconstrainedVertices.size()*2);
  MATRIX inverse;

  //step 1: compute K
  K.setZero();
  computeStiffnessMatrix(K);

  //step 2: compute internal material forces, R(uq)
  computeMaterialForces();

  // step 3: convert R(Uq) to R'
  VECTOR reducedR = _U.transpose() * _f;

  //step 4: external forces transform
  VECTOR reducedF = _U.transpose() * _fExternal;

  //step 5: form the residual (r = F + E)
  VECTOR r2 = -1*(reducedR + reducedF);

  MATRIX k_reduced = _U.transpose() * (K * _U);

  //step 6: compute x = K .inverse().eval()  * r
  MATRIX inverse2 = k_reduced.inverse().eval();
  VECTOR x2 = inverse2*r2;

  //step 6: add solution x to displacement vector _u
  _q += x2;
  qTou();

  //step 7: update all node positions w new displacement vector
  for(int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    VEC2 displacement;
    displacement[0] = _u[2*x];
    displacement[1] = _u[2*x + 1];
    _vertices[_unconstrainedVertices[x]] = _restVertices[_unconstrainedVertices[x]] + displacement;
    // printf("displacement at vertex %d is (%f, %f)\n", _unconstrainedVertices[x], displacement[0], displacement[1]);
  }

  //reset forces to 0
  _fExternal.setZero();
  _f.setZero();

  static int counter = 0;
  cout << " Quasistatic step: " << counter << endl;
  counter++;
  return true;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStiffnessMatrix(MATRIX& K)
{
  //can assume K is the correct size, 2V x 2V

  int n_of_triangles = _triangles.size();

  //go through each triangle and map vertices to their global vertices. if
  //unrestrained, then add to the poisition in the global K.
  for(int x = 0; x < n_of_triangles; x++)
  {
    TRIANGLE current = getTriangle(x);

    //find the force vector that's being applied to this triangle
    MATRIX current_force = current.computeForceJacobian();

    for(int i = 0; i < 3; i++)
    {
      //get current triangle vertex y
      VEC2* current_iv = _triangles[x].vertex(i);

      //if the vertex we're on is unconstrained, add it to the global K
      if(_vertexToIndex.find(current_iv) != _vertexToIndex.end())
      {
        //find the global index using the vertexToIndex map
        int global_index_i = _vertexToIndex[current_iv];

        for(int j = 0; j < 3; j++)
        {
          VEC2* current_jv = _triangles[x].vertex(j);

          if(_vertexToIndex.find(current_jv) != _vertexToIndex.end())
          {
            int global_index_j = _vertexToIndex[current_jv];
            K(global_index_i, global_index_j) += current_force(2*i, 2*j);
            K(global_index_i, global_index_j + 1) += current_force(2*i, 2*j + 1);
            K(global_index_i + 1, global_index_j) += current_force(2*i + 1, 2*j);
            K(global_index_i + 1, global_index_j + 1) += current_force(2*i + 1, 2*j + 1);
          }
        }
      }
    }
  }
}
