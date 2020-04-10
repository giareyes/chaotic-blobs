#include "TENSOR3.h"
#include "../EXTRAFUNCTIONS.h"


TENSOR3::TENSOR3(int rows, int cols, int slabs)
{
  _rows = rows;
  _cols = cols;
  _slabs = slabs;

  // _tensor = new MATRIX[slabs];

  for (int x = 0; x < slabs; x++)
  {
    MATRIX m(rows, cols);
    m.setZero();
    _tensor.push_back(m);
  }
}

TENSOR3::TENSOR3(const std::vector<MATRIX>& slabs)
{
  assert(slabs.size() > 0);

  _rows = slabs[0].rows();
  _cols = slabs[0].cols();
  _slabs = slabs.size();

  // _tensor = new MATRIX[slabs.size()];

  for (int x = 0; x < _slabs; x++)
  {
    assert(slabs[x].rows() == _rows);
    assert(slabs[x].cols() == _cols);
    _tensor.push_back(slabs[x]);
  }
}

MATRIX TENSOR3::modeThreeProduct(const VECTOR& x)
{
  assert(_slabs == x.size());
  assert(_slabs > 0);
  MATRIX result(_rows, _cols);
  result.setZero();

  for(int i = 0; i < _slabs; i++)
  {
    result = result + x[i]*_tensor[i];
  }

  return result;
}


void TENSOR3::toString()
{
  for(int k = 0; k < _slabs; k++)
  {
    printf("slab %d: [ \n", k);
    for(int i = 0; i < _rows; i++)
    {
      printf("row %d: ( ", i);
      for(int j = 0; j < _cols; j++)
      {
        printf("%f, ", _tensor[k](i,j));
      }
      printf(")\n");
    }
    printf("]\n");
  }
}

TENSOR3 TENSOR3::modeThreeProduct(const MATRIX& x)
{
  // need to look at cols x slabs matrices ???
  assert( x.cols() == _slabs);

  vector<MATRIX> result;

  for(int i = 0; i < x.rows(); i++)
  {
    MATRIX m(_rows, _cols);
    m.setZero();
    for(int j = 0; j < x.cols(); j++)
    {
      MATRIX mid = _tensor[j] * x(i,j);
      m += mid;
    }
    result.push_back(m);
  }

  TENSOR3 result_tensor(result);

  // free result ?
  for (int i = 0; i < _rows; i++)
    result[i].resize(0,0);

  return result_tensor;
}

// int main(int argc, char** argv)
// {
//   vector<MATRIX> matrix_vector;
//   MATRIX x(6,4);
//   x.setZero();
//   // VECTOR vec(4);
//
//   for(int i = 0; i < 4; i++)
//   {
//     MATRIX m(4,4);
//     m.setIdentity();
//     matrix_vector.push_back(m);
//
//     x(i,i) = 2;
//     // vec[i] = 3;
//   }
//
//   TENSOR3 newTensor(matrix_vector);
//   TENSOR3 result = newTensor.modeThreeProduct(x);
//   result.toString();
//   return 1;
// }