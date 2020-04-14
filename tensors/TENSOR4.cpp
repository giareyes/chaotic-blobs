#include "TENSOR4.h"
#include "../EXTRAFUNCTIONS.h"
#include <stdlib.h>


TENSOR4::TENSOR4(int rows, int cols, int slab_rows, int slab_cols)
{
  _rows = rows;
  _cols = cols;
  _slab_rows = slab_rows;
  _slab_cols = slab_cols;

  for (int x = 0; x < slab_cols; x++)
  {
    TENSOR3 m(rows, cols, slab_rows);
    _tensor.push_back(m);
  }
}

TENSOR4::TENSOR4(const std::vector<TENSOR3>& slabs)
{
  assert(slabs.size() > 0);

  _rows = slabs[0].rows();
  _cols = slabs[0].cols();
  _slab_rows = slabs[0].slabs();
  _slab_cols = slabs.size();

  for (int x = 0; x < _slab_cols; x++)
  {
    assert(slabs[x].rows() == _rows);
    assert(slabs[x].cols() == _cols);
    assert(slabs[x].slabs() == _slab_rows);
    _tensor.push_back(slabs[x]);
  }
}

TENSOR3 TENSOR4::modeFourProduct(const VECTOR& x)
{
  assert(_slab_cols == x.size());
  assert(_slab_cols > 0);
  TENSOR3 result(_rows, _cols, _slab_rows);

  for(int i = 0; i < _slab_rows; i++)
  {
    for(int j = 0; j < _slab_cols; j++)
    {
      result._tensor[i] = result._tensor[i] + x[j]*_tensor[j]._tensor[i];
    }
  }

  return result;
}

void TENSOR4::toString()
{
  for(int k = 0; k < _slab_cols; k++)
  {
    printf("column %d: { \n", k);
    _tensor[k].toString();
    printf("}\n");
  }
}

TENSOR4 TENSOR4::modeFourProduct(const MATRIX& x)
{
  // must make sure that the columns of the matrix we are multiplying match the
  // number of slab columns in our tensor
  assert( x.cols() == _slab_cols);

  vector<TENSOR3> result;
  for(int p = 0; p < x.rows(); p++)
  {
    TENSOR3 temp(_rows, _cols, _slab_rows);
    result.push_back(temp);
  }

  // loop through all of the slabs in the tensor, multiply by every constant in the
  // matrix x, and add together.
  // Equation: result(a,b,c,d) = sum[i_n = 0 to i_n = x.cols()] (tensor(a,b,c,i_n)*x(d, i_n))
  for(int k = 0; k < _slab_rows; k++)
  {
    for(int i = 0; i < x.rows(); i++)
    {
      for(int j = 0; j < x.cols(); j++)
      {
        MATRIX mid = _tensor[j]._tensor[k] * x(i,j);
        result[i]._tensor[k] += mid;
      }
    }
  }

  TENSOR4 result_tensor(result);

  // free result ?
  // for (int i = 0; i < _rows; i++)
  //   result[i].resize(0,0);

  return result_tensor;
}

TENSOR4 TENSOR4::modeThreeProduct(const MATRIX& x)
{
  // need to make sure that the number of columns in the matrix x match
  // the number of rows in the tensor
  assert( x.cols() == _slab_rows);

  // initialize the vector of third order tensors that will become our
  // resulting fourth order tensor. We must initialize because we are
  // altering rows, not columns, so it works better this way.
  vector<TENSOR3> result;
  for(int p = 0; p < _slab_cols; p++)
  {
    TENSOR3 temp(_rows, _cols, x.rows());
    result.push_back(temp);
  }

  // loop through all of the slabs in the tensor, multiply by every constant in the
  // matrix x, and add together.
  // Equation: result(a,b,c,d) = sum[i_n = 0 to i_n = x.cols()] (tensor(a,b,i_n,d)*x(c, i_n))
  for(int k = 0; k < _slab_cols; k++)
  {
    for(int i = 0; i < x.rows(); i++)
    {
      for(int j = 0; j < x.cols(); j++)
      {
        MATRIX mid = _tensor[k]._tensor[j] * x(i,j);
        result[k]._tensor[i] += mid;
      }
    }
  }

  TENSOR4 result_tensor(result);

  // free result ?
  // for (int i = 0; i < _rows; i++)
  //   result[i].resize(0,0);

  return result_tensor;
}

TENSOR4& TENSOR4::operator+=(const TENSOR4& m)
{
  assert(_rows == m.rows());
  assert(_cols == m.cols());
  assert(_slab_rows == m.slab_rows());
  assert(_slab_cols == m.slab_cols());

  for (int z = 0; z < _slab_cols; z++)
  {
      _tensor[z] += m._tensor[z];
  }

  return *this;
}

TENSOR4 TENSOR4::modeTwoProduct(const MATRIX& x)
{
  assert( x.cols() == _cols);

  vector<TENSOR3> result;
  for(int p = 0; p < _slab_cols; p++)
  {
    TENSOR3 temp(_rows, x.rows(), _slab_rows);
    result.push_back(temp);
  }

  // Hi, I'm really sorry about the code below, it's shameful
  for(int k = 0; k < _slab_cols; k++)
  {
    for(int l = 0; l < _slab_rows; l++)
    {
      for(int i = 0; i < x.rows(); i++)
      {
        for(int m = 0; m < _rows; m++)
        {
          for(int j = 0; j < x.cols(); j++)
          {
            result[k]._tensor[l](m,i) += _tensor[k]._tensor[l](m,j) * x(i,j);
          }
        }
      }
    }
  }

  TENSOR4 result_tensor(result);

  // free result ?
  // for (int i = 0; i < _rows; i++)
  //   result[i].resize(0,0);

  return result_tensor;
}

TENSOR4 TENSOR4::modeOneProduct(const MATRIX& x)
{
  assert( x.cols() == _rows);

  vector<TENSOR3> result;
  for(int p = 0; p < _slab_cols; p++)
  {
    TENSOR3 temp(x.rows(), _cols, _slab_rows);
    result.push_back(temp);
  }

  // Hi, I'm really sorry about the code below, it's shameful
  for(int k = 0; k < _slab_cols; k++)
  {
    for(int l = 0; l < _slab_rows; l++)
    {
      for(int i = 0; i < x.rows(); i++)
      {
        for(int m = 0; m < _cols; m++)
        {
          for(int j = 0; j < x.cols(); j++)
          {
            result[k]._tensor[l](i,m) += _tensor[k]._tensor[l](j,m) * x(i,j);
          }
        }
      }
    }
  }

  TENSOR4 result_tensor(result);

  // free result ?
  // for (int i = 0; i < _rows; i++)
  //   result[i].resize(0,0);

  return result_tensor;
}

void printMatrix(MATRIX matrix)
{
  int rows = matrix.rows();
  int cols = matrix.cols();
  for(int i = 0; i < rows; i++)
  {
    printf("row %d: ( ", i+1);
    for(int j = 0; j < cols; j++)
    {
      printf("%f, ", matrix(i,j));
    }
    printf(")\n");
  }
}

int main(int argc, char** argv)
{
  vector<TENSOR3> vec_tens;

  MATRIX x(6,4);
  x.setZero();
  VECTOR vec(6);

  vec[0] = 5;
  vec[1] = 1;
  vec[2] = 4;
  vec[3] = 3;
  vec[4] = 2;
  vec[5] = 1;

  // make 4 slab columns
  for(int i = 0; i < 4; i++)
  {
    vector<MATRIX> matrix_vector;
    // make 4 slab rows
    for(int j = 0; j < 4; j++ )
    {
      MATRIX m(4,4);
      m.setZero();
      // for 4 rows and for cols
      for(int k = 0; k < 4; k ++)
      {
        for(int l = 0; l < 4; l++)
        {
          m(k,l) = rand() % 5;
        }
      }
      matrix_vector.push_back(m);
    }
    TENSOR3 temp(matrix_vector);
    vec_tens.push_back(temp);
  }

  for(int i = 0; i < 24; i++)
  {
    int rand_x = rand() % 6;
    int rand_y = rand() % 4;

    x(rand_x, rand_y) = rand() % 5;
  }

  TENSOR4 fourTensor(vec_tens);

  // result1 = C * df/dx * df/dx * x * x
  TENSOR4 fourth_res = fourTensor.modeFourProduct(x);
  fourth_res = fourth_res.modeThreeProduct(x);
  fourth_res = fourth_res.modeTwoProduct(x);
  fourth_res = fourth_res.modeOneProduct(x);
  TENSOR3 t3_1 = fourth_res.modeFourProduct(vec);
  MATRIX result1 = t3_1.modeThreeProduct(vec);

  // result1 = x * result1;
  // result1 = result1 * x.transpose();

  // result2 = C * [df/dx T  * x] * [df/dx T  * x]
  VECTOR result = x.transpose() * vec;
  TENSOR3 t3_2 = fourTensor.modeFourProduct(result);
  MATRIX result2 = t3_2.modeThreeProduct(result);

  result2 = x * result2;
  result2 = result2 * x.transpose();

  printf("\n\n\n------------------------- result 1 ----------------------------\n\n\n");
  // fourth_res.toString();
  printMatrix(result1);
  printf("\n\n\n------------------------- result 2 ----------------------------\n\n\n");
  printMatrix(result2);

  return 1;
}
