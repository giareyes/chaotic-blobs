#include "TENSOR3.h"


TENSOR3::TENSOR3(int rows, int cols, int slabs):
  _rows(rows), _cols(cols), _slabs(slabs)
{
  for (int x = 0; x < slabs; x++)
  {
    MATRIX m(rows, cols);
    m.setZero();
    _tensor.push_back(m);
  }
}

TENSOR3::TENSOR3(const vector<MATRIX>& slabs)
{
  assert(slabs.size() > 0);

  _rows = slabs[0].rows();
  _cols = slabs[0].cols();
  _slabs = slabs.size();

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

TENSOR3 TENSOR3::modeThreeProduct(const MATRIX& x)
{
  // need to look at cols x slabs matrices ???
  assert(x.rows() == _slabs || x.cols() == _cols);

  // i guess we can create a new tensor 3, and put all of the info inside through loops to reorder rows/cols/slabs
  TENSOR3 reordered(_cols, _slabs, _rows);
  
  vector<MATRIX> result;

  if(x.rows() == _slabs)
  {
    for(int i = 0; i < _rows; i++)
    {
      MATRIX mid = x * reordered[i];
      result.push_back(mid);
    }
  }
  else
  {
    for(int i = 0; i < _rows; i++)
    {
      MATRIX mid = reordered[i] * x;
      result.push_back(mid);
    }
  }

  TENSOR3 result_tensor(result);
  return result_tensor;
}
