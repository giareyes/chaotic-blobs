#include "../EXTRAFUNCTIONS.h"
#include <vector>

using namespace std;

class TENSOR3
{
public:
  TENSOR3();
  TENSOR3(int rows, int cols, int slabs);
  TENSOR3(const std::vector<MATRIX>& slabs);

  int rows() const { return _rows; };
  int cols() const { return _cols; };
  int slabs() const { return _slabs; };


  MATRIX modeThreeProduct(const VECTOR& x);
  TENSOR3 modeThreeProduct(const MATRIX& x);
  TENSOR3 multiply(double x);

  void toString();

  vector<MATRIX> _tensor;
protected:
  int _rows;
  int _cols;
  int _slabs;
};
