#ifndef ANN_MATRIX_H
#define ANN_MATRIX_H

#include <iostream>
#include <random>


/// namespace of this artificial neural networks course
namespace ann
{

/// Class ann::matrix implements the concept of a matrix.
/// Basic matrix operations are implemented, as well as pseudoinversion and
/// initialization with random values.
class matrix
{
  // - - - type definitions - - - - - - - - - - - - - - - - - - - - - - - - - - -

 public:
  /// Defines which (floating point) type matrix elements shall be of.
  typedef long double real;

 private:
  /// This container type is used to internally store matrix elements.
  typedef std::vector<real> data;

 public:
  /// The type used for array access (ideally) fits the utilized container.
  typedef data::size_type size;
  /// The signed type used for array access fits the utilized container.
  typedef data::difference_type diff;

 private:
  /// The type of random engine used by matrix is defined here for convenience.
  /// This engine implemnts a 64-bit Mersenne Twister as suggested by Matsumoto
  /// and Nishimura in 2000.
  typedef std::mt19937_64 random_engine;

  // - - - member variables - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /// number of rows
  size rows_m;
  /// number of columns
  size cols_m;

  /// data array holding matrix elements
  data data_m;

  /// This random device is used to seed the random engine.
  /// It is a non-deterministic uniform random number generator, if one is
  /// available on your system.
  /// All matrices in a thread share the same random device.
  static thread_local std::random_device random_device_m;
  /// Each matrix object has its own (instance of the) random engine.
  random_engine random_engine_m;

 public:
  /// primary constructor
  /// Creates a matrix with the given number of columns and rows.
  /// All elements are initialized to the third parameter, which defaults to 0.
  /// This is also the default constructor, creating a matrix of size 1x1.
  /// When giving only the first parameter, a column vector is created.
  matrix (size rows_a = 1, size cols_a = 1, real value_a = 0.0);

  /// copy constructor
  /// Constructs a matrix which is a copy of the given other matrix.
  matrix (const matrix& other_a);

  /// move constructor
  /// Constructs a matrix taking over resources of the given other matrix.
  matrix (matrix&& other_a);

  // There's no need to explicitly define a destructor. The implicitly defined
  // default is fine, because all member variables may be destroyed in a default
  // way (e.g. no explicit deallocation on the heap is needed).

  /// assignment operator
  /// Turns this matrix into a copies of the given other matrix.
  matrix& operator= (const matrix& other_a);

  /// move assignement operator
  /// Lets this matrix take over the resources of the given other matrix.
  matrix& operator= (matrix&& other_a);

  /// assignement operator
  /// Assigns the given scalar to each element of the matrix.
  matrix& operator= (real value_a);

  /// Returnes the number of rows in this matrix.
  size rows () const;
  /// Returnes the number of columns in this matrix.
  size cols () const;

  /// Gets element at the given row and column.
  /// Warning: Parameters are considered indices starting at 0.
  real get (size row_a, size col_a) const;

  /// Sets element at the given row and column.
  /// Warning: Parameters are considered indices starting at 0.
  void set (size row_a, size col_a, real value_a);

  /// Accesses element (read/write) in the given row and column of the matrix.
  /// This version of the operator returns a non-const reference to allow
  /// modification of the matrix element.
  /// Warning: Parameters are considered indices starting at 0.
  /// E.g.: To access the element in the first row and second column use (0, 1).
  real& operator() (size row_a, size col_a);
  /// Accesses element (read only) in the given row and column of the matrix.
  /// This version of the operator returns a const reference to allow access to
  /// elements in case this matrix is const.
  /// Warning: Parameters are considered indices starting at 0.
  /// E.g.: To access the element in the first row and second column use (0, 1).
  const real& operator() (size row_a, size col_a) const;

  /// Accesses element (read/write) at the given index.
  /// This is useful for vectors (matrices with only one row or column),
  /// allowing access without having to give the superfluous row or column
  /// index as argument.
  real& operator() (size index_a);

  /// Accesses element (read only) at the given index.
  /// This is useful for vectors (matrices with only one row or column),
  /// allowing access without having to give the superfluous row or column
  /// index as argument.
  const real& operator() (size index_a) const;

  /// scalar cast operator
  /// Access the first element as if the matrix was a scalar.
  /// Returns a reference to element (0,0) when using this matrix in a context
  /// where a scalar is required.
  operator real& ();
  /// scalar cast operator
  /// Access the first element as if the matrix was a scalar.
  /// Returns a const reference to element (0,0) when using this matrix in a
  /// context where a scalar is required.
  operator const real& () const;

  /// Extracts a column from this matrix.
  /// Returns a matrix with the same number of rows as this matrix and one
  /// column, where elements are copied from the column at the given index
  /// of this matrix.
  matrix column (size col_a) const;

  /// Extracts a row from this matrix.
  /// Returns a matrix with the same number of columns as this matrix and one
  /// row, where elements are copied from the row at the given index
  /// of this matrix.
  matrix row (size row_a) const;

  /// matrix addition
  /// Creates a new matrix as the sum of this and the given other matrix.
  matrix operator+ (const matrix& other_a) const;
  /// matrix addition, compound assignment
  /// Adds the given matrix to this matrix, storing the result in this.
  matrix& operator+= (const matrix& other_a);

  /// matrix subtraction
  /// Creates a new matrix as the difference of this and the given other matrix.
  matrix operator- (const matrix& other_a) const;
  /// matrix subtraction, compound assignment
  /// Subtracts the given matrix from this matrix, storing the result in this.
  matrix& operator-= (const matrix& other_a);

  /// matrix multiplication
  /// Creates a new matrix as the product of this and the given other matrix.
  matrix operator* (const matrix& other_a) const;
  /// matrix multiplication, compound assignment
  /// Multiplies the given matrix with this matrix, storing the result in this.
  /// Warning: This may change the extents of this matrix.
  matrix& operator*= (const matrix& other_a);

  /// scalar multiplication
  /// Creates a copy of this matrix where all elements are scaled by the given
  /// factor.
  matrix operator* (real b_a) const;
  /// scalar multiplication, compound assignment
  /// Scales this matrix by the given factor.
  matrix& operator*= (real b_a);

  /// scalar division
  /// Creates a copy of this matrix where all elements are divided by the given
  /// value.
  matrix operator/ (real b_a) const;
  /// scalar division, compound assignment
  /// Divides this matrix by the given value.
  matrix& operator/= (real b_a);

  /// transpose
  /// Creates a transposed copy of this matrix.
  matrix transpose() const;

  /// pseudoinverse
  /// Computes the pseudoinverse of this matrix using singular value decomposition.
  matrix inverse() const;

  /// Fills the matrix with pseudo-random values drawn from a uniform
  /// distribution over the interval defined by the given lower and upper bound,
  /// where the lower bound is inclusive and the upper is not.
  matrix& fill_with_uniform_samples (real lower_bound_a, real upper_bound_a);

  /// Calls fill_with_uniform_samples.
  matrix& randomize (real lower_bound_a = -1.0, real upper_bound_a = 1.0);

  /// Compares two matrices for equality.
  bool operator== (const matrix& other_a) const;
  /// Compares two matrices for inequality.
  bool operator!= (const matrix& other_a) const;

  /// Returns true, iff this is a Hermitian (or self-adjoint) matrix.
  bool is_hermitian () const;

 private:
  /// Computes the singular value decomposition of this matrix.
  /// This matrix is decomposed into U_a * W_a * Vc_a, where U and V* are unitary
  /// and W is a diagonal matrix.
  void singular_value_decomposition(matrix& U_a, matrix& Vc_a, data& W_a) const;

}; // class matrix

} // namespace ann

/// scalar multiplication (with scalar on the left hand side).
/// Creates a copy of the matrix where all elements are scaled by the given
/// factor.
ann::matrix operator* (ann::matrix::real factor_a, ann::matrix matrix_a);

/// Overload of operator<< for writing matrices to output streams.
/// Must be declared outside of class matrix to achieve writing to stream in the
/// convenient fashion (stream << stuff).
std::ostream& operator<< (std::ostream& stream_a, const ann::matrix& matrix_a);

#endif //ANN_MATRIX_H
