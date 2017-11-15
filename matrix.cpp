#include "matrix.h"

#include "floating.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>


/// Overload of operator<< for writing matrices to output streams.
/// Must be declared outside of class matrix to achieve writing to stream in the
/// convenient fashion (stream << stuff).
std::ostream&
operator<<
(
  std::ostream& stream_a,
  const ann::matrix& matrix_a
)
{
  // Draw row-wise, as is convenient for text output.
  for (ann::matrix::size row_l = 0; row_l < matrix_a.rows(); ++row_l)
  {
    stream_a << matrix_a(row_l, 0);
    for (ann::matrix::size col_l = 1; col_l < matrix_a.cols(); ++col_l)
    {
      stream_a << "\t" << matrix_a(row_l, col_l);
    }
    stream_a << "\n";
  }
  // Return the new position in the stream, as usual for overloads of operator<<.
  return stream_a;
}

/// scalar multiplication (with scalar on the left hand side).
/// Creates a copy of the matrix where all elements are scaled by the given
/// factor.
ann::matrix
operator*
(
  ann::matrix::real factor_a,
  ann::matrix matrix_a
)
{
  matrix_a *= factor_a;

  return matrix_a;
}


// Definition of class matrix - - - - - - - - - - - - - - - - - - - - - - - - - -
thread_local std::random_device ann::matrix::random_device_m;

/// primary constructor
/// Creates a matrix with the given number of columns and rows.
/// All elements are initialized to the third parameter, which defaults to 0.
/// This is also the default constructor, creating a matrix of size 0.
ann::matrix
::
matrix
(
  size rows_a,
  size cols_a,
  real value_a
)
: rows_m(rows_a),
  cols_m(cols_a),
  data_m(rows_m * cols_m, value_a),
  random_engine_m(random_device_m())
{
}

/// copy constructor
/// Constructs a matrix which is a copy of the given other matrix.
ann::matrix
::
matrix
(
  const matrix& other_a
)
: rows_m(other_a.rows_m),
  cols_m(other_a.cols_m),
  data_m(other_a.data_m),
  random_engine_m(random_device_m())
{
}

/// move constructor
/// Constructs a matrix taking over resources of the given other matrix.
ann::matrix
::
matrix
(
 matrix&& other_a
)
: rows_m(other_a.rows_m),
  cols_m(other_a.cols_m),
  data_m(std::move(other_a.data_m)),
  random_engine_m(random_device_m())
{
}


/// assignment operator
/// Turns this matrix into a copies of the given other matrix.
ann::matrix&
ann::matrix
::
operator=
(
 const matrix& other_a
)
{
  rows_m = other_a.rows_m;
  cols_m = other_a.cols_m;
  data_m = other_a.data_m;

  return *this;
}

/// move assignement operator
/// Lets this matrix take over the resources of the given other matrix.
ann::matrix&
ann::matrix
::
operator=
(
  matrix&& other_a
)
{
  rows_m = other_a.rows_m;
  cols_m = other_a.cols_m;
  data_m = std::move(other_a.data_m);

  return *this;
}

/// assignement operator
/// Assigns the given scalar to each element of the matrix.
ann::matrix&
ann::matrix
::
operator=
(
  const real value_a
)
{
  for (real& element_l : data_m)  element_l = value_a;

  return *this;
}


/// Returnes the number of rows in this matrix.
ann::matrix::size
ann::matrix
::
rows
() const
{
  return rows_m;
}

/// Returnes the number of columns in this matrix.
ann::matrix::size
ann::matrix
::
cols
() const
{
  return cols_m;
}


/// Gets element at the given row and column.
ann::matrix::real
ann::matrix
::
get
(
  size row_a,
  size col_a
) const
{
  return operator()(row_a, col_a);
}

/// Sets element at the given row and column.
void
ann::matrix
::
set
(
  size row_a,
  size col_a,
  real value_a
)
{
  operator()(row_a, col_a) = value_a;
}

/// Accesses element (read/write) in the given row and column of the matrix.
/// This version of the operator returns a non-const reference to allow
/// modification of the matrix element.
/// Warning: Parameters are considered indices starting at 0.
/// E.g.: To access the element in the first row and second column use (0, 1).
ann::matrix::real&
ann::matrix
::
operator()
(
  const size row_a,
  const size col_a
)
{
  // Assert that indices are in bounds.
  assert(row_a < rows() && col_a < cols());
  // Return the specified element, retreived by 1D index.
  return data_m[row_a + rows() * col_a];
}

/// Accesses element (read only) in the given row and column of the matrix.
/// This version of the operator returns a const reference to allow access to
/// elements in case this matrix is const.
/// Warning: Parameters are considered indices starting at 0.
/// E.g.: To access the element in the first row and second column use (0, 1).
const ann::matrix::real&
ann::matrix
::
operator()
(
  const size row_a,
  const size col_a
) const
{
  // Assert that indices are in bounds.
  assert(row_a < rows() && col_a < cols());
  // Return the specified element, retreived by 1D index.
  return data_m[row_a + rows() * col_a];
}


/// Accesses element (read/write) at the given index.
/// This is useful for vectors (matrices with only one row or column),
/// allowing access without having to give the superfluous row or column
/// index as argument.
ann::matrix::real&
ann::matrix
::
operator()
(
  size index_a
)
{
  return data_m.at(index_a);
}

/// Accesses element (read only) at the given index.
/// This is useful for vectors (matrices with only one row or column),
/// allowing access without having to give the superfluous row or column
/// index as argument.
const ann::matrix::real&
ann::matrix
::
operator()
(
  size index_a
) const
{
  return data_m.at(index_a);
}


/// scalar cast operator
/// Access the first element as if the matrix was a scalar.
/// Returns a reference to element (0,0) when using this matrix in a context
/// where a scalar is required.
ann::matrix
::
operator ann::matrix::real&
()
{
  return data_m.front();
}

/// scalar cast operator
/// Access the first element as if the matrix was a scalar.
/// Returns a const reference to element (0,0) when using this matrix in a
/// context where a scalar is required.
ann::matrix
::
operator const ann::matrix::real&
() const
{
  return data_m.front();
}


/// Extracts a column from this matrix.
/// Returns a matrix with the same number of rows as this matrix and one
/// column, where the elements are copied from the column at the given index
/// of this matrix.
ann::matrix
ann::matrix
::
column
(
  size col_a
) const
{
  matrix vector_l(rows(), 1);
  for (size row_l = 0; row_l < rows(); ++row_l)
  {
    vector_l(row_l) = get(row_l, col_a);
  }

  return vector_l;
}

/// Extracts a row from this matrix.
/// Returns a matrix with the same number of columns as this matrix and one
/// row, where the elements are copied from the row at the given index
/// of this matrix.
ann::matrix
ann::matrix
::
row
(
  size row_a
) const
{
  matrix vector_l(1, cols());
  for (size col_l = 0; col_l < cols(); ++col_l)
  {
    vector_l(col_l) = get(row_a, col_l);
  }

  return vector_l;
}



/// matrix addition
/// Creates a new matrix as the sum of this and the given other matrix.
ann::matrix
ann::matrix
::
operator+
(
  const matrix& other_a
) const
{
  // Create resulting matrix as copy of this.
  matrix result_l(*this);
  // Add the other matrix using addition assignment operator, which also asserts
  // matching matrix extents.
  result_l += other_a;

  return result_l;
}

/// matrix addition compound assignment
/// Adds the given matrix to this matrix, storing the result in this.
ann::matrix&
ann::matrix
::
operator+=
(
  const matrix& other_a
)
{
  // Matrix addition is only possible for matrices of matching extents.
  assert(rows() == other_a.rows() && cols() == other_a.cols());

  // Add matrix elements utilizing transform from header algorithm.
  std::transform
  ( data_m.begin(),
    data_m.end(),
    other_a.data_m.begin(),
    data_m.begin(),
    std::plus<real>()
  );

  // Return a reference to this to allow chaining of operations.
  return *this;
}


/// matrix subtraction
/// Creates a new matrix as the difference of this and the given other matrix.
ann::matrix
ann::matrix
::
operator-
(
  const matrix& other_a
) const
{
  // Create resulting matrix as copy of this.
  matrix result_l(*this);
  // Add the other matrix using addition assignment operator, which also asserts
  // matching matrix extents.
  result_l -= other_a;

  return result_l;
}

/// matrix subtraction, compound assignment
/// Subtracts the given matrix from this matrix, storing the result in this.
ann::matrix&
ann::matrix
::
operator-=
(
  const matrix& other_a
)
{
  // Matrix subtraction is only possible for matrices of matching extents.
  assert(rows() == other_a.rows() && cols() == other_a.cols());

  // Subtract matrix elements utilizing transform from header algorithm.
  std::transform
  ( data_m.begin(),
    data_m.end(),
    other_a.data_m.begin(),
    data_m.begin(),
    std::minus<real>()
  );

  // Return a reference to this to allow chaining of operations.
  return *this;
}


/// matrix multiplication
/// Creates a new matrix as the product of this and the given other matrix.
ann::matrix
ann::matrix
::
operator*
(
  const matrix& other_a
) const
{
  // Matrix multiplication is only possible, if extents match crosswise.
  assert(cols() == other_a.rows());

  // Create matrix of proper size to hold the result.
  matrix result_l(rows(), other_a.cols(), 0.0);

  // Fill the resulting matrix.
  for (size col_l = 0; col_l < result_l.cols(); ++col_l)
  {
    for (size row_l = 0; row_l < result_l.rows(); ++row_l)
    {
      // Traverse the "other" dimension and sum up the products.
      for (size cow_l = 0; cow_l < cols(); ++ cow_l)
      {
        result_l(row_l, col_l) += (*this)(row_l, cow_l) * other_a(cow_l, col_l);
      }
    }
  }

  // Return the resulting matrix.
  return result_l;
}

/// matrix multiplication, compound assignment
/// Multiplies the given matrix with this matrix, storing the result in this.
/// Warning: This may change the extents of this matrix.
ann::matrix&
ann::matrix
::
operator*=
(
 const matrix& other_a
)
{
  // Mutliply matrices using multiplication operator, which also asserts fitting
  // extents of the arguments, and assigning result to this.
  *this = *this + other_a;

  // Return a reference to this to allow chaining of operations.
  return *this;
}


/// scalar multiplication
/// Creates a copy of this matrix where all elements are scaled by the given
/// factor.
ann::matrix
ann::matrix
::
operator*
(
  const real b_a
) const
{
  // Create a copy of this matrix.
  matrix result_l(*this);
  // Use compound assignment to execute the computation.
  result_l *= b_a;

  return result_l;
}


/// scalar multiplication, compound assignment
/// Scales this matrix by the given factor.
ann::matrix&
ann::matrix
::
operator*=
(
  const real b_a
)
{
  // There's no need to traverse the matrix by rows and columns for this
  // operation. Traversing the container data_m directly allows the usage of a
  // nice type of loop added with C++11.
  // Note the '&' declaring this element_l as a reference, thus allowing
  // modification of the elements.
  for (real& element_l : data_m)  element_l *= b_a;
  // For all loops ommiting the curly braces is allowed, if the body of the loop
  // contains not more than one line of code. This is considered bad style by
  // some and it's best to be careful with it.


  // Return a reference to this to allow chaining of operations.
  return *this;
}


/// scalar division
/// Creates a copy of this matrix where all elements are divided by the given
/// value.
ann::matrix
ann::matrix
::
operator/
(
  real b_a
) const
{
  return operator*(1.0 / b_a);
}

/// scalar division, compound assignment
/// Divides this matrix by the given value.
ann::matrix&
ann::matrix
::
operator/=
(
  real b_a
)
{
  return operator*=(1.0 / b_a);
}


/// transpose
/// Creates a transposed copy of this matrix.
ann::matrix
ann::matrix
::
transpose () const
{
  matrix transpose_l(cols(), rows());
  for (size row_l = 0; row_l < rows(); ++row_l)
    for (size col_l = 0; col_l < cols(); ++col_l)
      transpose_l(col_l, row_l) = operator()(row_l, col_l);

  return transpose_l;
}


/// pseudoinversion
/// Inverts this matrix using singular value decomposition.
/// All eigenvalues < epsilon_a are set to 0.
ann::matrix
ann::matrix
::
inverse () const
{
  // This matrix is called A in this function.
  // Copy A to U, as expected by function singular_value_decomposition.
  matrix U_l(*this);
  // Note: Vc_l here stands for V's complex conjugate V* (which for real V is
  // its transpose). U and V are unitary matrices, so their complex cojugate is
  // also their inverse matrix.
  matrix Vc_l(cols(), cols());
  // W (also known as Sigma) is a diagonal matrix, so space can be saved by only
  // storing the elements on the main diagonal (all others are 0).
  data W_l(cols());

  // Fill U, V* and W with the singular value decomposition matrices of this.
  singular_value_decomposition(U_l, Vc_l, W_l);
  // So now: A == U W V*

  // Next, matrices U, V* and W must be inverted, so A+ (pseudoinverse of A) can
  // be assembled. First, though, an epsilon must be chosen to crop singular
  // values.

  // According to Wikipedia article on Moore-Penrose pseudoinverse, this
  // tolerance is used by MATLAB, GNU Octave and NumPy when setting small
  // elements of W_l to 0.
  // It appears to be a reasonable choice for an epsilon of a magnitude fitting
  // the context, intuitively. I have not thought through this in full depth,
  // though.
  // In general, choosing a proper epsilon (whenever an epsilon is needed) can
  // be far from trivial and should not be taken lightly. A poor choice, not
  // fitting the context, may lead to arbitrarily wrong results.
  const real epsilon_l
    = std::numeric_limits<real>::epsilon()
    * std::max(rows(), cols())
    * *std::max_element(W_l.begin(), W_l.end());
  //const real epsilon_l = 1e-12;

  // Inversion of a diagonal matrix can be achieved by transposing it and
  // inverting all non-zero entries.
  // Also, delete all singular values (i.e. the eigenvalues less than epsilon_l).
  for (real& wi_l : W_l)
  {
    if (wi_l < epsilon_l) // Note: All wi_l are non-negative, no abs needed.
      wi_l = 0.0;
    else
      wi_l = 1.0 / wi_l;
  }
  // So W_l now holds W+, the pseudoinverse of W.

  // Create the matrix soon to hold A+, the pseudoinverse of this matrix.
  matrix inverse_l(cols(), rows(), 0.0);

  // Now fill the pseudoinverse, i.e.: Compute  A+ = V W+ U*
  // V and U* (that is the transpose of U_l and Vc_l) are implictly computed by
  // switching indices when accessing the matrix variables.
  // The need of only three loops is due to the fact that W+ is a diagonal matrix.
  for (size col_l = 0; col_l < inverse_l.cols(); ++col_l)
    for (size row_l = 0; row_l < inverse_l.rows(); ++row_l)
      for (size cow_l = 0; cow_l < W_l.size(); ++cow_l)
      {
        inverse_l(row_l, col_l)
          += Vc_l(row_l, cow_l) * W_l[cow_l] * U_l(col_l, cow_l);
      }

  // Return the finally computed pseudoinverse of this matrix.
  return inverse_l;
}


/// Fills this matrix with pseudo-random values drawn from a uniform
/// distribution over the interval defined by the given lower and upper bound,
/// where the lower bound is inclusive and the upper is not.
ann::matrix&
ann::matrix
::
fill_with_uniform_samples
(
  real lower_bound_a,
  real upper_bound_a
)
{
  // Correct the arguments, if the user failed to pass them in the right order.
  if (lower_bound_a > upper_bound_a)  std::swap(lower_bound_a, upper_bound_a);

  // Create a distribution object for the utilized floating point type.
  std::uniform_real_distribution<real>
    distribution_l(lower_bound_a, upper_bound_a);

  // Use the random engine of this matrix to draw samples.
  for (real& element_l : data_m)  element_l = distribution_l(random_engine_m);

  return *this;
}

/// Calls fill_with_uniform_samples.
ann::matrix&
ann::matrix
::
randomize
(
  real lower_bound_a,
  real upper_bound_a
)
{
  return fill_with_uniform_samples(lower_bound_a, upper_bound_a);
}


/// Compares two matrices for equality.
bool
ann::matrix
::
operator==
(
  const matrix& other_a
) const
{
  // Start with a cheap test: Matrices of different sizes can't be equal.
  if (cols() != other_a.cols() || rows() != other_a.rows())  return false;

  // Now commpare the individual elements.
  auto it1_l = data_m.cbegin();
  auto it2_l = other_a.data_m.cbegin();
  auto end1_l = data_m.cend();
  for (; it1_l != end1_l; ++it1_l, ++it2_l)
  {
    if (!ann::equal(*it1_l, *it2_l))  return false;
  }
  return true;
}

/// Compares two matrices for inequality.
bool
ann::matrix
::
operator!=
(
  const matrix& other_a
) const
{
  return !operator==(other_a);
}

/// Returns true, iff this is a Hermitian (or self-adjoint) matrix.
bool
ann::matrix
::
is_hermitian () const
{
  if (rows() != cols())  return false;

  for (size col_l = 1; col_l < cols(); ++col_l)
  {
    for (size row_l = 0; row_l < col_l; ++row_l)
    {
      if (!ann::equal((*this)(row_l, col_l), (*this)(col_l, row_l)))
        return false;
    }
  }
  return true;

}


// private methods - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//////////////////////////////////////////////////////
// SVD stuff copied from Shark (and slightly adapted)
// no more good documentation from here on
//
/// Computes the singular value decomposition of this matrix.
/// This matrix is decomposed into U_a * W_a * Vc_a, where U and V* are unitary
/// and W is a diagonal matrix.
void
ann::matrix
::
singular_value_decomposition
(
  matrix& u,
  matrix& v,
  data& w // also known as Sigma
) const
{
  diff m = rows(), n = cols();

  w.resize(n);
  data rv1(n);

  int flag;
  diff i, its, j, jj, k, l, nm(0);
  real anorm, c, f, g, h, p, s, scale, x, y, z;

  // Returns a number with the same magnitude as a and the same sign as b.
  auto adapt_sign
    = [](real a, real b)
      {
        return b > 0 ? std::fabs(a) : -std::fabs(a);
      };

  // householder reduction to bidiagonal form
  g = scale = anorm = 0.0;

  for (i = 0; i < n; i++){
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m){
      for (k = i; k < m; k++){
        scale += std::fabs(u(k, i));
      }

      if (scale != 0.0){
        for (k = i; k < m; k++){
          u(k, i) /= scale;
          s += u(k, i) * u(k, i);
        }

        f = u(i, i);
        g = -adapt_sign(std::sqrt(s), f);
        h = f * g - s;
        u(i, i) = f - g;

        for (j = l; j < n; j++){
          s = 0.0;
          for (k = i; k < m; k++){
            s += u(k, i) * u(k, j);
          }

          f = s / h;
          for (k = i; k < m; k++){
            u(k, j) += f * u(k, i);
          }
        }

        for (k = i; k < m; k++){
          u(k, i) *= scale;
        }
      }
    }

    w[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m && i != n - 1){
      for (k = l; k < n; k++){
        scale += std::fabs(u(i, k));
      }

      if (scale != 0.0){
        for (k = l; k < n; k++){
          u(i, k) /= scale;
          s += u(i, k) * u(i, k);
        }

        f = u(i, l);
        g = -adapt_sign(std::sqrt(s), f);
        h = f * g - s;
        u(i, l) = f - g;

        for (k = l; k < n; k++){
          rv1[k] = u(i, k) / h;
        }

        for (j = l; j < m; j++){
          s = 0.0;
          for (k = l; k < n; k++){
            s += u(j, k) * u(i, k);
          }

          for (k = l; k < n; k++){
            u(j, k) += s * rv1[k];
          }
        }

        for (k = l; k < n; k++){
          u(i, k) *= scale;
        }
      }
    }

    anorm = std::max(anorm, std::fabs(w[i]) + std::fabs(rv1[i]));
  }

  // accumulation of right-hand transformations
  for (l = i = n; i--; l--){
    if (l < n){
      if (g != 0.0){
        for (j = l; j < n; j++){
          // double division avoids possible underflow
          v(j, i) = (u(i, j) / u(i, l)) / g;
        }

        for (j = l; j < n; j++){
          s = 0.0;
          for (k = l; k < n; k++){
            s += u(i, k) * v(k, j);
          }

          for (k = l; k < n; k++){
            v(k, j) += s * v(k, i);
          }
        }
      }

      for (j = l; j < n; j++){
        v(i, j) = v(j, i) = 0.0;
      }
    }

    v(i, i) = 1.0;
    g = rv1[i];
  }

  // accumulation of left-hand transformations
  for (l = i = std::min(m, n); i--; l--){
    g = w[i];

    for (j = l; j < n; j++){
      u(i, j) = 0.0;
    }

    if (g != 0.0){
      g = 1.0 / g;

      for (j = l; j < n; j++){
        s = 0.0;
        for (k = l; k < m; k++){
          s += u(k, i) * u(k, j);
        }

        // double division avoids possible underflow
        f = (s / u(i, i)) * g;

        for (k = i; k < m; k++){
          u(k, j) += f * u(k, i);
        }
      }

      for (j = i; j < m; j++){
        u(j, i) *= g;
      }
    }
    else{
      for (j = i; j < m; j++){
        u(j, i) = 0.0;
      }
    }

    u(i, i)++;
  }

  // diagonalization of the bidiagonal form
  for (k = n; k--;){
    for (its = 1; its <= 30; its++){
      flag = 1;

      // test for splitting
      for (l = k + 1; l--;){
        // rv1 [0] is always zero, so there is no exit
        nm = l - 1;

        if (std::fabs(rv1[l]) + anorm == anorm){
          flag = 0;
          break;
        }

        if (std::fabs(w[nm]) + anorm == anorm){
          break;
        }
      }

      if (flag){
        // cancellation of rv1 [l] if l greater than 0
        c = 0.0;
        s = 1.0;

        for (i = l; i <= k; i++){
          f = s * rv1[i];
          rv1[i] *= c;

          if (std::fabs(f) + anorm == anorm){
            break;
          }

          g = w[i];
          h = std::hypot(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g * h;
          s = -f * h;

          for (j = 0; j < m; j++){
            y = u(j, nm);
            z = u(j, i);
            u(j, nm) = y * c + z * s;
            u(j, i) = z * c - y * s;
          }
        }
      }

      // test for convergence
      z = w[k];

      if (l == k){
        if (z < 0.0){
          w[k] = -z;
          for (j = 0; j < n; j++){
            v(j, k) = -v(j, k);
          }
        }
        break;
      }

      if (its == 30){
        throw k;
      }

      // shift from bottom 2 by 2 minor
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = std::hypot(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + adapt_sign(g, f))) - h)) / x;

      // next qr transformation
      c = s = 1.0;

      for (j = l; j < k; j++){
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g *= c;
        z = std::hypot(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;

        for (jj = 0; jj < n; jj++)
        {
          x = v(jj, j);
          z = v(jj, i);
          v(jj, j) = x * c + z * s;
          v(jj, i) = z * c - x * s;
        }

        z = std::hypot(f, h);
        w[j] = z;

        // rotation can be arbitrary if z is zero
        if (z != 0.0){
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }

        f = c * g + s * y;
        x = c * y - s * g;

        for (jj = 0; jj < m; jj++){
          y = u(jj, j);
          z = u(jj, i);
          u(jj, j) = y * c + z * s;
          u(jj, i) = z * c - y * s;
        }
      }

      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  //////////////////////////////////////////////
  //sort the eigenvalues in descending order
  for (i = 0; i < n - 1; i++){
    p = w[k = i];

    for (j = i + 1; j < n; j++){
      if (w[j] >= p){
        p = w[k = j];
      }
    }

    if (k != i){
      w[k] = w[i];
      w[i] = p;

      for (j = 0; j < n; j++){
        p = v(j, i);
        v(j, i) = v(j, k);
        v(j, k) = p;
      }

      for (j = 0; j < m; j++){
        p = u(j, i);
        u(j, i) = u(j, k);
        u(j, k) = p;
      }
    }
  }
}
