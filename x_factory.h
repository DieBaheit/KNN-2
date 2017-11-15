#ifndef ANN_X_FACOTRY_H
#define ANN_X_FACOTRY_H

#include "matrix.h"

#include <random>

namespace ann
{

class x_factory
{
 public:
  // Declare some types used by this class.
  typedef ann::matrix matrix;
  typedef matrix::real real;
  typedef matrix::size size;

  /// Constructs a factory for producing 2D vectors (2x1) with elements drawn
  /// from two different random distributions.
  x_factory
  (
    /// inclusive lower bound of the distribution of the first vector element
    real lower0_a,
    /// exclusive upper bound of the distribution of the first vector element
    real upper0_a,
    /// inclusive lower bound of the distribution of the second vector element
    real lower1_a,
    /// exclusive upper bound of the distribution of the second vector element
    real upper1_a
  );

  /// copy constructor
  x_factory (const x_factory& other_a);
  /// move constructor
  x_factory (x_factory&& other_a);

  /// Returns a random vector (according to the configured distributions).
  matrix operator() ();

 private:
  /// The type of random engine used by matrix is defined here for convenience.
  /// This engine implemnts a 64-bit Mersenne Twister as suggested by Matsumoto
  /// and Nishimura in 2000.
  typedef std::mt19937_64 random_engine;
  typedef std::uniform_real_distribution<real> distribution;

  /// This random device is used to seed the random engines.
  /// It is a non-deterministic uniform random number generator, if one is
  /// available on your system.
  static thread_local std::random_device random_device_m;

  /// This random engine is used to create values for the first vector component.
  random_engine random_engine0_m;
  /// This random engine is used to create values for the second vector component.
  random_engine random_engine1_m;

  /// Values of the first vector component are drawn from this distribution.
  distribution distribution0_m;
  /// Values of the second vector component are drawn from this distribution.
  distribution distribution1_m;
};

}

#endif //ndef ANN_X_FACOTRY_H
