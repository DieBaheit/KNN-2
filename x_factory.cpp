#include "x_factory.h"

thread_local std::random_device ann::x_factory::random_device_m;

/// Constructs a factory for producing 2D vectors (2x1) with elements drawn
/// from two different random distributions.
ann::x_factory
::
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
)
: random_engine0_m(random_device_m()),
  random_engine1_m(random_device_m()),
  // @task: Call the constructors of the two distributions with the passed values.
  //        It's straight forward, but if unsure, you can take a look at
  //        function fill_with_uniform_samples of class matrix.
  //distribution0_m(  todo  ),
  //distribution1_m(  todo  )
{
}

/// copy constructor
ann::x_factory
::
x_factory
(
  const x_factory& other_a
)
: random_engine0_m(random_device_m()),
  random_engine1_m(random_device_m()),
  distribution0_m(other_a.distribution0_m),
  distribution1_m(other_a.distribution1_m)
{
}

/// move constructor
ann::x_factory
::
x_factory
(
  x_factory&& other_a
)
: random_engine0_m(random_device_m()),
  random_engine1_m(random_device_m()),
  distribution0_m(std::move(other_a.distribution0_m)),
  distribution1_m(std::move(other_a.distribution1_m))
{
}


/// Returns a random vector (according to the configured distributions).
ann::x_factory::matrix
ann::x_factory
::
operator() ()
{
  // @task: Create a 2x1 matrix; set the elements to values from the respective
  //        distributions using the respective engines; return the matrix.
  //        Function fill_with_uniform_samples of class matrix shows how to get
  //        a random value from a distribution and an engine.
}

