#include "y_factory.h"

/// Constructs a factory for creating random vectors, which are averaged over
/// N samples of vectors created by the given factory.
ann::y_factory
::
y_factory
(
  x_factory x_factory_a,
  size N_a
)
: x_factory_m(std::move(x_factory_a)),
  N_m(N_a)
{
}

/// Returns a random vector (according to the used factory and sample size).
ann::y_factory::matrix
ann::y_factory
::
operator() ()
{
  // @task: Create a 2x1 matrix y and assign to it the average of N_m matrices x
  //        created by x_factory_m. Return the thus create vector y.
  //        You get a matrix from the factory by calling:   x_factory_m()
    matrix acc(2,1,0); //mit 0 initialisiert
    for(int n = 0 ; i < N_m ; i++){
        acc += x_factory_m();
    }
    return (1.0/N_m) * acc;
}


