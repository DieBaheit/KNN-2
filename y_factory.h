#ifndef ANN_Y_FACOTRY_H
#define ANN_Y_FACOTRY_H

#include "x_factory.h"

namespace ann
{

class y_factory
{
 public:
  // Declare some types used by this class.
  typedef x_factory::matrix matrix;
  typedef matrix::real real;
  typedef matrix::size size;

  /// Constructs a factory for creating random vectors, which are averaged over
  /// N samples of vectors created by the given factory.
  y_factory (x_factory x_factory_a, size N_a);

  /// Returns a random vector (according to the used factory and sample size).
  matrix operator() ();

 private:
  /// factory for creating sample vectors to average over
  x_factory x_factory_m;
  /// number of vectors to average over
  size N_m;
};

}

#endif //ndef ANN_Y_FACOTRY_H
