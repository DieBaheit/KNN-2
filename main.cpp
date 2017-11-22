/*  Exercise 2 in Artificial Neural Networks

  Please find @task markers throughout this file, x_factory.cpp and
  y_factory.cpp to help fit your solution to the exercise into these code
  templates. Compilation of this code in its delivery state is futile.
  After the exercise has been corrected you might find further tags:
    @error: Indicates errors that led to score penalties.
    @warning: Marks errors not being penalized yet.
    @note: Advertises harmless information (mostly safe to read).

  Score:
          | 2.2b  + 2.2c  + 2.2d  + 2.2e   =   total 
  _____________________________________________________
  maximum |   2   +   4   +   3   +   5    =     14   |
  -------------------------------------------===========
  reached |   0   +   0   +   0   +   0    = ||   0   ||
  -------------------------------------------===========

  Authors: Viktor Studenyak, Fabian Albert
*/

// @note: _USE_MATH_DEFINES is now directly defined in CMakeLists.txt.

#include "matrix.h"
#include "y_factory.h"
#include "x_factory.h"

#include <cmath>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Declare short hands to the utilized types.
typedef ann::matrix matrix;
typedef matrix::real real;
typedef matrix::size size;
typedef std::vector<matrix> samples;
typedef ann::x_factory x_factory;
typedef ann::y_factory y_factory;

/*
   @task: These function signatures are suggestions (see actual tasks below).
          Your welcome to use these, but may also write your own functions.
          If you desperately want to, you may also write an extra class for
          the histogram, but it's completely over the top for this purpose.

matrix
compute_mean (const samples& samples_a)
{
}

matrix
compute_standard_deviation (const samples& samples_a, const matrix& mean_a)
{
}

matrix
create_histogram
(
  const samples& samples_a,
  const real lower0_a,
  const real upper0_a,
  const real lower1_a,
  const real upper1_a,
  const real step_size_a
)
{
}

real gauss (const matrix& z_a, const matrix& mue_a, const matrix& sigma_a)
{
}

real
compute_histogram_error
( const matrix& histogram_a,
  const matrix& mue_a, const matrix& sigma_a,
  const real lower0_a,
  const real lower1_a,
  const real step_size_a
)
{
}
*/

matrix computeStandardDeviation(const samples &samps, matrix mean){
  unsigned int N = samps.size();
  real acc_x1 = 0;
  real acc_x2 = 0;
  //computes the Sum of (xi - mean)^2 for x1 and x2
  for (auto s = samps.begin(); s != samps.end(); s++) {
    acc_x1 += pow((s->get(0,0) - mean.get(0,0)), 2);
    acc_x2 += pow((s->get(1,0) - mean.get(1,0)), 2);
  }
  matrix res(2,1);
  //calc. and store the stdDev (sqrt(sum(x1-mean)^2 / n))
  res.set(0,0,sqrt(acc_x1/N));
  res.set(1,0,sqrt(acc_x2/N));
  return res;
}

matrix computeMean(const samples &samps){
  unsigned int N = samps.size();
  matrix acc(2,1,0);
  //Compute the sum of the samples
  for (auto s = samps.begin(); s != samps.end(); s++) {
    acc += *s;
  }
  //...and divide it by the amount of samples
  return (1.0/ (real) N) * acc;
}

real normalDistribution(matrix z, matrix mean, matrix stdDev){
  //extract all necessary values from the Parameters
  real z1 = z.get(0,0);
  real z2 = z.get(1,0);
  real m1 = mean.get(0,0);
  real m2 = mean.get(1,0);
  real d1 = stdDev.get(0,0);
  real d2 = stdDev.get(1,0);
  //use them in the formula of the normal distribution
  return (1.0/(2*M_PI*d1*d2))*exp(-0.5*((pow(z1-m1,2)/pow(d1,2))+(pow(z2-m2,2)/pow(d2,2))));
}


int main(int argc, char** argv)
{
  // Take command line arguments or use default values.
  // number of y vectors
  const size S_l = (argc > 1) ? atoi(argv[1]) : 100000;
  // number of x vectors to create one y vector
  const size N_l = (argc > 2) ? atoi(argv[2]) : 10;
  // step size of the histogram
  const real step_size_l = (argc > 3) ? atof(argv[3]) : 0.1;
  // lower and upper bound for first component of x vectors
  const real lower0_l = (argc > 4) ? atof(argv[4]) : 1.0;
  const real upper0_l = (argc > 5) ? atof(argv[5]) : 3.0;
  // lower and upper bound for second component of x vectors
  const real lower1_l = (argc > 6) ? atof(argv[6]) : -2.0;
  const real upper1_l = (argc > 7) ? atof(argv[7]) : 1.0;


  // @task: Create x_factory for producing x vectors.
  x_factory xFact(lower0_l, upper0_l, lower1_l, upper1_l);
  // @task: Create y_factory (utilizing x_factory) for producing y vectors.
  y_factory yFact(xFact, N_l);

  // @task: Generate S_l y vectors with sample size N_l.
  //        To store these, you can use the supplied type samples, which
  //        is a std::vector<matrix>. This is the easier way, but your
  //        program will be slower and consume more memory (still significantly
  //        below 1GB).
  //        If you'd rather like to use a single matrix holding all the samples,
  //        your free to do so. This just involves a bit more fiddling with
  //        individual elements.
  //        (If you'd like to get all elaborate, you can also ommit storing the
  //        data by computing mean and variance online and directly fill the
  //        histogram with the samples. This allows for larger S, of course.)
  samples samps;
  //creates Samples and store them in samps
  for (unsigned int i = 0; i<S_l ; i++){
    samps.push_back(yFact());
  }
  // @task: Compute mean and standard deviation. You may do this here or define
  //        functions for it. When using a function, make sure to pass your
  //        data as const &.

  //computes the mean and stdDev
  matrix mean = computeMean(samps);
  matrix standardDeviation = computeStandardDeviation(samps, mean);




  // @task: Create a 2D histogram (as matrix) and feed it with the S_l
  //        vectors created above.
  //        You may do this here or in a separate function (using const &).
  //        The width of a histogram bin shall be step_size_l and the lower and
  //        upper bounds must correspond to those of the x vectors.
  //        Don't forget to normalize the histogram as demanded on the sheet.

  //calculate the dimensions of the histogram
  unsigned int x1_space = floor((upper0_l-lower0_l)/step_size_l);
  unsigned int x2_space = floor((upper1_l-lower1_l)/step_size_l);

  //assign each sample to a slot in the histogram
  matrix histogram(x1_space, x2_space, 0);
  for (auto s = samps.begin(); s != samps.end(); s++){
    //computes the slot coordinates for a given sample
    unsigned int x1_slot = floor((s->get(0,0) - lower0_l)/step_size_l);
    unsigned int x2_slot = floor((s->get(1,0) - lower1_l)/step_size_l);
    // if the computed slot isn't in the histograms range then ignore the sample
    // (if step_size doesn't divide (upper - lower) then the last part(that is smaller then
    // step_size) isn't part of the histogram. Through the sample might fall in this area, so we
    // must test it)
    if(/*x1_slot >= 0 && */ x1_slot < x1_space && /*x2_slot >= 0 &&*/ x2_slot < x2_space){
      histogram.set(x1_slot, x2_slot, histogram.get(x1_slot,x2_slot) + 1);
    }
  }
  // normalize the histogram
  histogram *= 1.0/(S_l * step_size_l * step_size_l);



  // @task: Compare the histogram to the normal distribution by computing the
  //        sum of squared errors (SSE).
  //        You most likely want to write a separate function for this, as
  //        sepcified in 2.2e. Use the mean and standard deviation computed
  //        above and the centers of the histogram bins as parameters to the
  //        normal distribution. Iterate over all histogram bins and sum up the
  //        squared distance between the value in the histogram bin and the
  //        theoretical value of the normal distribution.
  //        Write the single computed error value to std::cout. Be astonished.

  // Computes a second histogram, with the values of the normal distribution (under
  // consideration of the memory used to save all samples, this little waste of memory
  // can be accepted :) )
  matrix normDistrHistogr(x1_space,x2_space);
  for (unsigned int x1_idx = 0; x1_idx<x1_space; x1_idx++){
    real x1_centre = lower0_l + x1_idx*step_size_l + step_size_l/2;
    for (unsigned int x2_idx = 0; x2_idx<x2_space; x2_idx++){
      real x2_centre = lower1_l + x2_idx*step_size_l + step_size_l/2;
      matrix z(2,1);
      z.set(0,0,x1_centre);
      z.set(1,0,x2_centre);
      normDistrHistogr.set(x1_idx,x2_idx, normalDistribution(z, mean, standardDeviation));
    }
  }

  // Sum up the squared error values for each slot of the histogram
  real errorAcc = 0;
  for (unsigned int x1_idx = 0; x1_idx<x1_space; x1_idx++) {
    for (unsigned int x2_idx = 0; x2_idx < x2_space; x2_idx++) {
      real hVal = histogram.get(x1_idx, x2_idx);
      real nVal = normDistrHistogr.get(x1_idx, x2_idx);
      real er = pow(hVal-nVal,2);
      errorAcc += er;
    }
  }

  std::cout<<"Error Sum: " << errorAcc <<std::endl;


  return 0;
}
