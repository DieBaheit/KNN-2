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

  Authors: @task: Put your names here.
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
  real acc_X = 0;
  real acc_Y = 0;
  for (auto s = samps.begin(); s != samps.end(); s++) {
    acc_X += pow((s->get(0,0) - mean.get(0,0)), 2);
    acc_Y += pow((s->get(1,0) - mean.get(1,0)), 2);
  }
  matrix res(2,1);
  res.set(0,0,sqrt(acc_X/N));
  res.set(1,0,sqrt(acc_Y/N));
  return res;
}

matrix computeMean(const samples &samps){
  unsigned int N = samps.size();
  matrix acc(2,1,0);
  for (auto s = samps.begin(); s != samps.end(); s++) {
    acc += *s;
  }
  return 1.0/N * acc;
}

real normalDistribution(matrix z, matrix mean, matrix stdDev){
  real z1 = z.get(0,0);
  real z2 = z.get(1,0);
  real m1 = mean.get(0,0);
  real m2 = mean.get(1,0);
  real d1 = stdDev.get(0,0);
  real d2 = stdDev.get(1,0);
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
  for (int i = 0; i<S_l ; i++){
    samps.push_back(yFact());
  }
  // @task: Compute mean and standard deviation. You may do this here or define
  //        functions for it. When using a function, make sure to pass your
  //        data as const &.
  matrix mean = computeMean(samps);
  matrix standardDeviation = computeStandardDeviation(samps, mean);




  // @task: Create a 2D histogram (as matrix) and feed it with the S_l
  //        vectors created above.
  //        You may do this here or in a separate function (using const &).
  //        The width of a histogram bin shall be step_size_l and the lower and
  //        upper bounds must correspond to those of the x vectors.
  //        Don't forget to normalize the histogram as demanded on the sheet.
  unsigned int xSpace = (int)((upper0_l-lower0_l)/step_size_l);
  unsigned int ySpace = (int)((upper1_l-lower1_l)/step_size_l);
  matrix histogram(xSpace, ySpace, 0);
  for (auto s = samps.begin(); s != samps.end(); s++){
    unsigned int xBin = floor((s->get(0,0) - lower0_l)/step_size_l);
    unsigned int yBin = floor((s->get(1,0) - lower1_l)/step_size_l);
    if(xBin >= 0 && xBin < xSpace && yBin >= 0 && yBin < ySpace){
      histogram.set(xBin, yBin, histogram.get(xBin,yBin) + 1);
    }
  }
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
  matrix normDistrHistogr(xSpace,ySpace);
  for (int xIdx = 0; xIdx<xSpace; xIdx++){
    real xCentre = lower0_l + xIdx*step_size_l + step_size_l/2;
    for (int yIdx = 0; yIdx<ySpace; yIdx++){
      real yCentre = lower1_l + yIdx*step_size_l + step_size_l/2;
      matrix z(2,1);
      z.set(0,0,xCentre);
      z.set(1,0,yCentre);
      normDistrHistogr.set(xIdx,yIdx, normalDistribution(z, mean, standardDeviation));
    }
  }

  real errorAcc = 0;
  matrix error(xSpace,ySpace); //Nur  zum testen
  for (int xIdx = 0; xIdx<xSpace; xIdx++) {
    for (int yIdx = 0; yIdx < ySpace; yIdx++) {
      real hVal = histogram.get(xIdx, yIdx);
      real nVal = normDistrHistogr.get(xIdx, yIdx);
      real er = pow(hVal-nVal,2);
      error.set(xIdx, yIdx, er);
      errorAcc += er;
    }
  }

  std::cout<<"Error Sum: " << errorAcc <<std::endl;


  return 0;
}
