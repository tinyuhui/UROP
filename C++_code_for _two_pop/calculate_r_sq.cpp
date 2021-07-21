#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double calculate_r_sq(IntegerMatrix m){

    /* inputs: a 2xn binary matrix
       output: r^2 for the matrix 
    */

    int n = m.ncol();
    NumericVector p (4); // vector of frequencies

    //calculating counts
    IntegerMatrix::iterator i = m.begin();
    while(i < m.end()) {
        if ( (*i == 0) & (*(i+1) == 0)) p[0] += 1;
        else if ((*i == 0) & (*(i+1) == 1)) p[1] += 1;
        else if ((*i == 1) & (*(i+1) == 0)) p[2] += 1;
        else p[3] += 1;
        ++i;
        ++i;
    }

    // transforming counts into frequencies
    p =  p/n;


    double raw_LD = p[0]*p[3] - p[1]*p[2];
    double r_sq = (pow(raw_LD, 2))/((p[0] + p[1])*(p[2] + p[3])*(p[0] + p[2])*(p[1] + p[3]));
    return r_sq;
}