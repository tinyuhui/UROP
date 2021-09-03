#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
IntegerVector migration_vec(IntegerVector N, NumericMatrix m, int index1, int index2){

    /* inputs:  N = vector of population sizes for generations 1, 2, 3,...
                m = migration matrix
                index1, index2 indicate which entry of m you want. i.e. migration rate from population "index1" (which has size N) to population "index2"
        output: a vector like N, that indicates at each generation the number of individuals going 
                from population "index1" to population "index2"
    */

    IntegerVector out (N.length());
    double val = m(index1, index2);

    IntegerVector::iterator i = N.begin();
    for (IntegerVector::iterator j = out.begin(); j != out.end(); ++j){
        *j = floor(val * (*i));
        ++i;
    }


    return out;
}