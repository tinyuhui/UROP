#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix merge_on_indices (IntegerMatrix m1, IntegerMatrix m2, IntegerVector index1, IntegerVector index2){

    /* inputs: m1, m2 matrices to merge
       inputs: index1 indicates which individuals from m1 to use
               for ex. if 2 is in index1, then the two gametes corresponding to individual 2 (i.e  m[:, 2*2], m[:, 2*2 +1])
                                          will be put in the output matrix
        output: a matrix which takes the individuals from m1, m2 according to index1, index2
    */     

    int len1 = index1.length();
    int len2 = index2.length();
    int n = len1 + len2;

    // we have n individuals, so matrix must have ncol = 2*n
    IntegerMatrix out (2, 2*n);

    // using pointers for better performance, first 2*len1 columns will come from m1
    // last 2*len2 columns will come from m2 
    IntegerMatrix::iterator i = out.begin();
    for (IntegerVector::iterator j = index1.begin(); j != index1.end(); ++j){
        *i = m1(0, 2*(*j));
        ++i;
        *i = m1(1, 2*(*j));
        ++i;
        *i = m1(0, 2*(*j) + 1);
        ++i;
        *i = m1(1, 2*(*j) + 1);
        ++i;
    }
    for(IntegerVector::iterator j = index2.begin(); j != index2.end(); ++j){
        *i = m2(0, 2*(*j));
        ++i;
        *i = m2(1, 2*(*j));
        ++i;
        *i = m1(0, 2*(*j) + 1);
        ++i;
        *i = m1(1, 2*(*j) + 1);
        ++i;
    }
    return out;    
}