#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix generate_child_matrix(IntegerMatrix m, double c, int size){

    /* inputs: m = 2x2*n parent matrix
               c = recombination rate in [0, 0.5]
               size = number of individuals you want in child generation
        output: 2x2*size child matrix
    */

    // 2*size vector of uniform r.v.
    NumericVector unif = Rcpp::runif(2*size, 0.0, 1.0 );  //unif[i] will tell us if allele is coming from parents[i], or needs to recombine
    
    // number of parents
    int n = m.ncol()/2;  

    // parents[i] is the parent of the i-th allele in the output matrix; we are obvs sampling with replacement
    IntegerVector v = seq(0, n -1);
    IntegerVector parents = sample(v, 2*size, true);

    // this is our output matrix    
    IntegerMatrix children (2, 2*size);

    // pointer at the start of the children matrix
    IntegerVector::iterator i = children.begin();

    // pointer at the start of the sampled parents vector
    IntegerVector::iterator j = parents.begin();

    // main loop, using pointers for better performance
    for(NumericVector::iterator k = unif.begin(); k != unif.end(); ++k){

        // recombination occurs Ab  (wlog allele 0 of parent is AB, allele 1 of parent is ab)
        if(*k < c/2.0){
            *i = m(0, 2*(*j));
            ++i;
            *i = m(1, 2*(*j) + 1);
            ++i;
            ++j;
        }

        // recombination occurs aB
        else if(*k < c){
            *i = m(0, 2*(*j) + 1);
            ++i;
            *i = m(1, 2*(*j));
            ++i;
            ++j;
        }

        // no recombination AB
        else if(*k < (1 + c)/2.0){
            *i = m(0, 2*(*j));
            ++i;
            *i = m(1, 2*(*j));
            ++i;
            ++j;
        }

        // no recombination ab
        else{
            *i = m(0, 2*(*j) + 1);
            ++i;
            *i = m(1, 2*(*j) + 1);
            ++i;
            ++j;
        }
        
    }

    return children;
}