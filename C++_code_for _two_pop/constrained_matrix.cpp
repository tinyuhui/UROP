#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix constrained_matrix(NumericVector p, int n){

    /* inputs: p = c(p1, p2, p3, p4) in [0,1], represent the (0,0), (0,1), (1,0) frequencies you want the simulated matrix to have
               n = number of diploid individuals
        output: 2x2n matrix with the required haplotype frequencies
    */

    // q0 is the number of alleles that will have form {0, 0}
    int q0 = 2*n*p[0];
    int q1 = 2*n*p[1];
    int q2 = 2*n*p[2];
    int q3 = 2*n - q0 - q1 - q2;

    // possible warning if n is too low
    if(abs(q3 - floor(2*n*(1 - p[0] - p[1] -p[2])))/ (2*n) > 0.1){
        Rcout << "imprecise approximation due to discretization error \n";
    }

    // sampling without replacement alleles {0, ..., 2n -1} 
    IntegerVector v = seq(0, 2*n -1);
    IntegerVector parents = sample(v, 2*n, false);
    
    // initializing empty matrix, which will contain the output 
    IntegerMatrix m(2 , 2*n);

    // unordered_map<int, IntegerVector> umap; //maps column index of the future matrix to the 2 alleles it will contain

    // giving the first q0 columns alleles (0, 0), then up to q0 + q1 alleles (0, 1), etc.
    int count = 0;
    for (IntegerVector::iterator i = parents.begin(); i != parents.end(); ++i){
        count += 1;
        if (count <= q0){
            m(_, *i) = IntegerVector::create(0, 0);
        }
        else if (count <= q0 + q1){
            m(_, *i) = IntegerVector::create(0, 1);    
        }
        else if (count <= q0 + q1 + q2){
            m(_, *i) = IntegerVector::create(1, 0);    
        }
        else{
            m(_, *i) = IntegerVector::create(1, 1);    
        }
    }
    
    return m;   
    
}