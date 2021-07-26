#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;

/* Main idea here was, after sampling the parents, instead that every time having to lookup their corresponding 
   value ex. matr[:,12] which takes o(n);    storing in a hashtable  parent -> allele.  Hence when building 
   children matrix lookup time is o(1)

/*
// [[Rcpp::export]]
IntegerMatrix constrained_matrix(NumericVector p, int n){
    int q0 = 2*n*p[0];
    int q1 = 2*n*p[1];
    int q2 = 2*n*p[2];
    int q3 = 2*n - q0 - q1 - q2;

    if(abs(q3 - floor(2*n*(1 - p[0] - p[1] -p[2])))/ (2*n) > 0.1){
    Rcout << "imprecise approximation due to discretization error \n";
    
    // IntegerVector v = seq(0, 2*n -1);
    IntegerVector parents = sample(2*n, 2*n);

    IntegerMatrix m2(2 , 2*n);

    unordered_map<int, IntegerVector> umap; //maps column index of the future matrix to the 2 alleles it will contain
    int count = 0;
    for (IntegerVector::iterator i = parents.begin(); i != parents.end(); ++i){
        count += 1;
        if (count <= q1){
            const int ptr = *i ;
            IntegerVector umap[1] = IntegerVector::create(0, 0); 
        }
    }
    
    return m2;   
    }
}
*/

// [[Rcpp::export]]
int pp(int n){
    IntegerVector parents = sample(2*n, 2*n);
    unordered_map<int, IntegerVector> umap;

    int count = 0;
    for (IntegerVector::iterator i = parents.begin(); i != parents.end(); ++i){
        count += 1;
        if (count <= 6){
            const int ptr = *i ;
            umap[ptr] = IntegerVector::create(0, 0); 
            Rcout << umap[ptr];   // tried same thing with *i directly but still didn't work. Initially I was getting Visual Studio 
            //to give me a warning that the index of umap had to be constant. Then it disappeared but I still get errors when compiling it using Rcpp
        }
    return 1;
}