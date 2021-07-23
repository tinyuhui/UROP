#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;

/*
// [[Rcpp::export]]

IntegerVector samp(int n, int size){
    IntegerVector vec = seq(0, n -1);
    IntegerVector out = sample(vec, size, true);
    double un = R::runif(0,1);
    if (un < 0.99){
        Rcout << "The value of v : " << un << "\n";
    }

    return out;
}

// [[Rcpp::export]]

vector<vector<int>> provetta(int n){
    vector<vector<int>> vec{ { 1, 2, 3 }, 
                         { 4, 5, 6 }, 
                         { 7, 8, 9, 4 } }; 
    return vec;
}


// [[Rcpp::export]]
double rcpp_sum1(NumericMatrix x) {
  double total = 0;
  for(NumericMatrix::iterator i = x.begin(); i != x.end(); ++i) {
    //Rprintf("%i", *i); 
    Rcout << "p" << *i;
  }
  return total;
}

// [[Rcpp::export]]
double rcpp_sum2(NumericMatrix x) {
  double total = 0;
  NumericMatrix::iterator i = x.begin();
  while(i < x.end()) {
    Rcout << "f" << *i;
    i++;
    Rcout << "s" << *i;
    i++;
  }
  return total;
}
*/

// [[Rcpp::export]]
IntegerMatrix constrained_matrix(NumericVector p, int n){
    int q0 = 2*n*p[0];
    int q1 = 2*n*p[1];
    int q2 = 2*n*p[2];
    int q3 = 2*n - q0 - q1 - q2;
    Rcout << q0 <<" "<< q1 <<" "<< q2 << " " << q3 << "\n";
    if(abs(q3 - floor(2*n*(1 - p[0] - p[1] -p[2])))/ (2*n) > 0.1){
        Rcout << "imprecise approximation due to discretization error \n";
    }
    IntegerVector v = seq(0, 2*n -1);
    IntegerVector parents = sample(v, 2*n, false);
    Rcout << parents << "\n";
    IntegerMatrix m(2 , 2*n);

    // unordered_map<int, IntegerVector> umap; //maps column index of the future matrix to the 2 alleles it will contain
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


