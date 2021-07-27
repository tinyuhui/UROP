#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;



// [[Rcpp::export]]
List average_r_sq_individual_based(NumericVector p0, NumericVector q0, IntegerVector N, IntegerVector M, double c, NumericMatrix migration, int simulations){
    
    // This function simply takes the average r^2 over a given number of simulations
    // Returns a list of the two average r^2 for pop 1 and pop 2
    
    
    NumericVector r_sq_1 (N.length());
    NumericVector r_sq_2 (N.length());


    for (int n = 0; n < simulations; ++n){

        List l_temp = simulation(p0, q0, N, M, c, migration);
        NumericVector out1 = l_temp[4];
        NumericVector out2 = l_temp[5];    
           
        r_sq_1 += out1;
        r_sq_2 += out2;
 
        
    }