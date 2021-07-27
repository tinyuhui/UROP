#include <R.h>
#include <Rcpp.h>
#include <cmath>    /* needed for floor() */
using namespace Rcpp;
using namespace std;





// [[Rcpp::export]]
List simulation(NumericVector p0, NumericVector q0, IntegerVector N, IntegerVector M, double c, NumericMatrix migration){

    /*

    inputs: p0: 4x1 vector; initial haplotypes frequencies at generation 0 for population 1
    inputs: q0: 4x1 vector; initial haplotypes frequencies at generation 0 for population 2
    inputs: N: vector representing the size of population 1 at generation i, BEFORE MIGRATION OCCURS
    inputs: M: vector representing the size of population 2 at generation i, BEFORE MIGRATION OCCURS
            M and N should have the same size, which is the number of generations to simulate
    inputs: c is the recombination rate for the two loci
    inputs: migration: FIXED migration matrix; migration[0,1] represents the proportion of population 1 that migrates
                       to population 2 at each generation
  
    output: list of size 6 containing: simulated haplotypes frequencies for the two populations for all generations
                                     simulated D and r^2 values for the two populations for all generations

    */

    // number of generations
    int t = N.length();

    // initializing matrices of haplotype frequencies (to be filled during the loop)
    NumericMatrix p (4, t);
    NumericMatrix q (4, t);

    // initializing vectors of raw LD measures (length = number of generations)
    NumericVector D1 (t);
    NumericVector D2 (t);

    // initializing vectors of r^2 measures (length = number of generations)
    NumericVector r_sq1 (t);
    NumericVector r_sq2 (t);

    // calculating the migration flows at each generation
    // ex. size2[i] indicates how many individuals will migrate from pop1 to pop2 between generation i and generation i+1
    IntegerVector size1 = migration_vec(N, migration, 0, 0);   //stay in pop1
    IntegerVector size2 = N - size1;
    IntegerVector size3 = migration_vec(M, migration, 1, 1);   //stay in pop2
    IntegerVector size4 = M - size3;
    /*
    Rcout << size3 << "\n";
    Rcout << size1 << "\n";
    */

    // initialization: old_matrix1 will contain the parent generation 1 through the loop, new_matrix will containg the child generation
    IntegerMatrix old_matrix1 = constrained_matrix(p0, N[0]);
    IntegerMatrix old_matrix2 = constrained_matrix(q0, M[0]);
    IntegerMatrix new_matrix1( 2 , 3 );  // this is just a random initialization, to avoid having multiple objects 
    IntegerMatrix new_matrix2( 2 , 2);   // with same name but different pointers (which is what happens if we define new_matrix within the loop,
                                         // at every iteration, we would define a new object; in this way at each iteration we just modifiy this object)
    
    /*
    Rcout << "before loop" << "\n";
    Rcout << old_matrix1 << "\n";
    Rcout << old_matrix2 << "\n";
    */
    
    // calculating the 3 statistics for generation 0
    List stats_pop1 = calculate_r_sq(old_matrix1);
    List stats_pop2 = calculate_r_sq(old_matrix2);
     
    // stats_pop1[0];

    // updating the 3 statistics for generation 0
    p(_, 0) = static_cast<NumericVector>(stats_pop1[0]);
    q(_, 0) = static_cast<NumericVector>(stats_pop2[0]);  

    D1[0] = stats_pop1[1];
    D2[0] = stats_pop2[1];

    r_sq1[0] = stats_pop1[2];
    r_sq2[0] = stats_pop2[2];
    

    for (int i = 1; i < t; ++i){
        
        // sampling which individuals will stay in pop1, which will migrate to pop2
        IntegerVector v = seq(0, N[i - 1] - 1 );
        IntegerVector temp1 = sample(v, v.length(), false);
        IntegerVector pop_1_stays = temp1[Range(0, size1[i- 1] - 1)];
        IntegerVector pop_1_migrates = temp1[Range(size1[i - 1], v.length() - 1)];

        // same thing for pop2
        v = seq(0, M[i - 1] - 1 );                                             
        IntegerVector temp2 = sample(v, v.length(), false);
        IntegerVector pop_2_stays = temp2[Range(0, size3[i - 1] - 1)];
        IntegerVector pop_2_migrates = temp2[Range(size3[i - 1], v.length() - 1)];

        /*
        Rcout << "migration vecs" << "\n";
        Rcout << pop_1_stays << "\n";
        Rcout << pop_1_migrates << "\n";
        Rcout << pop_2_stays << "\n";
        Rcout << pop_2_migrates << "\n";

        Rcout << "old matrix" << "\n";
        Rcout << old_matrix1 << "\n";
        Rcout << old_matrix2 << "\n";
        */

        // creating the new populations after migration occurs, i.e. putting pop_1_stays and pop_2_migrates together etc.
        new_matrix1 = merge_on_indices(old_matrix1, old_matrix2, pop_1_stays, pop_2_migrates);
        new_matrix2 = merge_on_indices(old_matrix1, old_matrix2, pop_1_migrates, pop_2_stays);

        /*
        Rcout << "After migration, before reproduction" << "\n";
        Rcout << new_matrix1 << "\n";
        Rcout << new_matrix2 << "\n";
        */

        // making the new populations reproduce, size of the populations are specified by N[i], M[i]
        new_matrix1 = generate_child_matrix(new_matrix1, c, N[i]);    
        new_matrix2 = generate_child_matrix(new_matrix2, c, M[i]);

        /*
        Rcout << "After migration, after reproduction" << "\n";
        Rcout << new_matrix1 << "\n";
        Rcout << new_matrix2 << "\n";
        */

        // updating old matrix to generation i
        old_matrix1 = new_matrix1;
        old_matrix2 = new_matrix2;

        // calculating 3 statistics for generation i
        List stats_pop1 = calculate_r_sq(old_matrix1);
        List stats_pop2 = calculate_r_sq(old_matrix2);

        // updating the 3 statistics in the output vectors
        p(_,i) = static_cast<NumericVector>(stats_pop1[0]);
        q(_,i) = static_cast<NumericVector>(stats_pop2[0]);

        D1[i] = stats_pop1[1];
        D2[i] = stats_pop2[1];

        r_sq1[i] = stats_pop1[2];
        r_sq2[i] = stats_pop2[2];


    }
    
    List out = List::create(p, q, D1, D2, r_sq1, r_sq2);
    return out;

}