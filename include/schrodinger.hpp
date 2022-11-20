//
//  schrodinger.hpp
//  
//
//  Created by Felix Aarekol Forseth on 20/11/2022.
//

#ifndef schrodinger_hpp
#define schrodinger_hpp

#include <stdio.h>
#include <armadillo>

class Schrodinger{
public:
    int M; // Size of box.
    int T; // Time of evolving.
    arma::cx_mat A;
    arma::cx_mat B;
    
    // Declaration function.
    Schrodinger(int M, int T);
    
    // Switch between matrix indexes ij and vector position k.
    int matvec(int i, int j);
    
    // Initilizing A matrix.
    void initialize_A(double r, arma::cx_vec a);
    
    // Initilizing B matrix. 
    void initialize_B(double r, arma::cx_vec b);
};


#endif /* schrodinger_hpp */
