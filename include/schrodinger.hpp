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
    arma::cx_mat u;
    
    // Declaration function.
    Schrodinger(int M, int T);
    
    // Switch between matrix indexes ij and vector position k.
    int matvec(int i, int j);
    
    // Evolve in time.
    void evolve();
    
    // Initilizing space.
    void initialize_u(double xc, double yc, double sigmax, double sigmay, double px, double py);
    
    // Creating starting conditions.
    void create_AB(double h, double dt, arma::mat V);
    
    // Initilizing A matrix.
    void initialize_A(double r, arma::cx_vec a);
    
    // Initilizing B matrix. 
    void initialize_B(double r, arma::cx_vec b);
};


#endif /* schrodinger_hpp */
