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
    arma::mat V;
    arma::cx_mat u;
    double h;
    double dt;
    
    // Declaration function.
    Schrodinger(int T, double h, double dt);
    
    // Switch between matrix indexes ij and vector position k.
    int matvec(int i, int j);
    
    // Evolve in time.
    void evolve();
    
    // Initilizing space.
    void initialize_u(double xc, double yc, double sigmax, double sigmay, double px, double py);
    
    // Creating starting conditions.
    void create_AB();
    
    // Create potential matrix.
    void initialize_V_double(double wdx, double wx, double sdy, double soy, double v0);
    
    // Initilizing A matrix.
    void initialize_A(double r, arma::cx_vec a);
    
    // Initilizing B matrix. 
    void initialize_B(double r, arma::cx_vec b);
    
    // Writing matrix to file
    void writematrixtofile(arma::mat M, std::string direc);
};


#endif /* schrodinger_hpp */
