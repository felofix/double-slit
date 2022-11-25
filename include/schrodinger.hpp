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
    double T; // Time of evolving.
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    arma::mat V;
    arma::cx_mat u;
    double h;
    double dt;
    
    // Declaration function.
    Schrodinger(double T, double h, double dt);
    
    // Solver.
    void solve(double xc, double yc, double sigmax, double sigmay, double px, double py, double v0,int slit);
    
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
    void initialize_A(arma::cx_double r, arma::cx_vec a);
    
    // Initilizing B matrix. 
    void initialize_B(arma::cx_double r, arma::cx_vec b);
    
    // Writing matrix to file
    void writematrixtofile(arma::cx_mat M, std::string direc);
};


#endif /* schrodinger_hpp */
