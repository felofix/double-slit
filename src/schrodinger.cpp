//
//  schrodinger.cpp
//  
//
//  Created by Felix Aarekol Forseth on 20/11/2022.
//

#include "../include/schrodinger.hpp"

Schrodinger::Schrodinger(int time, double hh, double dtt){
    T = time;
    h = hh;
    dt = dtt;
    M = 1/hh - 2;
}

// Create potential matrix.
void Schrodinger::initialize_V_double(double wdx, double wx, double sdy, double soy, double v0){
    V = V.zeros(M+2, M+2);
    
    for (double i = 1; i < M + 1; i++){
        for (double j = 1; j < M + 1; j++){
            if (j*h > wx-wdx && j*h < wx + wdx){
                V(i, j) = v0;
            }
            
            if (i*h > 0.5 - sdy/2 - soy && i*h < 0.5 - sdy/2){
                V(i, j) = 0;
            }
            
            if (i*h > 0.5 + sdy/2 && i*h < 0.5 + sdy/2 + soy){
                V(i, j) = 0;
            }
        }
    }
}

int Schrodinger::matvec(int i, int j){
    int k = i + j*M;
    if (k > M*M){ // Test.
        std::cout << "Indecies doesnt match matrix size.";
        return 0;
    }
    else{
        return k;
    }
}

void Schrodinger::initialize_u(double xc, double yc, double sigmax, double sigmay, double px, double py){
    u = u.zeros(M+2,M+2);
    
    for (double i = 1; i < M + 1; i++){
        for (double j = 1; j < M + 1; j ++){
            
            arma::cx_double alpha(exp(-pow((j*h - xc),2)/(2*pow(sigmax,2))),0);
            
            arma::cx_double beta(exp(-pow((i*h - yc),2)/(2*pow(sigmay,2))),0);

            arma::cx_double gamma(cos(px*(j*h - xc)), sin(px*(j*h - xc)));

            arma::cx_double sigma(cos(py*(i*h - yc)), sin(py*(i*h - yc)));
            
            u(i, j) = alpha*beta*gamma*sigma;
        }
    }
    
    u = u/cdot(u, u);
}

void Schrodinger::evolve(){
    arma::cx_vec vecu((M+2)*(M+2), arma::fill::zeros);
    for (int i = 1; i < M + 1; i++){
        for (int j = 1; j < M + 1; j++){
            vecu(matvec(i, j)) = u(i, j);
        }
    }
    arma::cx_vec b = B*u;
    arma::cx_vec a = arma::inv(A)*b;
    
    for (int i = 1; i < (M+1)*(M+1); i++){
        u(i) = a(i);
    }
}

void Schrodinger::create_AB(){
    double r = dt/(2*pow(h, 2));
    arma::cx_vec ak(M*M, arma::fill::zeros);
    arma::cx_vec bk(M*M, arma::fill::zeros);
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            ak(matvec(i, j)) = arma::cx_double(1 + 4*r, dt/2*V(i, j));
            bk(matvec(i, j)) = arma::cx_double(1 - 4*r, -dt/2*V(i, j));
        }
    }
    initialize_A(r, ak);
    initialize_B(r, bk);
}

void Schrodinger::initialize_A(double r, arma::cx_vec a){
    A = A.zeros(M*M,M*M);
    arma::cx_mat OF(M*M, M*M, arma::fill::zeros); // Off-diagonal sadly.
    for (int i = 0; i < M*M - 1; i++){
        A(i,i) = -a[i];
        
        if ((i+1) % M != 0 && i < (M*M-1)){
            OF(i+1, i) = -r;
        }
        
        if (i < M*M - M){
            OF(i+M,i) = -r;
        }
    }
    A += OF + OF.t();
}

void Schrodinger::initialize_B(double r, arma::cx_vec b){
    B = B.zeros(M*M,M*M);
    arma::cx_mat OF(M*M, M*M, arma::fill::zeros); // Off-diagonal sadly.
    for (int i = 0; i < M*M - 1; i++){
        B(i,i) = -b[i];
        
        if ((i+1) % M != 0 && i < (M*M-1)){
            OF(i+1, i) = r;
        }
        
        if (i < M*M - M){
            OF(i+M,i) = r;
        }
    }
    B += OF + OF.t();
}

void Schrodinger::writematrixtofile(arma::mat M, std::string direc){
    // Writing to file with float values.
    std::fstream fw;
    fw.open(direc, std::fstream::app);
    if (fw.is_open())
    {
      for (int i = 0; i < M.row(0).n_elem; i++) {
          fw << M.row(i) << "\n";
      }
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
}
