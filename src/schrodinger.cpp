//
//  schrodinger.cpp
//  
//
//  Created by Felix Aarekol Forseth on 20/11/2022.
//

#include "../include/schrodinger.hpp"

Schrodinger::Schrodinger(int size, int time){
    M = size;
    T = time;
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
            
            arma::cx_double alpha(exp(-pow((i/M - xc),2)/(2*pow(sigmax,2))),0);
            
            arma::cx_double beta(exp(-pow((j/M - yc),2)/(2*pow(sigmay,2))),0);

            arma::cx_double gamma(cos(px*(i/M - xc)), sin(px*(i/M - xc)));

            arma::cx_double sigma(cos(py*(j/M - yc)), sin(py*(j/M - yc)));
            
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

void Schrodinger::create_AB(double h, double dt, arma::mat V){
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
