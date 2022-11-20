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
