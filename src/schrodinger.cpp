    //
//  schrodinger.cpp
//  
//
//  Created by Felix Aarekol Forseth on 20/11/2022.
//

#include "../include/schrodinger.hpp"

Schrodinger::Schrodinger(double time, double hh, double dtt){
    T = time;
    h = hh;
    dt = dtt;
    M = 1/hh;
}

// Solver.
void Schrodinger::solve(double xc, double yc, double sigmax, double sigmay, double px, double py, double v0, int slit){
    
    // if (slit == 1){}; single slit
    
    if (slit == 2){
        initialize_V_double(0.02, 0.5, 0.05, 0.05, v0);
    }
    
    // if (slit == 3){}; single slit
    
    create_AB();
    initialize_u(xc, yc, sigmax, sigmay, px, py);
    double temptime = 0;
    while (temptime < T){
        writematrixtofile(u, "matrix" + std::to_string(slit) + ".txt");
        evolve();
        temptime += dt;
    }
}

// Create potential matrix.
void Schrodinger::initialize_V_double(double wdx, double wx, double sdy, double soy, double v0){
    V = V.zeros(M, M);
    for (double i = 1; i < M -1; i++){
        for (double j = 1; j < M -1; j++){
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
    if (k > (M)*(M)){ // Test.
        std::cout << "Indecies doesnt match matrix size.";
        return 0;
    }
    else{
        return k;
    }
}

void Schrodinger::initialize_u(double xc, double yc, double sigmax, double sigmay, double px, double py){
    u = u.zeros(M,M);
    
    for (double i = 1; i < M - 1; i++){
        for (double j = 1; j < M - 1; j ++){
            
            arma::cx_double alpha(exp(-pow((j*h - xc),2)/(2*pow(sigmax,2))),0);
            
            arma::cx_double beta(exp(-pow((i*h - yc),2)/(2*pow(sigmay,2))),0);

            arma::cx_double gamma(cos(px*(j*h - xc)), sin(px*(j*h - xc)));

            arma::cx_double sigma(cos(py*(i*h - yc)), sin(py*(i*h - yc)));
            
            u(i, j) = alpha*beta*gamma*sigma;
        }
    }
    u = u/arma::accu(arma::abs(u));
}

void Schrodinger::evolve(){
    arma::cx_vec vecu(M*M, arma::fill::zeros);
    
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            vecu(matvec(i, j)) = u(i, j);
        }
    }
    
    arma::cx_vec b = B*vecu;
    arma::superlu_opts opts;
    opts.symmetric = true;
    arma::cx_vec a = arma::spsolve(A, b, "superlu", opts);
    
    for (int i = 0; i < (M)*(M); i++){
        u(i) = a(i);
        u.col(0).fill(arma::cx_double(0,0));
        u.row(0).fill(arma::cx_double(0,0));
        u.col(M-1).fill(arma::cx_double(0,0));
        u.row(M-1).fill(arma::cx_double(0,0));
    }
}

void Schrodinger::create_AB(){
    double r = dt/(2*pow(h, 2));
    arma::cx_vec ak(M*M);
    arma::cx_vec bk(M*M);
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            ak(matvec(i, j)) = arma::cx_double(1, 4*r+dt/2*V(i, j));
            bk(matvec(i, j)) = arma::cx_double(1, -4*r-dt/2*V(i, j));
        }
    }
    initialize_A(arma::cx_double(0, r), ak);
    initialize_B(arma::cx_double(0, r), bk);
}

void Schrodinger::initialize_A(arma::cx_double r, arma::cx_vec a){
    arma::umat sub(2, (M-1)*M);
    arma::umat subsub(2, M*M - M);
    arma::umat sup(2, (M-1)*M);
    arma::umat supsup(2, M*M - M);
    arma::umat main(2, M*M);;
    int totalr = (M-1)*M + (M-1)*M + M*M - M + M*M - M;
    arma::cx_vec rvalues(totalr, arma::fill::value(-r));
    int c = 0;   // counter one.
    int cc = 0;  // counter two.
    
    for (int i = 0; i < M*M; i++){
        main(0,i) = i;
        main(1,i) = i;
        
        if ((i+1) % M != 0 && i < (M*M-1)){
            sub(0, c) = i;
            sub(1, c) = i+1;
            sup(0, c) = i+1;
            sup(1, c) = i;
            c++;
        }
        
        if (i < M*M - M){
            subsub(0, cc) = i;
            subsub(1, cc) = i + M;
            supsup(0, cc) = i + M;
            supsup(1, cc) = i;
            cc++;
        }
    }
    
    arma::umat positionsr = arma::join_rows(arma::join_rows(sub, subsub), arma::join_rows(sup, supsup));
    arma::sp_cx_mat OD(positionsr, rvalues, M*M, M*M);
    arma::sp_cx_mat D(main, a, M*M, M*M);
    A = OD + D;
}

void Schrodinger::initialize_B(arma::cx_double r, arma::cx_vec b){
    arma::umat sub(2, (M-1)*M);
    arma::umat subsub(2, M*M - M);
    arma::umat sup(2, (M-1)*M);
    arma::umat supsup(2, M*M - M);
    arma::umat main(2, M*M);
    int totalr = (M-1)*M + (M-1)*M + M*M - M + M*M - M;
    arma::cx_vec rvalues(totalr, arma::fill::value(r));
    int c = 0;   // counter one.
    int cc = 0;  // counter two.
    
    for (int i = 0; i < M*M; i++){
        main(0,i) = i;
        main(1,i) = i;
        
        if ((i+1) % M != 0 && i < (M*M-1)){
            sub(0, c) = i;
            sub(1, c) = i+1;
            sup(0, c) = i+1;
            sup(1, c) = i;
            c++;
        }
        
        if (i < M*M - M){
            subsub(0, cc) = i;
            subsub(1, cc) = i + M;
            supsup(0, cc) = i + M;
            supsup(1, cc) = i;
            cc++;
        }
    }
    
    arma::umat positionsr = arma::join_rows(arma::join_rows(sub, subsub), arma::join_rows(sup, supsup));
    arma::sp_cx_mat OD(positionsr, rvalues, M*M, M*M);
    arma::sp_cx_mat D(main, b, M*M, M*M);
    B = OD + D;
}

void Schrodinger::writematrixtofile(arma::cx_mat M, std::string direc){
    // Writing to file with float values.
    arma::mat MM = arma::abs(M);
    std::fstream fw;
    fw.open(direc, std::fstream::app);
    if (fw.is_open())
    {
      for (int i = 0; i < MM.row(0).n_elem; i++) {
          fw << MM.row(i) << "\n";
      }
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
    
    fw.open("prob" + direc, std::fstream::app);
    if (fw.is_open())
    {
      fw << arma::accu(MM) << "\n";
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
}
