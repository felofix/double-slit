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
    M = 1/hh + 1;
    I = M - 2; // Internal points.
}

// Solver.
void Schrodinger::solve(double xc, double yc, double sigmax, double sigmay, double px, double py, double v0, int slit){
    
    if (slit == 0){
        V = V.zeros(I, I);
    }
    
    if (slit == 1){
        initialize_V_single(0.02, 0.5, 0.05, 0.05, v0);
    }
    
    if (slit == 2){
        initialize_V_double(0.02, 0.5, 0.05, 0.05, v0);
    }
    
    if (slit == 3){
        initialize_V_triple(0.02, 0.5, 0.05, 0.05, v0);
    }
    create_AB();
    initialize_u(xc, yc, sigmax, sigmay, px, py);
    double temptime = 0;
    while (temptime < T){
        writematrixtofile(u, "matrix" + std::to_string(slit) + ".txt");
        writerandcmatrixtofile(u, "matrix" + std::to_string(slit) + ".txt");
        evolve();
        temptime += dt;
    }
}

void Schrodinger::initialize_V_single(double wdx, double wx, double sdy, double soy, double v0){
    V = V.zeros(I, I);
    for (double i = 0; i < I; i++){
        for (double j = 0; j < I; j++){
            if (j*h > wx-wdx && j*h < wx + wdx){ // Fills all x values.
                V(i, j) = v0;
            }
            
            if (i*h > 0.5 - sdy/2 && i*h < 0.5 + sdy/2){ // removes y values.
                V(i, j) = 0;
            }
        }
    }
}

// Create potential matrix.
void Schrodinger::initialize_V_double(double wdx, double wx, double sdy, double soy, double v0){
    V = V.zeros(I, I);
    for (double i = 0; i < I; i++){
        for (double j = 0; j < I; j++){
            if (j*h > wx-wdx && j*h < wx + wdx){ // Fills all x values.
                V(i, j) = v0;
            }
            
            if (i*h > 0.5 - sdy/2 - soy && i*h < 0.5 - sdy/2){ // Removse y valeues.
                V(i, j) = 0;
            }
            
            if (i*h > 0.5 + sdy/2 && i*h < 0.5 + sdy/2 + soy){ // removes y values.
                V(i, j) = 0;
            }
        }
    }
}

// Create triple slit potential.
void Schrodinger::initialize_V_triple(double wdx, double wx, double sdy, double soy, double v0){
    V = V.zeros(I, I);
    for (double i = 0; i < I; i++){
        for (double j = 0; j < I; j++){
            if (j*h > wx-wdx && j*h < wx + wdx){
                V(i, j) = v0;
            }
            
            if (i*h > 0.5 - sdy/2 && i*h < 0.5 + sdy/2){ // removes y values.
                V(i, j) = 0;
            }
            
            if (i*h < 0.5 - sdy/2 - sdy && i*h > 0.5 - sdy/2 - 2*sdy){
                V(i, j) = 0;
            }
            
            if (i*h > 0.5 + sdy/2 + sdy && i*h < 0.5 + sdy/2 + 2*sdy){
                V(i, j) = 0;
            }
        }
    }
}

int Schrodinger::matvec(int i, int j){
    int k = j + i*I;
    return k;
}

void Schrodinger::initialize_u(double xc, double yc, double sigmax, double sigmay, double px, double py){
    u = u.zeros(I,I);

    for (double i = 0; i < I; i++){
        for (double j = 0; j < I; j ++){

            arma::cx_double alpha(exp(-pow((j*h - xc),2)/(2*pow(sigmax,2))),0);
            
            arma::cx_double beta(exp(-pow((i*h - yc),2)/(2*pow(sigmay,2))),0);

            arma::cx_double gamma(cos(px*(j*h - xc)), sin(px*(j*h - xc)));

            arma::cx_double sigma(cos(py*(i*h - yc)), sin(py*(i*h - yc)));
            
            u(i, j) = alpha*beta*gamma*sigma;

        }
    }
    
    u /= std::sqrt(arma::cdot(u,u));
}

void Schrodinger::evolve(){    
    arma::cx_vec vecu(I*I, arma::fill::zeros);
    
    for (int i = 0; i < I; i++){
        for (int j = 0; j < I; j++){
            vecu(matvec(i, j)) = u(i,j);
        }
    }
    
    arma::superlu_opts opts;
    opts.symmetric = true;
    arma::cx_vec a = arma::spsolve(A, B*vecu, "superlu");
    
    
    for (int i = 0; i < I; i++){
        for (int j = 0; j < I; j++){
            u(i, j) = a(matvec(i, j));
        }
    }
}

void Schrodinger::create_AB(){
    double r = dt/(2*pow(h, 2));
    arma::cx_vec ak(I*I);
    arma::cx_vec bk(I*I);
    for (int i = 0; i < I; i++){
        for (int j = 0; j < I; j++){
            ak(matvec(i, j)) = arma::cx_double(1, 4*r+dt/2*V(i, j));
            bk(matvec(i, j)) = arma::cx_double(1, -4*r-dt/2*V(i, j));
        }
    }
    initialize_A(arma::cx_double(0, r), ak);
    initialize_B(arma::cx_double(0, r), bk);
    
}

void Schrodinger::initialize_A(arma::cx_double r, arma::cx_vec a){
    arma::umat sub(2, (I-1)*I);
    arma::umat subsub(2, I*I - I);
    arma::umat sup(2, (I-1)*I);
    arma::umat supsup(2, I*I - I);
    arma::umat main(2, I*I);
    int totalr = (I-1)*I + (I-1)*I + I*I - I + I*I - I;
    arma::cx_vec rvalues(totalr, arma::fill::value(-r));
    int c = 0;   // counter one.
    int cc = 0;  // counter two.
    
    for (int i = 0; i < I*I; i++){
        main(0,i) = i;
        main(1,i) = i;
        
        if ((i+1) % I != 0 && i < (I*I-1)){
            sub(0, c) = i;
            sub(1, c) = i+1;
            sup(0, c) = i+1;
            sup(1, c) = i;
            c++;
        }
        
        if (i < I*I - I){
            subsub(0, cc) = i;
            subsub(1, cc) = i + I;
            supsup(0, cc) = i + I;
            supsup(1, cc) = i;
            cc++;
        }
    }
    
    arma::umat positionsr = arma::join_rows(arma::join_rows(sub, subsub), arma::join_rows(sup, supsup));
    arma::sp_cx_mat OD(positionsr, rvalues, I*I, I*I);
    arma::sp_cx_mat D(main, a, I*I, I*I);
    A = OD + D;
}

void Schrodinger::initialize_B(arma::cx_double r, arma::cx_vec b){
    arma::umat sub(2, (I-1)*I);
    arma::umat subsub(2, I*I - I);
    arma::umat sup(2, (I-1)*I);
    arma::umat supsup(2, I*I - I);
    arma::umat main(2, I*I);
    int totalr = (I-1)*I + (I-1)*I + I*I - I + I*I - I;
    arma::cx_vec rvalues(totalr, arma::fill::value(r));
    int c = 0;   // counter one.
    int cc = 0;  // counter two.
    
    for (int i = 0; i < I*I; i++){
        main(0,i) = i;
        main(1,i) = i;
        
        if ((i+1) % I != 0 && i < (I*I-1)){
            sub(0, c) = i;
            sub(1, c) = i+1;
            sup(0, c) = i+1;
            sup(1, c) = i;
            c++;
        }
        
        if (i < I*I - I){
            subsub(0, cc) = i;
            subsub(1, cc) = i + I;
            supsup(0, cc) = i + I;
            supsup(1, cc) = i;
            cc++;
        }
    }
    
    arma::umat positionsr = arma::join_rows(arma::join_rows(sub, subsub), arma::join_rows(sup, supsup));
    arma::sp_cx_mat OD(positionsr, rvalues, I*I, I*I);
    arma::sp_cx_mat D(main, b, I*I, I*I);
    B = OD + D;
}

void Schrodinger::writerandcmatrixtofile(arma::cx_mat M, std::string direc){
    std::ofstream wf("plotting/datafiles/real" + direc, std::fstream::app);
    
    if (wf.is_open())
    {
      wf << real(M)<< "\n";
      wf.close();
    }
    
    std::ofstream wfi("plotting/datafiles/imag" + direc, std::fstream::app);
    
    if (wfi.is_open())
    {
      wfi << imag(M)<< "\n";
      wfi.close();
    }
}

void Schrodinger::writematrixtofile(arma::cx_mat M, std::string direc){
    // Writing to file with float values.
    
    double nice = real(std::sqrt(arma::cdot(u,u)));
    std::ofstream wf("plotting/datafiles/" + direc, std::fstream::app);
    
    if (wf.is_open())
    {
      wf << real(M) % real(M) + imag(M) % imag(M) << "\n";
      wf.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
     
    wf.open("plotting/datafiles/prob" + direc, std::fstream::app);
    if (wf.is_open())
    {
      wf << std::setprecision(20) << nice << "\n";
      wf.close();
    }
    else std::cout <<"The file couldnt be opened. " << std::endl;
}
