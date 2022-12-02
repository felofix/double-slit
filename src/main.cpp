//
//  main.cpp
//  
//
//  Created by Felix Aarekol Forseth on 20/11/2022.
//

#include <stdio.h>
#include <iostream>
#include <armadillo>
#include "../include/schrodinger.hpp"

int main(int argc, const char * argv[]) {
    if (argc != 12){ // Checking if there is enough command-line.argumets.
        std::cout << "You have entered to few arguments." << std::endl;
        return 0;
    }
    
    double size = atof(argv[1]);
    double dt = atof(argv[2]);
    double time = atof(argv[3]);
    double xc = atof(argv[4]);
    double yc = atof(argv[5]);
    double sigmax = atof(argv[6]);
    double sigmay = atof(argv[7]);
    double px =  atof(argv[8]);
    double py = atof(argv[9]);
    double v0 = atof(argv[10]);
    int slit = atoi(argv[11]);
    Schrodinger test(time, size, dt);
    test.solve(xc, yc, sigmax, sigmay, px, py, v0, slit);
    return 0;
}

