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
    int size = 2; // When a M size is specified, you need to subtract 2.
    int time = 10;
    Schrodinger test(size, time);
    arma::cx_vec A = {1, 2, 3};
    test.initialize_B(1.0, A);
    return 0;
}

