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
    double size = 0.005; // When a M size is specified, you need to subtract 2.
    int time = 10;
    Schrodinger test(time, size, 0.1);
    test.initialize_V_double(0.02, 0.5, 0.05, 0.05, 1);
    test.initialize_u(0.25, 0.5, 0.05, 0.05, 200, 0);

    test.writematrixtofile(real(test.u), "test_matrix_real.txt");
    test.writematrixtofile(imag(test.u), "test_matrix_imag.txt");
    return 0;
}

