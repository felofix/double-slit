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
    double size = 0.005;
    double time = 0.008;
    double dt = 2.5e-5;
    Schrodinger test(time, size, dt);
    test.solve(0.25, 0.5, 0.05, 0.05, 200, 0, 0, 2);
    return 0;
}

