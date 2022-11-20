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
    int size = 3; // When a M size is specified, you need to subtract 2.
    int time = 10;
    Schrodinger test(size, time);
    test.initialize_u(0.5, 0.5, 0.05, 0.05, 200, 0);
    return 0;
}

