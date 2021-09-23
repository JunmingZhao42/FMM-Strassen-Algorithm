/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * test.cpp   : examples and tests of BLAS operation (a) and (d) implemented.
**/


#include "matrix1.h"
#include <iostream>


int main(){
    int m = 6;
    int n = 7;
    
    Matrix matrix1 = Matrix(m, n);
    matrix1.assign_random();
    matrix1.print_matrix();

    Matrix matrix2 = matrix1.slice(1, 4, 3, 0);
    // double ** subdata = new double*[3];
    
    // for(int i=0; i<3; i++){
    //     std::cout << matrix1.data[i] << std::endl;
    //     subdata[i] = matrix1.data[i];
    // }

    matrix2.print_matrix();

    //Matrix matrix3 = matrix1 + matrix2;
    //matrix3.print_matrix();

    return 0;
}
