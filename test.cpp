/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * test.cpp   : examples and tests of BLAS operation (a) and (d) implemented.
**/


#include "matrix.hpp"
#include <iostream>


int main(){
    int m = 4;
    int n = 7;

    Matrix matrix1 = Matrix(m, m);
    matrix1.assign_random();
    matrix1.print();

    Matrix matrix2 = Matrix(m, m);
    matrix2.assign_random();
    matrix2.print();

    Matrix matrix3 = strassen(matrix1, matrix2);
    matrix3.print();

    Matrix matrix4 = matrix1 * matrix2;
    matrix4.print();
    return 0;
}
