/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * test.cpp   : examples and tests of BLAS operation (a) and (d) implemented.
**/


#include "blas_routine.cpp"
#include <iostream>


int main(){
    int m = 4;
    // int n = 3;

    // Matrix matrix1 = Matrix(m, n);
    // matrix1.assign_random();
    // matrix1.print();

    // Matrix matrix2 = Matrix(n, m);
    // matrix2.assign_random();
    // matrix2.print();

    // Matrix matrix3 = Matrix::strassen(matrix1, matrix2);
    // matrix3.print();

    // Matrix matrix4 = matrix1 * matrix2;
    // matrix4.print();


    Matrix id3(m);
    id3.print();
    Matrix m1 = Matrix(m,m);
    m1.assign_random();
    m1.print();

    Matrix m2 = BLAS_3D(1, m1, id3);
    m2.print();
    return 0;
}
