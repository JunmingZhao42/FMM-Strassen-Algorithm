/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * test.cpp   : examples and tests of BLAS operation (a) and (d) implemented.
**/


#include "blas_routine.cpp"
#include <iostream>
#include <chrono>


int main(){

    for (int m=100; m < 500; m+=10){
        std::cout << m << ",";
        // create matrix A
        Matrix A = Matrix(m);
        A.assign_random();

        // create matrix B
        Matrix B = Matrix(m);
        B.assign_random();

        // create matrix C
        Matrix C(m);

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // C = 2.0*A*B + 2.0*C
        BLAS_3A(2.0, 2.0, A, B, C);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << ",";

        // standard matrix multiplication
        begin = std::chrono::steady_clock::now();
        C.assign_zeros();
        // C = 2.0*A*B + 2.0*C
        C *= 2.0;
        C += A*B;
        end = std::chrono::steady_clock::now();

        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
    }

    
    // Matrix matrix4 = matrix1 * matrix2;
    // matrix4.print();


    // blas 3d test
    // Matrix id3(m);
    // id3.print();
    // Matrix m1 = Matrix(m,m);
    // m1.assign_random();
    // m1.print();

    // Matrix m2 = BLAS_3D(1, m1, id3);
    // m2.print();
    return 0;
}
