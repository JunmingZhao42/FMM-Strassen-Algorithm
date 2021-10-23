/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * test.cpp   : examples and tests of BLAS operation (a) and (d) implemented.
**/


#include "blas_routine.cpp"
// #include "matrix.cpp"
#include <iostream>
#include <chrono>


const int REPEAT = 5;
const int SIZE = 1000;
const int STEP = 10;


int main(){
    // // test for Strassen algorithm
    // std::cout << "matrix dimension," << "Strassen," << "Standard" << std::endl;
    // for (int m=0; m<SIZE; m+=STEP){
    //     for (int i=0; i<REPEAT; i++){
    //         std::cout << m << ",";
    //         // create matrix A
    //         Matrix A = Matrix(m);
    //         A.assign_random();

    //         // create matrix B
    //         Matrix B = Matrix(m);
    //         B.assign_random();

    //         std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //         // C = A*B
    //         Matrix C = Matrix::strassen(A,B);
    //         std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    //         std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << ",";

    //         // standard matrix multiplication
    //         begin = std::chrono::steady_clock::now();
    //         C = A*B;
    //         end = std::chrono::steady_clock::now();

    //         std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
    //     }
    // }
    

    // test for BLAS3-A routine
    std::cout << "matrix dimension," << "Strassen," << "Standard" << std::endl;
    
    for (int m=0; m < SIZE; m+=STEP){
        for (int i=0; i<REPEAT; i++){
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
    }

    return 0;
}
