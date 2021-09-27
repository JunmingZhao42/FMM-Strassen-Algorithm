/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * BLAS.cpp - implementation of BLAS3 operation (a) and (d)
**/

#include "matrix.hpp"
#include <stdexcept>


/** BLAS level3 function a): C = alpha*A*B + beta*C
    ---- parameters ----
    [in] alpha : scalar
    [in] beta  : scalar
    [in] A : m x n matrix
    [in] B : n x p matrix
    [in,out] C : m x p matrix
**/
void BLAS_3A(double alpha, double beta,
            Matrix A, Matrix B, Matrix &C){
    C *= beta;
    if (A.m_rows < B.n_cols){
        // A = alpha * A
        C += strassen(alpha*A,B);
    }
    else{
        // B = alpha * B
        C += strassen(A, alpha*B);
    }
}


Matrix BLAS_3D(double alpha, Matrix T, Matrix B){
    int m = T.m_rows;
    int p = B.n_cols;

    if (m != T.n_cols){
        throw std::invalid_argument("T is not square matrix");
    } 
    if (m != B.m_rows){
        throw std::invalid_argument("dimensions not matching for BLAS3-d");
    }

    // base case
    if (m == 1){
        return (1/(T.data[0][0]))*B;
    }

    Matrix C(m,p);
    C.assign_zeros();

    int m2 = m/2;
    int p2 = p/2;

    Matrix T11 = T.slice(m2,0,m2,0);
    Matrix T12 = T.slice(m2,0,m2,m2);
    Matrix T22 = T.slice(m2,m2,m2,m2);

    Matrix B11 = B.slice(m2,0,p2,0);
    Matrix B12 = B.slice(m2,0,p2,p2);
    Matrix B21 = B.slice(m2,m2,p2,0);
    Matrix B22 = B.slice(m2,m2,p2,p2);

    Matrix C11 = C.slice(m2,0,p2,0);
    Matrix C12 = C.slice(m2,0,p2,p2);
    Matrix C21 = C.slice(m2,m2,p2,0);
    Matrix C22 = C.slice(m2,m2,p2,p2);

    C21 += BLAS_3D(1, T22, B21);
    C22 += BLAS_3D(1, T22, B22);
    C11 += BLAS_3D(1, T11, (B11 - strassen(T12, C21)));
    C12 += BLAS_3D(1, T11, (B12 - strassen(T12, C22)));

    return C;
}