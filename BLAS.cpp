/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * BLAS.cpp - implementation of BLAS3 operation (a) and (d)
**/

#include <iostream>
#include <vector>
#include <matrix.cpp>
using namespace std;


/** BLAS level3 function a): C = alpha*A*B + beta*C
    ---- parameters ----
    [in] alpha : scalar
    [in] beta  : scalar
    [in] m : int
    [in] n : int
    [in] p : int
    [in] A : m x n matrix
    [in] B : n x p matrix
    [in,out] C : m x p matrix
**/
void BLAS_3A(int alpha, int beta, int m, int n, int p,
            vector<vector <int>> A,
            vector<vector <int>> B,
            vector<vector <int>> &C){
    if (m < p){
        // A = alpha * A
    }
    else{
        // B = alpha * B
    }
    // Stressen algorith A*B

    // + beta * C
}


void BLAS_3D(int alpha, int m, int p, vector<vector <int> > T, vector<vector <int> > &B);