#include <iostream>
#include <vector>
#include <matrix.cpp>
using namespace std;

/**
 * @brief BLAS3 (a) subroutine C = alpha*A*B + beta*C 
 * 
 * @param alpha 
 * @param beta 
 * @param m 
 * @param n 
 * @param p 
 * @param A 
 * @param B 
 * @param C 
 */
void BLAS_3A(int alpha, int beta, int m, int n, int p,
            vector<vector <int>> A,
            vector<vector <int>> B,
            vector<vector <int>> &C);



/**
 * @brief BLAS3 (d) subroutine B = alpha*(T^-1)*B
 * 
 * @param alpha 
 * @param m 
 * @param p 
 * @param T 
 * @param B 
 */
void BLAS_3D(int alpha, int m, int p, vector<vector <int> > T, vector<vector <int> > &B);