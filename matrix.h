#include <iostream>
#include <vector>
using namespace std;


/**
 * @brief Generate a random matrix
 * 
 * @param m 
 * @param n 
 * @return vector<vector <int> > 
 */
vector<vector <int> > generate(int m, int n);


/**
 * @brief Print out a matrix
 * 
 * @param m 
 * @param n 
 * @param A 
 */
void Matrix_Print(int m, int n, vector<vector <int> > A);


/**
 * @brief Standard matrix multiplication: C = A x B
 * 
 * @param m 
 * @param n 
 * @param p 
 * @param A 
 * @param B 
 * @return vector<vector <int> > 
 */
vector<vector <int> > Matrix_Multiply(int m, int n, int p,
                                    vector<vector <int> > A,
                                    vector<vector <int> > B);


/**
 * @brief Standard Matrix Addition: C = A + B
 * 
 * @param m 
 * @param n 
 * @param A 
 * @param B 
 * @return vector<vector <int> > 
 */
vector<vector <int> > Matrix_Add(int m, int n,
                                vector<vector <int> > A,
                                vector<vector <int> > B);



/**
 * @brief Standard Matrix Subtraction : C = A - B
 * 
 * @param m 
 * @param n 
 * @param A 
 * @param B 
 * @return vector<vector <int> > 
 */
vector<vector <int> > Matrix_Sub(int m, int n,
                                vector<vector <int> > A,
                                vector<vector <int> > B);


/**
 * @brief Get submatrix
 * 
 * @param m1 
 * @param m2 
 * @param n1 
 * @param n2 
 * @param A 
 * @return vector<vector <int> > 
 */
vector<vector <int> > Matrix_Slice(int m1, int m2, int n1, int n2, vector<vector <int> > A);



/**
 * @brief Compute inner porduct of two length-n vectors
 * 
 * @param n 
 * @param v1 
 * @param v2 
 * @return int 
 */
int Inner_Product(int n, vector <int> v1, vector <int> v2);



/**
 * @brief Get c-th column from matrix A
 * TODO: replace this inefficient method!
 * 
 * @param m 
 * @param c 
 * @param A 
 * @return vector<int> 
 */
vector<int> Get_Col(int m, int c, vector<vector <int> > A);




/**
 * @brief [(1xn) row vector] = [(1xm) row vector] x [(mxn) matrix]
 * 
 * @param m 
 * @param n 
 * @param v 
 * @param A 
 * @return vector<int> 
 */
vector<int> Vector_Matrix_Mul(int m, int n, vector <int> v, vector<vector <int> > A);


/**
 * @brief [(mx1) col vector] = [(mxn) matrix] x [(nx1) row vector]
 * 
 * @param m 
 * @param n 
 * @param A 
 * @param v 
 * @return vector<int> 
 */
vector<int> Matrix_Vector_Mul(int m, int n, vector<vector <int> > A, vector <int> v);


/** Strassen algorithm for matrix multiplication : C = A x B
    ---- input ----
    m : int
    n : int
    p : int

    A : m x n matrix
    B : n x p matrix

    ---- output ----
    C : m x p matrix
**/
//template <typename T>

/**
 * @brief 
 * 
 * @param m 
 * @param n 
 * @param p 
 * @param A 
 * @param B 
 * @return vector<vector <int> > 
 */
vector<vector <int> > Strassen(int m, int n, int p,
              vector<vector <int> > A,
              vector<vector <int> > B);