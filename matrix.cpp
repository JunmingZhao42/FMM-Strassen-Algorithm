// MATH3512 Matrix Computations, The Australian National University
// Supervisor: Professor Linda Stals
// Student : u6633756 Junming Zhao

// matrix.cpp - implementation of matrix helper functions (including strassen algorithm)

#include <iostream>
#include <vector>
using namespace std;


/** ------ Definition Section ------ **/
// Size of the matrix
// #define N 8
const int N = 8;

// template<typename T>
// void Strassen(int n, T A[][N], T B[][N], T C[][N]);
//
// template<typename T>
// void output(int n, T C[][N]);



/** ------ Implementation Section ------ **/
// TODO: generalise the function input as type T


/** Generate a random matrix
    ---- input ----
    m : # of rows
    n : # of cols

    --- output ---
    A : m x n matrix to print
**/
//template<typename T>
vector<vector <int>> generate(int m, int n){

    vector<vector <int>> Matrix_Out;
    Matrix_Out.resize(m);

    for (int i=0; i<m; i++){
        Matrix_Out[i].resize(n);
        for (int j=0; j<n; j++){
            Matrix_Out[i][j] = rand() % 5;
        }
    }
    return Matrix_Out;
}


/** Print out a matrix
    ---- input ----
    m : # of rows
    n : # of cols
    A : m x n matrix to print
**/

//template<typename T>
void Matrix_Print(int m, int n, vector<vector <int>> A) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            cout << A[i][j];
            cout << " ";
        }
        cout << "\n";
    }
}


/** Standard matrix multiplication: C = A x B
    ---- input ----
    m : # of rows of matrix A
    n : # of cols of matrix A
        # of rows of matrix B
    p : # of cols of matrix B
    A : m x k matrix
    B : k x n matrix

    ---- output ----
    C : m x p matrix
**/
//template<typename T>
vector<vector <int>> Matrix_Multiply(int m, int n, int p,
                                    vector<vector <int>> A,
                                    vector<vector <int>> B) {
    vector<vector <int>> Matrix_Out;
    Matrix_Out.resize(m);

    for (int i=0; i<m; i++){
        Matrix_Out[i].resize(p);
        for (int j=0; j<p; j++){
            Matrix_Out[i][j] = 0;
            for (int k=0; k<n; k++){
                Matrix_Out[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return Matrix_Out;
}


/** Standard Matrix Addition: C = A + B
    ---- input ----
    m : # of rows of matrix A, B
    n : # of cols of matrix A, B
    A : m x n matrix
    B : m x n matrix

    ---- output ----
    [in, out] C : m x n matrix
**/
//template <typename T>
vector<vector <int>> Matrix_Add(int m, int n,
                                vector<vector <int>> A,
                                vector<vector <int>> B) {
    vector<vector <int>> Matrix_Out;
    Matrix_Out.resize(m);

    for (int i=0; i<m; i++) {
        Matrix_Out[i].resize(n);
        for (int j=0; j<n; j++) {
           Matrix_Out[i][j] = A[i][j] + B[i][j];
        }
    }
    return Matrix_Out;
}

/** Standard Matrix Subtraction : C = A - B
    ---- input ----
    m : # of rows of matrix A, B
    n : # of cols of matrix A, B
    A : m x n matrix
    B : m x n matrix

    --- output ----
    C : m x n matrix
**/
//template <typename T>
vector<vector <int>> Matrix_Sub(int m, int n,
                                vector<vector <int>> A,
                                vector<vector <int>> B) {
    vector<vector <int>> Matrix_Out;
    Matrix_Out.resize(m);

    for (int i=0; i<m; i++) {
        Matrix_Out[i].resize(n);
        for(int j=0; j<n; j++) {
           Matrix_Out[i][j] = A[i][j] - B[i][j];
       }
    }
    return Matrix_Out;
}



/** Get submatrix
    ---- input ----
    m1 : start row index of the submatrix
    m2 : end row index (exclude) of the submatrix
    n1 : start column index of the submatrix
    n2 : end column index (exclude) of the submatrix
    A  : the original (parent) matrix

    ---- output ----
    Submatrix = A[m1:m2, n1:n2]
**/
vector<vector <int> > Matrix_Slice(int m1, int m2, int n1, int n2, vector<vector <int> > A){
    vector<vector <int> > Matrix_Out;

    for (int i=m1; i<m2; i++){
        Matrix_Out.emplace_back(A[i].begin()+n1, A[i].begin()+n2);
    }
    return Matrix_Out;
}


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
vector<vector <int>> Strassen(int m, int n, int p,
              vector<vector <int>> A,
              vector<vector <int>> B) {

    // if demension small enough, apply standard matrix multiplication
    // TODO : empirically test out this threshold
    if(n <= 2) {
        return Matrix_Multiply(m, p, n, A, B);
    }
    // otherwise apply Strassen algorithm recursively
    else {
        // 1. deal with odd dimension
        // if ()
        // 1. Get submatrix
        vector<vector <int>> A11 = Matrix_Slice(0, m/2, 0, n/2, A);
        vector<vector <int>> A12 = Matrix_Slice(0, m/2, n/2, n, A);
        vector<vector <int>> A21 = Matrix_Slice(m/2, m, 0, n/2, A);
        vector<vector <int>> A22 = Matrix_Slice(m/2, m, n/2, n, A);

        vector<vector <int>> B11 = Matrix_Slice(0, n/2, 0, p/2, B);
        vector<vector <int>> B12 = Matrix_Slice(0, n/2, p/2, p, B);
        vector<vector <int>> B21 = Matrix_Slice(n/2, n, 0, p/2, B);
        vector<vector <int>> B22 = Matrix_Slice(n/2, n, p/2, p, B);


        // 2. Calculate intermediate result
        // TODO: optimize the memory allocations

        // P1 = (A11 + A22) × (B11 + B22)
        vector<vector <int>> P1 = Strassen(m/2, n/2, p/2, Matrix_Add(m/2, n/2, A11, A22), Matrix_Add(n/2, p/2, B11, B22));

        // P2 = (A21 + A22) × B11
        vector<vector <int>> P2 = Strassen(n/2, n/2, n/2, Matrix_Add(m/2, n/2, A21, A22), B11);

        // P3 = A11 × (B12 - B22)
        vector<vector <int>> P3 = Strassen(n/2, n/2, n/2, A11, Matrix_Sub(n/2, n/2, B12, B22));

        // P4 = A22 × (B21 - B22)
        vector<vector <int>> P4 = Strassen(n/2, n/2, n/2, A22, Matrix_Sub(n/2, n/2, B21, B11));

        // P5 = (A11 + A12) × B22
        vector<vector <int>> P5 = Strassen(n/2, n/2, n/2, Matrix_Add(n/2, n/2, A11, A12), B22);

        // P6 = (A21 - A11) × (B11 + B12)
        vector<vector <int>> P6 = Strassen(n/2, n/2, n/2, Matrix_Sub(n/2, n/2, A21, A11), Matrix_Add(n/2, n/2, B11, B12));

        // P7 = (A12 - A22) × (B21 + B22)
        vector<vector <int>> P7 = Strassen(n/2, n/2, n/2, Matrix_Sub(n/2, n/2, A12, A22), Matrix_Add(n/2, n/2, B21, B22));


        // 3. Combine intermediates to get result C
        // C11 = P1 + P4 - P5 + P7
        vector<vector <int>> C11 = Matrix_Add(n/2, n/2, Matrix_Sub(n/2, n/2, Matrix_Add(n/2, n/2, P1, P4), P5), P7);

        // C12 = P3 + P5
        vector<vector <int>> C12 = Matrix_Add(n/2, N/2, P3, P5);

        // C21 = P2 + P4
        vector<vector <int>> C21 = Matrix_Add(n/2, n/2, P2, P4);

        // C22 = P1 - P2 + P3 + P6
        vector<vector <int>> C22 = Matrix_Add(n/2, n/2, P1, Matrix_Add(n/2, n/2, Matrix_Sub(n/2, n/2, P3, P2), P6));


        // 4. Copy C11, C12, C21, C22 to C
        // TODO: do the previous with pointer so no need to copy again
        vector<vector <int>> C;
        C.resize(p);
        for(int i=0; i<n/2; i++) {
            C[i].resize(m);
            C[i+n/2].resize(m);
            for(int j=0; j<n/2; j++) {
                C[i][j]         = C11[i][j];
                C[i][j+n/2]     = C12[i][j];
                C[i+n/2][j]     = C21[i][j];
                C[i+n/2][j+n/2] = C22[i][j];
           }
        }
        return C;
    }
}


int main() {
    // generate matrices A and B
    vector <vector <int>> A = generate(2, 3);
    vector <vector <int>> B = generate(3, 4);

    cout << "Print A :" << endl;
    Matrix_Print(2,3,A);

    cout << "Print B :" << endl;
    Matrix_Print(3,4,B);

    // Strassen algorithm: C=A*B
    // vector <vector <int>> C = Strassen(N, N, N, A, B);

    // print result C
    // cout << "Print result of C=A*B strassen :" << endl;
    // Matrix_Print(N, N, C);

    cout << "Print result of C=A*B standard :" << endl;
    vector <vector <int>> D = Matrix_Multiply(2, 3, 4, A, B);
    Matrix_Print(2, 4, D);

    return 0;
}
