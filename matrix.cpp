// MATH3512 Matrix Computations, The Australian National University
// Supervisor: Professor Linda Stals
// Student : u6633756 Junming Zhao

// matrix.cpp - implementation of matrix helper functions (including strassen algorithm)
#include <iostream>
#include <vector>
using namespace std;


/** ------ Definition Section ------ **/

/** ------ Implementation Section ------ **/
// TODO:
// 1. generalise the function input as type T
// 2. implement error handling
// 3. optimize the code



/** Generate a random matrix
    ---- input ----
    m : # of rows
    n : # of cols

    --- output ---
    A : m x n matrix to print
**/
//template<typename T>
vector<vector <int> > generate(int m, int n){

    vector<vector <int> > Matrix_Out;
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
void Matrix_Print(int m, int n, vector<vector <int> > A) {
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
vector<vector <int> > Matrix_Multiply(int m, int n, int p,
                                    vector<vector <int> > A,
                                    vector<vector <int> > B) {
    vector<vector <int> > Matrix_Out;
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
vector<vector <int> > Matrix_Add(int m, int n,
                                vector<vector <int> > A,
                                vector<vector <int> > B) {
    vector<vector <int> > Matrix_Out;
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
vector<vector <int> > Matrix_Sub(int m, int n,
                                vector<vector <int> > A,
                                vector<vector <int> > B) {
    vector<vector <int> > Matrix_Out;
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




/** Compute inner porduct of 2 length-n vectors
    ---- input ----
    n : int
    v1 : vector
    v2 : vector

    ---- output ----
    s = v1^T x v2
**/
int Inner_Product(int n, vector <int> v1, vector <int> v2){
    int s = 0;
    for (int i=0; i<n; i++){
        s += v1[i]*v2[i];
    }
    return s;
}


/** Temperary solution to get a column from a matrix
    TODO: replace this inefficient method!
**/
vector<int> Get_Col(int m, int c, vector<vector <int> > A){
    vector<int> result;
    result.resize(m);
    for (int i=0; i<m; i++){
        result[i] = A[i][c];
    }
    return result;
}


/**
vector : 1xm
matrix : mxn
result : 1xn
**/
vector<int> Vector_Matrix_Mul(int m, int n, vector <int> v, vector<vector <int> > A){
    vector <int> result;
    result.resize(n);
    for (int i=0; i<n; i++){
        result[i] = Inner_Product(m, v, Get_Col(m, i, A));
    }
    return result;
}


/**
vector : mxn
matrix : nx1
result : mx1
**/
vector<int> Matrix_Vector_Mul(int m, int n, vector<vector <int> > A, vector <int> v){
    vector <int> result;
    result.resize(m);
    for (int i=0; i<m; i++){
        result[i] = Inner_Product(n, A[i], v);
    }
    return result;
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
vector<vector <int> > Strassen(int m, int n, int p,
              vector<vector <int> > A,
              vector<vector <int> > B) {

    // if demension small enough, apply standard matrix multiplication
    // TODO : empirically test out this threshold
    if(m <= 2 || n <=2 || p <= 2) {
        cout << "base case of recursion" << endl;
        return Matrix_Multiply(m, n, p, A, B);
    }
    // otherwise apply Strassen algorithm recursively
    else {
        // 0. deal with odd dimension
        int m2 = m/2;
        int n2 = n/2;
        int p2 = p/2;

        vector <int> s1 (2*p2, 0);
        vector <int> s2 (m, 0);

        if (m2 * 2 < m){
            // cout << "m is odd" << endl;
            // if row is odd, chop last row
            s1 = Vector_Matrix_Mul(n, 2*p2, A[m-1], Matrix_Slice(0,n,0,2*p2,B));
        }
        if (p2 * 2 < p){
            // cout << "p is odd" << endl;
            // if col is odd, chop last col
            s2 = Matrix_Vector_Mul(m, n, A, Get_Col(n,p-1,B));
        }


        // 1. Get submatrix
        // cout << "strassen method" << m << n << p << endl;
        vector<vector <int> > A11 = Matrix_Slice(0, m2, 0, n2, A);
        vector<vector <int> > A12 = Matrix_Slice(0, m2, n2, n2*2, A);
        vector<vector <int> > A21 = Matrix_Slice(m2, m2*2, 0, n2, A);
        vector<vector <int> > A22 = Matrix_Slice(m2, m2*2, n2, n2*2, A);

        vector<vector <int> > B11 = Matrix_Slice(0, n2, 0, p2*2, B);
        vector<vector <int> > B12 = Matrix_Slice(0, n2, p2, p*2, B);
        vector<vector <int> > B21 = Matrix_Slice(n2, n2*2, 0, p2, B);
        vector<vector <int> > B22 = Matrix_Slice(n2, n2*2, p2, p2*2, B);


        // 2. Calculate intermediate result
        // TODO: optimize the memory allocations

        // P1 = (A11 + A22) × (B11 + B22)
        vector<vector <int> > P1 = Strassen(m2, n2, p2, Matrix_Add(m2, n2, A11, A22), Matrix_Add(n2, p2, B11, B22));

        // P2 = (A21 + A22) × B11
        vector<vector <int> > P2 = Strassen(m2, n2, p2, Matrix_Add(m2, n2, A21, A22), B11);

        // P3 = A11 × (B12 - B22)
        vector<vector <int> > P3 = Strassen(m2, n2, p2, A11, Matrix_Sub(n2, p2, B12, B22));

        // P4 = A22 × (B21 - B22)
        vector<vector <int> > P4 = Strassen(m2, n2, p2, A22, Matrix_Sub(n2, p2, B21, B11));

        // P5 = (A11 + A12) × B22
        vector<vector <int> > P5 = Strassen(m2, n2, p2, Matrix_Add(m2, n2, A11, A12), B22);

        // P6 = (A21 - A11) × (B11 + B12)
        vector<vector <int> > P6 = Strassen(m2, n2, p2, Matrix_Sub(m2, n2, A21, A11), Matrix_Add(n2, p2, B11, B12));

        // P7 = (A12 - A22) × (B21 + B22)
        vector<vector <int> > P7 = Strassen(m2, n2, p2, Matrix_Sub(m2, n2, A12, A22), Matrix_Add(n2, p2, B21, B22));


        // 3. Combine intermediates to get result C
        // C11 = P1 + P4 - P5 + P7
        vector<vector <int> > C11 = Matrix_Add(m2, p2, Matrix_Sub(m2, p2, Matrix_Add(m2, p2, P1, P4), P5), P7);

        // C12 = P3 + P5
        vector<vector <int> > C12 = Matrix_Add(m2, p2, P3, P5);

        // C21 = P2 + P4
        vector<vector <int> > C21 = Matrix_Add(m2, p2, P2, P4);

        // C22 = P1 - P2 + P3 + P6
        vector<vector <int> > C22 = Matrix_Add(m2, p2, Matrix_Add(m2, p2, Matrix_Sub(m2, p2, P1, P2), P3), P6);


        // 4. Copy C11, C12, C21, C22 to C
        // cout << "compose together to get C" << endl;
        // TODO: do the previous with pointer so no need to copy again
        vector<vector <int> > C (m, vector<int> (p, 0));

        for(int i=0; i<m2; i++) {
            for(int j=0; j<p2; j++) {
                C[i][j]       = C11[i][j];
                C[i][j+p2]    = C12[i][j];
                C[i+m2][j]    = C21[i][j];
                C[i+m2][j+p2] = C22[i][j];
           }
        }
        Matrix_Print(m,p,C);

        // 5. deal with odd dimensions
        // cout << "deal with odd dimensions" << endl;
        // if A odd cols (B odd rows)
        if (n2*2 < n){
            // cout << "n is odd" << endl;
            for (int i=0; i<m2*2; i++){
                for (int j=0; j<p2*2; j++){
                    // add entries to A1 x B1
                    // cout << C[i][j] << endl;
                    // cout << A[i][n-1]*B[n-1][j] << endl;
                    C[i][j] += A[i][n-1]*B[n-1][j];
                }
            }
        }

        // if A odd rows
        if (m2*2 < m){
            // cout << "m is odd" << endl;
            // add last row of C (except last entry)
            for (int j=0; j<p2*2; j++){
                C[m-1][j] = s1[j];
            }
        }

        // if B odd cols
        if (p2*2 < p){
            // cout << "p is odd" << endl;
            // add last col of C
            for (int i=0; i<m; i++){
                // cout << s2[i] << endl;
                C[i][p-1] = s2[i];
            }
        }
        // Matrix_Print(m,p,C);
        return C;
    }
}


int main() {
    // generate matrices A and B
    int m = 5;
    int n = 3;
    int p = 5;

    vector <vector <int> > A = generate(m, n);
    vector <vector <int> > B = generate(n, p);

    // cout << "Print A :" << endl;
    // Matrix_Print(m,n,A);

    // cout << "Print B :" << endl;
    // Matrix_Print(n,p,B);

    // Strassen algorithm: C=A*B
    vector <vector <int> > C = Strassen(m, n, p, A, B);

    // print result C
    cout << "Print result of C=A*B strassen :" << endl;
    Matrix_Print(m, p, C);

    cout << "Print result of C=A*B standard :" << endl;
    vector <vector <int>> D = Matrix_Multiply(m, n, p, A, B);
    Matrix_Print(m, p, D);

    return 0;
}



/*
(potential) TODO : create a self-defined matrix class similar to below
reference : https://stackoverflow.com/questions/15778377/get-the-first-column-of-a-matrix-represented-by-a-vector-of-vectors
template <class T>
class SimpleMatrix
{
public:
    SimpleMatrix( int rows, int cols, const T& initVal = T() );

    // Size and structure
    int NumRows() const                       { return m_rows; }
    int NumColumns() const                    { return m_cols; }
    int NumElements() const                   { return m_data.size(); }

    // Direct vector access and indexing
    operator const vector<T>& () const        { return m_data; }
    int Index( int row, int col ) const       { return row * m_cols + col; }

    // Get a single value
          T & Value( int row, int col )       { return m_data[Index(row,col)]; }
    const T & Value( int row, int col ) const { return m_data[Index(row,col)]; }
          T & operator[]( size_t idx )        { return m_data[idx]; }
    const T & operator[]( size_t idx ) const  { return m_data[idx]; }

    // Simple row or column slices
    vector<T> Row( int row, int colBegin = 0, int colEnd = -1 ) const;
    vector<T> Column( int row, int colBegin = 0, int colEnd = -1 ) const;

private:
    vector<T> StridedSlice( int start, int length, int stride ) const;

    int m_rows;
    int m_cols;

    vector<T> m_data;
};

template <class T>
vector<T> SimpleMatrix<T>::StridedSlice( int start, int length, int stride ) const
{
    vector<T> result;
    result.reserve( length );
    const T *pos = &m_data[start];
    for( int i = 0; i < length; i++ ) {
        result.push_back(*pos);
        pos += stride;
    }
    return result;
}

template <class T>
SimpleMatrix<T>::SimpleMatrix( int rows, int cols, const T& initVal )
    : m_data( rows * cols, initVal )
    , m_rows( rows )
    , m_cols( cols )
{
}

template <class T>
vector<T> SimpleMatrix<T>::Row( int row, int colBegin, int colEnd ) const
{
    if( colEnd < 0 ) colEnd = m_cols-1;
    if( colBegin <= colEnd )
        return StridedSlice( Index(row,colBegin), colEnd-colBegin+1, 1 );
    else
        return StridedSlice( Index(row,colBegin), colBegin-colEnd+1, -1 );
}

template <class T>
vector<T> SimpleMatrix<T>::Column( int col, int rowBegin, int rowEnd ) const
{
    if( rowEnd < 0 ) rowEnd = m_rows-1;
    if( rowBegin <= rowEnd )
        return StridedSlice( Index(rowBegin,col), rowEnd-rowBegin+1, m_cols );
    else
        return StridedSlice( Index(rowBegin,col), rowBegin-rowEnd+1, -m_cols );
}
*/
