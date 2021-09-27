#include "matrix.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>


// ---- class functions implementation ----


// ---- constructors ----
/**
 * @brief Construct a new mxn matrix
 * 
 * @param m number of rows
 * @param n number of cols
 */
Matrix::Matrix(int m, int n) : 
m_rows(m), n_cols(n)
{
    alloc_space();
}


/**
 * @brief Construct a new matrix by copying 
 * 
 * @param matrix 
 */
Matrix::Matrix(const Matrix& matrix) : 
m_rows(matrix.m_rows), n_cols(matrix.n_cols) 
{
    alloc_space();
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] = matrix.data[i][j];
        }
    }
}


/**
 * @brief Construct a new Matrix:: Matrix object
 * data referes to the input refd
 * @param m number of rows
 * @param n number of cols
 * @param data pointer of double array data
 */
Matrix::Matrix(int m, int n, double** refd) : 
m_rows(m), n_cols(n), data(refd) {}



/**
 * @brief Construct a new Matrix:: Matrix object
 * matrix is identity matrix
 * 
 * @param m number of rows, number of cols
 */
Matrix::Matrix(int m) : m_rows(m), n_cols(m)
{
    alloc_space();
    for (int i=0; i<m; i++){
        for (int j=0; j<m; j++){
            if (i==j) data[i][j] = 1;
            else data[i][j] = 0;
            //data[i][j] = (i==j)? : 1, 0; 
        }
    }
}


// ---- destructor ----
/**
 * @brief Destroy the Matrix:: Matrix object
 * 
 */
Matrix::~Matrix() 
{
    if (!is_submatrix){
        for (int i=0; i<m_rows; i++){
            delete[] data[i];
        }
    }
    delete[] data;
}


// ---- standard matrix operations ----
/**
 * @brief assign matrix
 * 
 * @param matrix 
 * @return Matrix& 
 */
Matrix& Matrix::operator=(const Matrix& matrix) 
{   
    // if pointers are the same
    if (this == &matrix){
        return *this;
    }

    // if dimension different, destroy original data
    if ((m_rows != matrix.m_rows) || (n_cols != matrix.n_cols)){
        for (int i = 0; i < m_rows; ++i){
            delete[] data[i];
        }
        delete[] data;

        m_rows = matrix.m_rows;
        n_cols = matrix.n_cols;
        alloc_space();
    }

    // copy value from input matrix
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] = matrix.data[i][j];
        }
    }
    return *this;
}


/**
 * @brief Matrix addition
 * 
 * @param matrix 
 * @return Matrix& 
 */
Matrix& Matrix::operator+=(const Matrix& matrix) 
{
    if ((m_rows != matrix.m_rows) || (n_cols != matrix.n_cols)){
        throw std::invalid_argument("dimensions not matching for matrix addition");
    }

    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] += matrix.data[i][j];
        }
    }
    return *this;
}


/**
 * @brief Matrix subtraction
 * 
 * @param matrix 
 * @return Matrix& 
 */
Matrix& Matrix::operator-=(const Matrix& matrix) 
{
    if ((m_rows != matrix.m_rows) || (n_cols != matrix.n_cols)){
        throw std::invalid_argument("dimensions not matching for matrix subraction");
    }

    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] -= matrix.data[i][j];
        }
    }
    return *this;
}


/**
 * @brief Matrix multiplication
 * 
 * @param matrix 
 * @return Matrix& 
 */
Matrix& Matrix::operator*=(const Matrix& matrix) 
{
    if (n_cols != matrix.m_rows){
        throw std::invalid_argument("dimensions not matching for matrix multiplication");
    }

    Matrix result(m_rows, matrix.n_cols);
    result.assign_zeros();

    for (int i=0; i<result.m_rows; i++){
        for (int j=0; j<result.n_cols; j++){
            for (int k=0; k<n_cols; k++){
                result.data[i][j] += (data[i][k] * matrix.data[k][j]);
            }
        }
    }
    return (*this = result);
}


/**
 * @brief Matrix scalar multiplication
 * 
 * @return Matrix& 
 */
Matrix& Matrix::operator*=(double alpha) 
{
    for (int i=0; i<m_rows; i++) {
        for (int j=0; j<n_cols; j++) {
            data[i][j] *= alpha;
        }
    }
    return *this;
}


/**
 * @brief Matrix transpose
 * 
 * @return Matrix 
 */
Matrix Matrix::transpose() 
{
    Matrix result(n_cols, m_rows);
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            result.data[j][i] = data[i][j];
        }
    }
    return result;
}


/**
 * @brief Create submatrix = Matrix[m1:m1+m][n1:n1+n]
 * excluding row (m1+m) and col (n1+n)
 * 
 * @param m  number of rows
 * @param m1 row starting point
 * @param n  number of cols
 * @param n2 col starting point
 * @return Matrix with data pointer cooresponds to the parent matrix
 */
Matrix Matrix::slice(int m, int m1, int n, int n1) 
{
    if ((m+m1 > m_rows) || (n+n1 > n_cols)){
        throw std::invalid_argument("submatrix index out of the original one");
    }

    double ** sub_data = new double*[m];

    for(int i=0; i<m; i++){
        sub_data[i] = data[i+m1]+n1;
    }
    Matrix result(m, n, sub_data);
    result.is_submatrix = true;
    return result;
}


/**
 * @brief assign random value (0-4) to entries 
 * 
 */
void Matrix::assign_random() 
{
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] = (i+1)*(j+1);
        }
    }
}


/**
 * @brief assign zeros to all entries
 * 
 */
void Matrix::assign_zeros(){
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] = 0;
        }
    }
}


/**
 * @brief print out matrix data
 * 
 */
void Matrix::print() 
{
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


/**
 * @brief Get jth column of a matrix
 * 
 * @param j column index
 * @return double*, pointer of column array
 */
double* Matrix::get_column(int j){
    double* result = new double[m_rows];
    for (int i=0; i<m_rows; i++){
        result[i] = data[i][j];
    }
    return result;
}


/**
 * @brief allocate double array space by m_rows and n_cols
 * 
 */
void Matrix::alloc_space() 
{
    data = new double*[m_rows];
    for (int i=0; i < m_rows; i++){
        data[i] = new double[n_cols];
    }
}



// ---- public static function implementations ----
Matrix operator+(const Matrix& m1, const Matrix& m2) 
{
    Matrix result(m1);
    return (result += m2);
}


Matrix operator-(const Matrix& m1, const Matrix& m2) 
{
    Matrix result(m1);
    return (result -= m2);
}


Matrix operator*(const Matrix& m1, const Matrix& m2) 
{
    Matrix result(m1);
    return (result *= m2);
}


Matrix operator*(const Matrix& matrix, double alpha) 
{
    Matrix result(matrix);
    return (result *= alpha);
}


Matrix operator*(double alpha, const Matrix& matrix){
    return (matrix * alpha);
}


double Matrix::inner_product(int n, double* v1, double* v2){
    double s = 0;
    for (int i=0; i<n; i++){
        s += v1[i]*v2[i];
    }
    return s;
}


/**
 * @brief vector matrix multiplication
 * 
 * @param v 1xm vector
 * @param A mxn matrix
 * @param v_out 1xn row vector
 */
void Matrix::vector_matrix_mul(double* v, Matrix A, double * & v_out){
    for (int i=0; i<A.n_cols; i++){
        v_out[i] = inner_product(A.m_rows, v, A.get_column(i));
    }
}


/**
 * @brief matrix vector multiplication
 * 
 * @param A mxn matrix
 * @param v nx1 vector
 * @return double* mx1 col vector
 */
double* Matrix::matrix_vector_mul(Matrix A, double* v){
    double * result = new double[A.m_rows];
    for (int i=0; i<A.m_rows; i++){
        result[i] = inner_product(A.n_cols, A.data[i], v);
    }
    return result;
}


/**
 * @brief Strassen algorithm
 * 
 * @param A mxn matrix
 * @param B nxp matrix
 * @return  mxp matrix
 */
Matrix Matrix::strassen(Matrix A, Matrix B){
    int m = A.m_rows;
    int n = A.n_cols;
    int p = B.n_cols;

    if (n!= B.m_rows){
        std::invalid_argument("dimensions not matching for matrix multiplication");
    }

    if (m <= 1 || n <= 1 || p <= 1){
        // std::cout << "base case of recursion" << std::endl;
        return (A*B);
    }

    // 0. prepare
    int m2 = m/2;
    int n2 = n/2;
    int p2 = p/2;

    Matrix C(m,p);
    C.assign_zeros();

    // 1. get submatrix of A, B
    Matrix A11 = A.slice(m2,0,n2,0);
    Matrix A12 = A.slice(m2,0,n2,n2);
    Matrix A21 = A.slice(m2,m2,n2,0);
    Matrix A22 = A.slice(m2,m2,n2,n2);

    Matrix B11 = B.slice(n2,0,p2,0);
    Matrix B12 = B.slice(n2,0,p2,p2);
    Matrix B21 = B.slice(n2,n2,p2,0);
    Matrix B22 = B.slice(n2,n2,p2,p2);


    // 2. calculate intermediate result
    Matrix P1 = strassen((A11 + A22),(B11 + B22));
    Matrix P2 = strassen((A21 + A22), B11);
    Matrix P3 = strassen(A11, (B12 - B22));
    Matrix P4 = strassen(A22, (B21 - B11));
    Matrix P5 = strassen((A11 + A12), B22);
    Matrix P6 = strassen((A21 - A11), (B11 + B12));
    Matrix P7 = strassen((A12 - A22), (B21 + B22));

    // 3. compose intermediate result to get C
    Matrix C11 = C.slice(m2,0,p2,0);
    Matrix C12 = C.slice(m2,0,p2,p2);
    Matrix C21 = C.slice(m2,m2,p2,0);
    Matrix C22 = C.slice(m2,m2,p2,p2);

    // C11
    C11 += P1;
    C11 += P4;
    C11 -= P5;
    C11 += P7;

    // C12
    C12 += P3;
    C12 += P5;

    // C21
    C21 += P2;
    C21 += P4;

    // C22
    C22 += P1;
    C22 += P3;
    C22 -= P2;
    C22 += P6;
    

    // 4. deal with odd dimensions
    if (m2*2 < m){
        // m is odd
        // fill last row of C
        vector_matrix_mul(A.data[m-1], B.slice(n,0,2*p2,0), C.data[m-1]);
    }

    if (p2*2 < p){
        // p is odd
        // fill last col of C
        double * v = B.get_column(p-1);
        for (int i=0; i<m; i++){
            // std::cout << s2[i] << std::endl;
            C.data[i][p-1] = inner_product(A.n_cols, A.data[i], v);
        }
        delete[] v;
    }

    if (n2*2 < n){
        // n is odd
        // add entries to A1 x B1
        for (int i=0; i<m2*2; i++){
            for (int j=0; j<p2*2; j++){
                C.data[i][j] += A.data[i][n-1]*B.data[n-1][j];
            }
        }
    }

    return C;
}