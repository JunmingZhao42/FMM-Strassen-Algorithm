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

    // TODO: this doesn't have to assign 0s
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            data[i][j] = 0;
        }
    }

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
 * @param m 
 * @param n 
 * @param data 
 */
Matrix::Matrix(int m, int n, double** refd) : 
m_rows(m), n_cols(n), data(refd) {}


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
 * @brief Matrix addition: add input matrix into data
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
 * @brief Matrix subtraction: subtract input matrix from data
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
 * @brief matrix scalar multiplication
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
 * @brief return a transposed new matrix
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
 * @brief Get a pointer submatrix = Matrix[m1:m1+m][n1:n1+n]
 * excluding row (m1+m) and col (n1+n)
 * 
 * @param m  number of rows
 * @param m1 row starting point
 * @param n  number of cols
 * @param n2 col starting point
 * @return Matrix
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
 * @brief assign random value (0-5) to entries 
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



// ---- implemenation of class Vector  ----

/**
 * @brief Construct a new Vector:: Vector object
 * 
 * @param n 
 */
Vector::Vector(int n) : n_len(n)
{
    alloc_space();
    for (int i=0; i<n; i++){
        data[i] = 0;
    }
}


/**
 * @brief Construct a new Vector:: Vector object
 * 
 * @param n 
 * @param refd 
 */
Vector::Vector(int n, double * refd) : 
n_len(n), data(refd)
{}


/**
 * @brief Destroy the Vector:: Vector object
 * 
 */
Vector::~Vector() 
{
    delete[] data;
}


/**
 * @brief assign random value to the vector
 * 
 */
void Vector::assign_random() 
{
    for (int i=0; i<n_len; i++){
        data[i] = (i+1);
    }
}


/**
 * @brief print out the vector horizontally
 * 
 */
void Vector::print() 
{
    for (int i=0; i<n_len; i++){
       std::cout << data[i] << " ";
    }
    std::cout << std::endl;
}


/**
 * @brief allocate data space to the pointer
 * 
 */
void Vector::alloc_space() 
{
    data = new double[n_len];
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



Matrix strassen(Matrix A, Matrix B){
    int m = A.m_rows;
    int n = A.n_cols;
    int p = B.n_cols;

    if (n!= B.m_rows){
        std::invalid_argument("dimensions not matching for matrix multiplication");
    }

    if (m <= 2 || n <= 2 || p <= 2){
        // std::cout << "base case of recursion" << std::endl;
        return (A*B);
    }

    int m2 = m/2;
    int n2 = n/2;
    int p2 = p/2;

    // assign temp variables

    // get submatrix
    Matrix A11 = A.slice(m2,0,n2,0);
    Matrix A12 = A.slice(m2,0,n2,n2);
    Matrix A21 = A.slice(m2,m2,n2,0);
    Matrix A22 = A.slice(m2,m2,n2,n2);

    Matrix B11 = B.slice(n2,0,p2,0);
    Matrix B12 = B.slice(n2,0,p2,p2);
    Matrix B21 = B.slice(n2,n2,p2,0);
    Matrix B22 = B.slice(n2,n2,p2,p2);


    // calculate intermediate result
    Matrix P1 = strassen((A11 + A22),(B11 + B22));
    Matrix P2 = strassen((A21 + A22), B11);
    Matrix P3 = strassen(A11, (B12 - B22));
    Matrix P4 = strassen(A22, (B21 - B11));
    Matrix P5 = strassen((A11 + A12), B22);
    Matrix P6 = strassen((A21 - A11), (B11 + B12));
    Matrix P7 = strassen((A12 - A22), (B21 + B22));

    // compse to get result C
    Matrix C(m,p);
    Matrix C11 = C.slice(m2,0,p2,0);
    Matrix C12 = C.slice(m2,0,p2,p2);
    Matrix C21 = C.slice(m2,m2,p2,0);
    Matrix C22 = C.slice(m2,m2,p2,p2);

    C11 += P1;
    C11 += P4;
    C11 -= P5;
    C11 += P7;

    C12 += P3;
    C12 += P5;

    C21 += P2;
    C21 += P4;

    C22 += P1;
    C22 += P3;
    C22 -= P2;
    C22 += P6;
    

    return C;
}