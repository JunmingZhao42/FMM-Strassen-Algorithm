#include "matrix1.h"
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
        std::cout << "delete row" << std::endl;
        for (int i=0; i<m_rows; i++){
            delete[] data[i];
        }
    }
    std::cout << "delete data" << std::endl;
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
void Matrix::print_matrix() 
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

    int m2 = m/2;
    int n2 = n/2;
    int p2 = p/2;

    // assign temp variables
    Matrix P1(m2,p2);
    Matrix P2(m2,p2);
    Matrix P3(m2,p2);
    Matrix P4(m2,p2);
    Matrix P5(m2,p2);
    Matrix P6(m2,p2);
    Matrix P7(m2,p2);

    // get submatrix
    Matrix A11(m2,n2,A.data);
    Matrix A12(m2,n2,A.data);
    

    
}