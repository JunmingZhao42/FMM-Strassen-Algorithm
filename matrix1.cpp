#include "matrix1.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>


// ---- class functions implementation ----


// ---- constructors ----
/**
 * @brief Construct a new zero matrix
 * 
 * @param m number of rows
 * @param n number of cols
 */
Matrix::Matrix(int m, int n) : m_rows(m), n_cols(n)
{
    alloc_space();
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
Matrix::Matrix(const Matrix& matrix) : m_rows(matrix.m_rows), n_cols(matrix.n_cols) 
{
    alloc_space();
    for (int i=0; i<m_rows; i++){
        for (int j=0; j<n_cols; j++){
            data[i][j] = matrix.data[i][j];
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
    for (int i=0; i<m_rows; i++){
        delete[] data[i];
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
        return;
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
        return;
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
        return;
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
 * @brief Get a submatrix = Matrix[m1:m2][n1:n2]
 * excluding row m2 and col n2
 * 
 * @param m1 
 * @param m2 
 * @param n1 
 * @param n2 
 * @return Matrix
 */
Matrix Matrix::slice(int m1, int m2, int n1, int n2) 
{
    Matrix result(m2-m1, n2-n1);
    for (int i=0; i<result.m_rows; i++){
        for (int j=0; j<result.n_cols; j++){
            result.data[i][j] = data[i+m1][j+m2];
        }
    }
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
            data[i][j] = rand() % 5;
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
