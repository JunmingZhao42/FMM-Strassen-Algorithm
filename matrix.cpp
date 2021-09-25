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
 * @brief Get the column object
 * 
 * @param j 
 * @return double* 
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
 * @brief 
 * 
 * @param vector 
 * @return Vector& 
 */
Vector& Vector::operator=(const Vector& vector) 
{
    if (this == &vector){
        return *this;
    }

    // if dimension different, destroy original data
    if (n_len != vector.n_len){
        delete [] data;
        n_len = vector.n_len;
        alloc_space();
    }

    // copy data
    for (int i=0; i<n_len; i++){
        data[i] = vector.data[i];
    }

    return *this;
}


/**
 * @brief Vector addition
 * 
 * @return Vector& 
 */
Vector& Vector::operator+=(const Vector& vector) 
{
    if (n_len != vector.n_len){
        throw std::invalid_argument("dimensions not matching for vector addition");
    }

    for (int i=0; i<n_len; i++){
        data[i] += vector.data[i];
    }
    return *this;
}


/**
 * @brief Vector subtraction
 * 
 * @return Vector& 
 */
Vector& Vector::operator-=(const Vector& vector) 
{
    if (n_len != vector.n_len){
        throw std::invalid_argument("dimensions not matching for vector subraction");
    }

    for (int i=0; i<n_len; i++){
        data[i] -= vector.data[i];
    }
    return *this;
}


/**
 * @brief scalar vector multipliation
 * 
 * @return Vector& 
 */
Vector& Vector::operator*=(double alpha) 
{
    for (int i=0; i<n_len; i++){
        data[i] *= alpha;
    }
    return *this;
}


/**
 * @brief inner product of two vectors
 * 
 * @return double 
 */
double Vector::inner_product(const Vector& vector) {
    if (n_len != vector.n_len){
        throw std::invalid_argument("dimensions not matching for inner product");
    }

    double result = 0;
    for (int i=0; i<n_len; i++){
        result += data[i]*vector.data[i];
    }
    return result;
}


/**
 * @brief 2_norm of the vector
 * 
 * @return double 
 */
double Vector::norm2(){
    return inner_product(*this);
}


/**
 * @brief get subvector of the current vector
 * 
 * @param m 
 * @param m1 
 * @return Vector 
 */
Vector Vector::slice(int m, int m1){
    Vector subvector(m);
    return subvector;
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


double inner_product(int n, double* v1, double* v2){
    double s = 0;
    for (int i=0; i<n; i++){
        s += v1[i]*v2[i];
    }
    return s;
}


/**
vector : 1xm
matrix : mxn
result : 1xn
**/
double* vector_matrix_mul(double* v, Matrix A){
    double * result = new double[A.n_cols];
    for (int i=0; i<A.n_cols; i++){
        result[i] = inner_product(A.m_rows, v, A.get_column(i));
    }
    return result;
}


/**
vector : mxn
matrix : nx1
result : mx1
**/
double* matrix_vector_mul(Matrix A, double* v){
    double * result = new double[A.m_rows];
    for (int i=0; i<A.m_rows; i++){
        result[i] = inner_product(A.n_cols, A.data[i], v);
    }
    return result;
}


Matrix strassen(Matrix A, Matrix B){
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

    // 0. deal with odd dimension
    int m2 = m/2;
    int n2 = n/2;
    int p2 = p/2;

    Matrix C(m,p);

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
    

    // 5. deal with odd dimensions
    // TODO: pass C pointer
    if (m2*2 < m){
        // A has odd number of rows
        // fill last row of C
        C.data[m-1] = vector_matrix_mul(A.data[m-1], B.slice(n,0,2*p2,0));
    }
    if (p2*2 < p){
        double * s2 = matrix_vector_mul(A, B.get_column(p-1));
        // fill last col of C
        for (int i=0; i<m; i++){
            // std::cout << s2[i] << std::endl;
            C.data[i][p-1] = s2[i];
        }
        delete[] s2;
    }
    if (n2*2 < n){
        // n is odd
        for (int i=0; i<m2*2; i++){
            for (int j=0; j<p2*2; j++){
                // add entries to A1 x B1
                C.data[i][j] += A.data[i][n-1]*B.data[n-1][j];
            }
        }
    }

    return C;
}



