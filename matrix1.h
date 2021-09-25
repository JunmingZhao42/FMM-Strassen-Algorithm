
class Matrix{
    public:
        // constructors
        Matrix(int m, int n);
        Matrix(const Matrix&);
        Matrix(int m, int n, double** refd);

        // destructor
        ~Matrix();


        // standard matrix operations
        Matrix& operator=(const Matrix&);
        Matrix& operator+=(const Matrix&);
        Matrix& operator-=(const Matrix&);
        Matrix& operator*=(const Matrix&);
        Matrix& operator*=(double);


        // other matrix operations
        Matrix transpose();
        Matrix slice(int m, int m1, int n, int n1);
        

        // helper functions
        // inline int get_m_rows() {return m_rows;}
        // inline int get_n_cols() {return n_cols;}
        inline double operator()(int i, int j) {return data[i][j];}
        void assign_random();
        void print();

        int m_rows;
        int n_cols;
        double ** data;
        void alloc_space();
        bool is_submatrix = false;
};


class Vector{
    public:
        // constructors
        Vector(int n);
        Vector(int n, double * refd);

        // destructor
        ~Vector();

        // helper function
        inline double operator()(int i) {return data[i];}
        void assign_random();
        void print();
        void alloc_space();
        int n_len;
        double * data;
};



// some static methods for Matrix Class
Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Matrix strassen(Matrix A, Matrix B);
