class Vector{
    public:
        // constructors
        Vector(int n);
        Vector(int n, double * refd);

        // destructor
        ~Vector();

        // standard vector operations
        Vector& operator=(const Vector&);
        Vector& operator+=(const Vector&);
        Vector& operator-=(const Vector&);
        Vector& operator*=(double);

        // other vector opersations
        double inner_product(const Vector&);
        double norm2();
        Vector slice(int m, int m1);

        // helper function
        inline double operator()(int i) {return data[i];}
        void assign_random();
        void print();
        void alloc_space();
        int n_len;
        double * data;
};


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
        Vector get_row(int i);
        Vector get_column(int j);

        int m_rows;
        int n_cols;
        double ** data;
        void alloc_space();
        bool is_submatrix = false;
};






// some static methods for Matrix Class
Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Matrix strassen(Matrix A, Matrix B);


// methods for vector and matrix
Matrix operator*(const Matrix&, const Vector&);
Matrix operator*(const Vector&, const Matrix&);