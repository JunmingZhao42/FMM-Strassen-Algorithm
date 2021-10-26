class Matrix{
    public:
        // constructors
        Matrix(int m, int n);
        Matrix(const Matrix&);
        Matrix(int m, int n, double** refd);
        Matrix(int m);

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
        inline int get_rows() {return m_rows;}
        inline int get_cols() {return n_cols;}
        inline double operator()(int i, int j) {return data[i][j];}

        void assign_random();
        void assign_triangular();
        void assign_zeros();
        void print();
        double* get_column(int j);

        // static helper functions
        static double inner_product(int n, double* v1, double* v2);
        static void vector_matrix_mul(double* v, Matrix A, double * & v_out);
        static double* matrix_vector_mul(Matrix A, double* v);
        static Matrix strassen(Matrix A, Matrix B);

    private:
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
