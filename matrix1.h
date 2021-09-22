class Matrix{
    public:
        // constructors
        Matrix(int, int);
        Matrix(const Matrix&);


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
        Matrix slice(int, int, int, int);
        

        // helper functions
        inline int get_m_rows() {return m_rows;}
        inline int get_n_cols() {return n_cols;}
        inline double operator()(int i, int j) {return data[i][j];}
        void assign_random();
        void print_matrix();

    private:
        int m_rows;
        int n_cols;
        double ** data;
        void alloc_space();
};



// some static methods for Matrix Class
Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
