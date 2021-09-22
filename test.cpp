/** 
 * MATH3512 Matrix Computations, The Australian National University
 * Supervisor : Professor Linda Stals
 * Student    : u6633756 Junming Zhao
 * test.cpp   : examples and tests of BLAS operation (a) and (d) implemented.
**/

#include <iostream>
using namespace std;

int main(){
    int m = 4;
    int n = 6;
    double **data;
    data = new double*[m];
    for (int i=0; i < m; i++){
        data[i] = new double[n];
    }
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            data[i][j] = i*j;
        }
    }

    for (unsigned int i=0; i<m; i++){
        for (unsigned int j=0; j<n; j++){
            cout << data[i][j] << " ";
        }
        cout << endl;
    }

    int m1 = 2;
    int m2 = 4;
    int n1 = 4;
    int n2 = 6;

    double **data1;
    data1 = new double*[m2-m1];

    for (int i=0; i < m; i++){
        data1[i] = new double[n2-n1];
    }

    for (int i=0; i<m2-m1; i++){
        for (int j=0; j<n2-n1; j++){
            cout << data[i+m1][j+m2] << endl;
            data1[i][j] = data[i+m1][j+m2];
        }
    }


    for (unsigned int i=0; i<m2-m1; i++){
        for (unsigned int j=0; j<n2-n1; j++){
            cout << data1[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
