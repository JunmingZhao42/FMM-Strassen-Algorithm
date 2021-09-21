#include <iostream>
#include <vector>
using namespace std;


vector<vector <int> > Matrix_Slice(int m1, int m2, int n1, int n2, vector<vector <int> > A){
    vector<vector <int> > Matrix_Out;

    for (int i=m1; i<m2; i++){
        Matrix_Out.emplace_back(A[i].begin()+n1, A[i].begin()+n2);
    }
    return Matrix_Out;
}


int main(){
    int m = 4;
    int n = 6;
    vector<vector <int> > A1 = {{1,2,3,4,5,6}, {4,5,6,7,8,9}, {1,3,4,6,8,0}, {0,5,7,8,1,2}};
    vector<vector <int> > A2 = Matrix_Slice(0, m, 0, n-1, A1);


    // A2.reserve(2);

    for (unsigned int i=0; i<A2.size(); i++){
        for (unsigned int j=0; j<A2[i].size(); j++){
            // cout << A1[i][j] << " ";
            cout << A2[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
