This is a one-semester project about implementing Strassen's algorithm in `C++` and exploiting it in BLAS3-A and BLAS3-D specificaiton.

--------
## Code structure
`test.cpp` includes three sets of test.
Set 1. Matrix multiplication AxB
Set 2. BLAS3-A routine using Strassen's algorithm
Set 3. BLAS3-D routine using Strassen's algorithm
Test will record the timing and generate `csv` file.

`plot.ipynb` provides visualisation from `csv` file.

## How to run
1. Compile
```
g++ test.cpp -std=c++11 -pthread -o test
```
2. Run
```
./test > result.csv
```

## References
Higham, N.J. (1990). Exploiting fast matrix multiplication within the level 3 BLAS. _ACM Transactions on Mathematical Software_, 16(4), pp.352–368.

Strassen, V. (1969). Gaussian elimination is not optimal. _Numerische Mathematik_, 13(4), pp.354–356.

Brent, R. (1970). _Algorithms for matrix multiplication_, Technical Report TR-CS-70-157, DCS, Stanford, 3+52 pp.

--------

**MATH3512 Matrix Computations** semester project
- Supervisor: Assoc Prof. Linda Stals
- Author: Junming Zhao
- The Australian National Unviersity
- Department of Mathematics