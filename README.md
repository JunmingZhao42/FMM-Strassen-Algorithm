# BLAS3 (a) (b)

This is a one-semester project about the BLAS3 (basic linear algebra subroutine level 3) operation (a) and (d).
The project focus on implementing and analysing the behaviour of operation (A) and operation (B) of BLAS3.

## (a)
1. Multiplication is done using Strassen algorithm
2. To handle odd dimension matrix, chopping method is used according to R. P. BRENT, _Algorithms for Matrix Multiplication_



Reference for the algorithm:
NICHOLAS J. HIGHAM. _Exploiting Fast Matrix Multiplication_


# TODO:
1. ~~structure the code with headers and proper comments~~ tidy up code comments
2. implement BLAS3D
3. ~~implement error handling~~ (done)
4. ~~generalise matrix entry type~~ (discard)
5. ~~optimise code memory allocation and complexity with pointers~~ (done)
6. create testing for correctness
7. create testing for timing
8. change matrix constructor (do not assign 0s)
