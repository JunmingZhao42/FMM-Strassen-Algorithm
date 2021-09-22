Report
1. simple algebraic proof about the correctness of the algorithm
2. total flops (with dimension = 2^r for some positive integer r)
3. total flops (general case)
4. empirical testing result
5. **memory analysis?**



Some Questions:
1. `vector` object in C++ is not the most efficient way. Should use `[][]` and pointer instead?
2. Memory handling for the intermediate matrices. Have a global variable of fix size instead of creating new one each time?
3. Should we use generic matrix entry type in this project?
