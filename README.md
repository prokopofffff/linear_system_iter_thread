# Finite Element Approximation with Sparse Matrices

This program implements a method for constructing function approximation by finite elements of degree 1 using the least squares method on parallel computers with shared memory.

## Method Details
- **Iterative Method**: Minimum error method with preconditioner
- **Preconditioner**: Jacobi preconditioner
- **Domain**: Rectangular region [a,b]×[c,d]
- **Approximation**: Continuous linear functions on triangles using the least squares method

## Building the Program

To build the program, run:

```bash
make
```

This will create an executable called `fem_approx`.

## Running the Program

The program takes 10 mandatory command-line arguments:

1. `a` - left boundary of the area (double)
2. `b` - right boundary of the area (double)
3. `c` - bottom boundary of the area (double)
4. `d` - top boundary of the area (double)
5. `n_x` - number of interpolation points on X axis (int)
6. `n_y` - number of interpolation points on Y axis (int)
7. `k` - function to be approximated (int, 0-7)
8. `ε` - accuracy of the solution (double)
9. `m_i` - maximum number of iterations (int)
10. `p` - number of computational threads (int)

### Example Usage

```bash
./fem_approx 0 1 0 1 10 10 2 1e-6 1000 4
```

This will approximate function k=2 (f(x,y) = y) on the unit square [0,1]×[0,1] with a 10×10 grid, accuracy of 1e-6, maximum 1000 iterations, and using 4 computational threads.

### Available Test Functions

The program supports 8 different functions to approximate:

1. For k = 0: f(x,y) = 1
2. For k = 1: f(x,y) = x
3. For k = 2: f(x,y) = y
4. For k = 3: f(x,y) = x + y
5. For k = 4: f(x,y) = sqrt(x^2 + y^2)
6. For k = 5: f(x,y) = x^2 + y^2
7. For k = 6: f(x,y) = e^(x^2 - y^2)
8. For k = 7: f(x,y) = 1/(25(x^2 + y^2) + 1)

### Output

The program outputs results in the following format:

```
<program_name> : Task = <task_id> R1 = <r1> R2 = <r2> R3 = <r3> R4 = <r4> T1 = <t1> T2 = <t2> It = <it> E = <eps> K = <k> Nx = <nx> Ny = <ny> P = <p>
```

Where:
- `r1` - max error at triangle centroids (C1 norm)
- `r2` - sum of weighted errors at triangle centroids (L1 norm)
- `r3` - max error at grid points (C2 norm)
- `r4` - sum of weighted errors at grid points (L2 norm)
- `t1` - time to build and solve the system (seconds)
- `t2` - time to calculate errors (seconds)
- `it` - number of iterations performed by the solver
- Other parameters match the input arguments

## Implementation Details

The program uses:
- Modified Sparse Row (MSR) format for sparse matrices
- pthreads for parallel computation
- Finite elements of degree 1 (linear functions on triangles)
- Least squares method for approximation
- Minimum error method with Jacobi preconditioner for solving the system
