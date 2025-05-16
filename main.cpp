#include <iostream>
#include <cmath>
#include <cstring>
#include <pthread.h>
#include <chrono>
#include "matrix.h"
#include "approximation.h"
#include "solver.h"
#include "functions.h"

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " a b c d n_x n_y k ε m_i p" << std::endl;
    std::cerr << "  a, b, c, d: boundaries of area [a,b]×[c,d] (double)" << std::endl;
    std::cerr << "  n_x, n_y: number of interpolation points on X and Y axes (int)" << std::endl;
    std::cerr << "  k: function to approximate (int, 0-7)" << std::endl;
    std::cerr << "  ε: accuracy of the solution (double)" << std::endl;
    std::cerr << "  m_i: maximum number of iterations (int)" << std::endl;
    std::cerr << "  p: number of computational threads (int)" << std::endl;
}

int main(int argc, char* argv[]) {
    const int TASK_ID = 1; // Task identifier

    // Check command line arguments
    if (argc != 11) {
        printUsage(argv[0]);
        return 1;
    }

    // Parse command line arguments
    double a = std::atof(argv[1]);
    double b = std::atof(argv[2]);
    double c = std::atof(argv[3]);
    double d = std::atof(argv[4]);
    int n_x = std::atoi(argv[5]);
    int n_y = std::atoi(argv[6]);
    int k = std::atoi(argv[7]);
    double eps = std::atof(argv[8]);
    int m_i = std::atoi(argv[9]);
    int p = std::atoi(argv[10]);

    // Validate arguments
    if (a >= b || c >= d || n_x <= 0 || n_y <= 0 || k < 0 || k > 7 || eps <= 0 || m_i <= 0 || p <= 0) {
        std::cerr << "Invalid argument values!" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Calculate grid step sizes
    double h_x = (b - a) / n_x;
    double h_y = (d - c) / n_y;

    // Create context for the problem
    ApproximationContext context;
    context.a = a;
    context.b = b;
    context.c = c;
    context.d = d;
    context.n_x = n_x;
    context.n_y = n_y;
    context.h_x = h_x;
    context.h_y = h_y;
    context.k = k;
    context.eps = eps;
    context.m_i = m_i;
    context.p = p;

    // Timing variables
    double t1 = 0.0, t2 = 0.0;
    int it = 0;
    double r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 0.0;

    // Start timing for coefficient calculation
    auto start_t1 = std::chrono::high_resolution_clock::now();

    // Build the sparse matrix and solve the system
    SparseMatrix matrix;
    buildMatrixStructure(matrix, context);
    calculateGramMatrix(matrix, context);

    double* b_vector = new double[(n_x + 1) * (n_y + 1)];
    calculateRightHandSide(b_vector, context);

    double* solution = new double[(n_x + 1) * (n_y + 1)];
    it = solveSystem(matrix, b_vector, solution, context);

    // End timing for coefficient calculation
    auto end_t1 = std::chrono::high_resolution_clock::now();
    t1 = std::chrono::duration<double>(end_t1 - start_t1).count();

    // Start timing for error calculation
    auto start_t2 = std::chrono::high_resolution_clock::now();

    // Calculate error metrics
    r1 = calculateC1Error(solution, context);
    r2 = calculateL1Error(solution, context);
    r3 = calculateC2Error(solution, context);
    r4 = calculateL2Error(solution, context);

    // End timing for error calculation
    auto end_t2 = std::chrono::high_resolution_clock::now();
    t2 = std::chrono::duration<double>(end_t2 - start_t2).count();

    // Print results in the required format
    printf(
        "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        argv[0], TASK_ID, r1, r2, r3, r4, t1, t2, it, eps, k, n_x, n_y, p
    );

    // Cleanup
    delete[] b_vector;
    delete[] solution;
    freeMatrix(matrix);

    return 0;
}
