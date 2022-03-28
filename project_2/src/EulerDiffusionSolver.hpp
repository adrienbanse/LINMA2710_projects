#ifndef EULERDIFFUSIONSOLVERHEADERDEF
#define EULERDIFFUSIONSOLVERHEADERDEF
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "DiffusionProblem.hpp"

// This file is a proposition, feel free to change it, but the constructor and "solve(sol)" should still be implemented

class EulerDiffusionSolver
{
    private:
        int nt;
        int nx;
        int ny;
        double dt;
        double dx;
        double dy;
        double alphax; // dt / (dx * dx)
        double alphay; // dt / (dy * dy)
        const DiffusionProblem &p;
        SparseMatrix *A;
        Vector *b;
        SparseMatrix* computeAb(Vector &b);// return A and change b in argument
    public:
        // construct a solver with parameters nx, ny, nt to solve a DiffusionProblem, precompute A and b
        EulerDiffusionSolver(int nt, int nx, int ny, const DiffusionProblem &p);

        // solve the problem and set the solution at each time step in sol (vectorized)
        void solve(Vector **sol);
};
#endif
