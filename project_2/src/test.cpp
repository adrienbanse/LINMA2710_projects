#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "EulerDiffusionSolver.hpp"
#include "DiffusionProblem.hpp"
#include "globals.hpp"
#include "AbstractVector.hpp"

bool isapprox(const Vector& x, const Vector& y, double rtol=1e-8, double atol=1e-8) {
    int n = length(x);
    if (length(y) != n) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        double a = x.Read(i);
        double b = y.Read(i);
        if (!isapprox(a, b, rtol, atol)) {
            printf(ANSI_COLOR_RED " x[%d] = %lf != %lf = y[%d]\n" ANSI_COLOR_RESET, i+1, a, b, i+1);
            return false;
        }
    }
    return true;
}

void testisapprox(double x, double y, char const *s) {
    printf(ANSI_COLOR_BLUE "%s" ANSI_COLOR_RESET ":", s);
    if (!isapprox(x, y)) {
        printf(ANSI_COLOR_RED "%lf != %lf\n" ANSI_COLOR_RESET, x, y);
        exit(EXIT_FAILURE);
    } else {
        printf(ANSI_COLOR_GREEN " OK\n" ANSI_COLOR_RESET);
    }
}

void testisapprox(Vector const &x, Vector const &y, char const *s) {
    printf(ANSI_COLOR_BLUE "%s" ANSI_COLOR_RESET ":", s);
    if (!isapprox(x, y)) {
        exit(EXIT_FAILURE);
    } else {
        printf(ANSI_COLOR_GREEN " OK\n" ANSI_COLOR_RESET);
    }
}

int main() {
    double Lt = 2., Lx = 1., Ly = 1., C = 1.5, D = 0.0;
    int nt = 8, nx = 4, ny = 6;
    double dt = Lt / nt;
    int N = (nx+1) * (ny+1);

    // constant problem
    // expected solution
    Vector *exp[nt+1], *sol[nt+1];
    for (int i = 0; i <= nt; i++) {
        exp[i] = new Vector(N);
        sol[i] = new Vector(N);
        for (int j = 0; j < N; j++)
            (*(exp[i]))(j) = C;
    }

    //solve and check
    #define SLEN 16
    char s[SLEN];
    std::function<double(double)> constantC = [C](double) {return C;};
    DiffusionProblem p(Lt, Lx, Ly, [D](double,double) {return D;}, [C](double,double) {return C;}, [](double, double) {return 0.;}, constantC, constantC, constantC, constantC);
    EulerDiffusionSolver constsolver(nt, nx, ny, p);
    snprintf(s, SLEN, "csol");
    constsolver.solve(sol);
    for (int i = 0; i <= nt; i++) {
        snprintf(s, SLEN, "csol%2d", i);
        testisapprox(*sol[i], *exp[i], s);
    }

    //free ressources
    for (int i = 0; i <= nt; i++) {
        delete exp[i];
        delete sol[i];
    }

    return EXIT_SUCCESS;
}
