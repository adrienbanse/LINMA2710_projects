#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "EulerDiffusionSolver.hpp"
#include "DiffusionProblem.hpp"
#include "globals.hpp"

using namespace std;

void write_tensor(const char fname[], Vector** A, int nt, int nx, int ny) {
    FILE* f = fopen(fname, "w");
    if (f == NULL) {
        perror("fopen");
    }
    int err = fprintf(f, "%d %d %d\n", nt+1, nx+1, ny+1);
    if (err < 0) fprintf(stderr, "Error when writing. Do you have enough space left ?");
    for (int k = 1; k <= nt+1; k++)
        for (int i = 1; i <= nx+1; i++)
            for (int j = 1; j <= ny+1; j++) {
                    err = fprintf(f, "%lf\n", (*A[k-1]).Read((j-1)*(nx+1)+i-1));
                    if (err < 0) fprintf(stderr, "Error when writing. Do you have enough space left ?");
                }
    err = fclose(f);
    if (err != 0) perror("fclose");
}

int main() {
    // Taken from: 
    // https://levelup.gitconnected.com/solving-2d-heat-equation-numerically-using-python-3334004aa01a
    double Lx = 50., Ly = 50.;
    double Lt = 300.;
    double alpha = 2.;

    int i = 3; int nt = 2400 * pow(2., i);
    int j = 4; int nx = 8 * pow(2., j);
    int k = 3; int ny = 8 * pow(2., k);

    double dx = Lx / nx;
    double dy = Ly / ny;
    double dt = Lt / nt;

    int N = (nx+1) * (ny+1);

    Vector *sol[nt+1];
    for (int i = 0; i <= nt; i++)
        sol[i] = new Vector(N);

    // Your problem here (we use lambda function but not mandatory)
    std::function<double(double,double)> getD = [](double, double) {return 2.;};
    std::function<double(double,double)> getInitial = [](double, double) {return 0.;};
    std::function<double(double,double)> getf = [](double, double) {return 0.;};

    std::function<double(double)> getBoundaryx0 = [](double) {return 100.;};
    std::function<double(double)> getBoundaryx1 = [](double) {return 0.;};
    std::function<double(double)> getBoundaryy0 = [](double) {return 0.;};
    std::function<double(double)> getBoundaryy1 = [](double) {return 0.;};

    DiffusionProblem p(Lt, Lx, Ly, getD, getInitial, getf, getBoundaryx0, getBoundaryx1, getBoundaryy0, getBoundaryy1);
    
    EulerDiffusionSolver solver(nt, nx, ny, p);
    solver.solve(sol);

    char filename[1024];
    sprintf(filename, "sol_%d_%d_%d.txt", nx, ny, nt);

    //write_tensor(filename, sol, nt, nx, ny);
    write_tensor(filename, sol, nt, nx, ny);
    for (int i = 0; i <= nt; i++) {
        delete sol[i];
    }

    return EXIT_SUCCESS;
}
