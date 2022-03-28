#ifndef DIFFUSIONPROBLEMHEADERDEF
#define DIFFUSIONPROBLEMHEADERDEF

#include <functional>

/*
 * Diffusion problem
 * rho(x, y) * cv(x, y) * du/dt = k(x, y) (d^2u/dx^2 + d^2u/dy^2) + f(t, x, y)
 * From t = 0 to t = Lt
 * in the domain [0, Lx] x [0, Ly]
 */
class DiffusionProblem
{
    private:
        double Lt;
        double Lx;
        double Ly;

    public:
        DiffusionProblem(
            double Lt, 
            double Lx, 
            double Ly, 
            std::function<double(double,double)> getD, 
            std::function<double(double,double)> getInitial, 
            std::function<double(double,double)> getf, 
            std::function<double(double)> getBoundaryx0, 
            std::function<double(double)> getBoundaryx1, 
            std::function<double(double)> getBoundaryy0, 
            std::function<double(double)> getBoundaryy1
        );
        double getLt() const;
        double getLx() const;
        double getLy() const;

        std::function<double(double,double)> getD; // double getD(double x, double y) -> Value of D(x, y)
        std::function<double(double,double)> getInitial; // double getInitial(double x, double y) -> Value of u(0, x, y)
        std::function<double(double,double)> getf; //  double getf(double x, double y) -> Value of f(x, y)
        std::function<double(double)> getBoundaryx0; // double getBoundaryx0(double y) -> Value at y of the boundary condition at x=0
        std::function<double(double)> getBoundaryx1; // double getBoundaryx1(double y) -> Value at y of the boundary condition at x=Lx
        std::function<double(double)> getBoundaryy0; // double getBoundaryy0(double x) -> Value at x of the boundary condition at y=0
        std::function<double(double)> getBoundaryy1; // double getBoundaryy1(double x) -> Value at x of the boundary condition at y=Ly
};
#endif