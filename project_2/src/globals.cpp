#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "globals.hpp"

double frand() {
    // srand needs to be called only once
    static bool firstcall = true;
    if (firstcall) {
        srand(time(NULL));
        firstcall = false;
    }
    // rand() is in the range [0, RAND_MAX]
    // the output of this function is in the range [-1, 1]
    return 2 * (rand() / RAND_MAX) - 1;
}

// rtol: relative tolerance
// atol: absolute tolerance
bool isapprox (double a, double b, double rtol, double atol) {
    return fabs(a - b) < atol + rtol*fmax(fabs(a), fabs(b));
}
