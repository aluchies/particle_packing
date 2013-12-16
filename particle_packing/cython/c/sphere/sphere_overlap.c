#include <math.h>

double sphere_overlap(double *rA, double radiiA, double *rB, double radiiB)
{

    double rAB[3];
    rAB[0] = rB[0] - rA[0];
    rAB[1] = rB[1] - rA[1];
    rAB[2] = rB[2] - rA[2];

    double F;
    F = (pow(rAB[0], 2) + pow(rAB[1], 2) + pow(rAB[2], 2)) / pow(radiiA + radiiB, 2);

    return F;

}