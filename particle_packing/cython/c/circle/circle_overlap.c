#include <math.h>

double circle_overlap(double *rA, double radiiA, double *rB, double radiiB)
{

    double rAB[2];
    rAB[0] = rB[0] - rA[0];
    rAB[1] = rB[1] - rA[1];

    double F;
    F = (pow(rAB[0], 2) + pow(rAB[1], 2)) / pow(radiiA + radiiB, 2);

    return F;

}