#include <math.h>

double boxcar_overlap(double rA, double radiiA, double rB, double radiiB)
{

    double F;
    F = pow(rB - rA, 2) / pow(radiiA + radiiB, 2);

    return F;

}