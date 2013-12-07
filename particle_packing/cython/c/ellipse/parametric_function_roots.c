#include "parametric_function_roots.h"
#include <math.h>

double f(double L, double a11, double a12, double a21, double a22, double b1, double b2) {

   return (pow(L, 3)*(a11*pow(b2, 2) - a12*b1*b2 - a21*b1*b2 + a22*pow(b1, 2) - pow(b1, 2) - pow(b2, 2)) + pow(L, 2)*(-2*a11*pow(b2, 2) + 2*a12*b1*b2 + 2*a21*b1*b2 - 2*a22*pow(b1, 2) + pow(b1, 2) + pow(b2, 2)) + L*(a11*pow(b2, 2) - a12*b1*b2 - a21*b1*b2 + a22*pow(b1, 2)))/(pow(L, 2)*(a11*a22 - a11 - a12*a21 - a22 + 1) + L*(-2*a11*a22 + a11 + 2*a12*a21 + a22) + a11*a22 - a12*a21);

}

double h0(double a11, double a12, double a21, double a22, double b1, double b2) {

   return pow(a11, 2)*a22*pow(b2, 2) - a11*a12*a21*pow(b2, 2) - a11*a12*a22*b1*b2 - a11*a21*a22*b1*b2 + a11*pow(a22, 2)*pow(b1, 2) + pow(a12, 2)*a21*b1*b2 + a12*pow(a21, 2)*b1*b2 - a12*a21*a22*pow(b1, 2);

}

double h1(double a11, double a12, double a21, double a22, double b1, double b2) {

   return -4*pow(a11, 2)*a22*pow(b2, 2) + 4*a11*a12*a21*pow(b2, 2) + 4*a11*a12*a22*b1*b2 + 4*a11*a21*a22*b1*b2 - 4*a11*pow(a22, 2)*pow(b1, 2) + 2*a11*a22*pow(b1, 2) + 2*a11*a22*pow(b2, 2) - 4*pow(a12, 2)*a21*b1*b2 - 4*a12*pow(a21, 2)*b1*b2 + 4*a12*a21*a22*pow(b1, 2) - 2*a12*a21*pow(b1, 2) - 2*a12*a21*pow(b2, 2);

}

double h2(double a11, double a12, double a21, double a22, double b1, double b2) {

   return 6*pow(a11, 2)*a22*pow(b2, 2) - pow(a11, 2)*pow(b2, 2) - 6*a11*a12*a21*pow(b2, 2) - 6*a11*a12*a22*b1*b2 + a11*a12*b1*b2 - 6*a11*a21*a22*b1*b2 + a11*a21*b1*b2 + 6*a11*pow(a22, 2)*pow(b1, 2) - 6*a11*a22*pow(b1, 2) - 6*a11*a22*pow(b2, 2) + a11*pow(b1, 2) + 6*pow(a12, 2)*a21*b1*b2 + 6*a12*pow(a21, 2)*b1*b2 - 6*a12*a21*a22*pow(b1, 2) + 5*a12*a21*pow(b1, 2) + 5*a12*a21*pow(b2, 2) + a12*a22*b1*b2 + a12*b1*b2 + a21*a22*b1*b2 + a21*b1*b2 - pow(a22, 2)*pow(b1, 2) + a22*pow(b2, 2);

}

double h3(double a11, double a12, double a21, double a22, double b1, double b2) {

   return -4*pow(a11, 2)*a22*pow(b2, 2) + 2*pow(a11, 2)*pow(b2, 2) + 4*a11*a12*a21*pow(b2, 2) + 4*a11*a12*a22*b1*b2 - 2*a11*a12*b1*b2 + 4*a11*a21*a22*b1*b2 - 2*a11*a21*b1*b2 - 4*a11*pow(a22, 2)*pow(b1, 2) + 6*a11*a22*pow(b1, 2) + 6*a11*a22*pow(b2, 2) - 2*a11*pow(b1, 2) - 2*a11*pow(b2, 2) - 4*pow(a12, 2)*a21*b1*b2 - 4*a12*pow(a21, 2)*b1*b2 + 4*a12*a21*a22*pow(b1, 2) - 4*a12*a21*pow(b1, 2) - 4*a12*a21*pow(b2, 2) - 2*a12*a22*b1*b2 - 2*a21*a22*b1*b2 + 2*pow(a22, 2)*pow(b1, 2) - 2*a22*pow(b1, 2) - 2*a22*pow(b2, 2);

}

double h4(double a11, double a12, double a21, double a22, double b1, double b2) {

   return pow(a11, 2)*a22*pow(b2, 2) - pow(a11, 2)*pow(b2, 2) - a11*a12*a21*pow(b2, 2) - a11*a12*a22*b1*b2 + a11*a12*b1*b2 - a11*a21*a22*b1*b2 + a11*a21*b1*b2 + a11*pow(a22, 2)*pow(b1, 2) - 2*a11*a22*pow(b1, 2) - 2*a11*a22*pow(b2, 2) + a11*pow(b1, 2) + 2*a11*pow(b2, 2) + pow(a12, 2)*a21*b1*b2 + a12*pow(a21, 2)*b1*b2 - a12*a21*a22*pow(b1, 2) + a12*a21*pow(b1, 2) + a12*a21*pow(b2, 2) + a12*a22*b1*b2 - a12*b1*b2 + a21*a22*b1*b2 - a21*b1*b2 - pow(a22, 2)*pow(b1, 2) + 2*a22*pow(b1, 2) + a22*pow(b2, 2) - pow(b1, 2) - pow(b2, 2);

}
