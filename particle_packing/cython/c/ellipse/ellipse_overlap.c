#include <math.h>
#include <stdio.h>
#include "../matrix_operations.c"
#include "parametric_function_roots.c"
#include "../polynomial.c"





void characteristic_ellipse_matrix(double *X, double *R, double phi, double exponent)
{


    // rotation matrix
    double Q[2][2] = {{ 0 }};
    double Qt[2][2] = {{ 0 }};
    Q[0][0] = cos(phi);
    Q[0][1] = sin(phi);
    Q[1][0] = -sin(phi);
    Q[1][1] = cos(phi);
    matrix_transpose(&Qt[0][0], &Q[0][0], 2, 2);


    // radii matrix
    double O[2][2] = {{ 0 }};
    R[0] = pow(R[0], exponent * -2.0);
    R[1] = pow(R[1], exponent * -2.0);
    set_diagonal(&O[0][0], R, 2, 2);


    // characteristic ellipse matrix
    double X_temp[2][2] = {{ 0 }};
    matrix_multiply(&X_temp[0][0], &O[0][0], &Q[0][0], 2, 2, 2);
    matrix_multiply(X, &Qt[0][0], &X_temp[0][0], 2, 2, 2);

}





double ellipse_overlap(double *rA, double *radiiA, double phiA, double *rB, double *radiiB, double phiB)
{


    // find XA^(-1) and XB^(1/2)
    double XA[2][2] = {{ 0 }};
    double XB[2][2] = {{ 0 }};
    characteristic_ellipse_matrix(&XA[0][0], &radiiA[0], phiA, -1.0);
    characteristic_ellipse_matrix(&XB[0][0], &radiiB[0], phiB, 0.5);


    // find A_AB
    double A_AB[2][2] = {{ 0 }};
    double A_temp[2][2] = {{ 0 }};
    matrix_multiply(&A_temp[0][0], &XA[0][0], &XB[0][0], 2, 2, 2);
    matrix_multiply(&A_AB[0][0], &XB[0][0], &A_temp[0][0], 2, 2, 2);

    // find r_AB
    double rAB[2];
    rAB[0] = rB[0] - rA[0];
    rAB[1] = rB[1] - rA[1];

    // find a_AB
    double a_AB[2];
    matrix_multiply(&a_AB[0], &XB[0][0], &rAB[0], 2, 2, 1);

    // extract elements of the matrix A_AB and a_AB
    double a11, a12, a21, a22;
    double b1, b2;
    a11 = A_AB[0][0];
    a12 = A_AB[0][1];
    a21 = A_AB[1][0];
    a22 = A_AB[1][1];
    b1 = a_AB[0];
    b2 = a_AB[1];




    // find coefficients for the parametric polynomial derivative used to find max
    double h[5];
    double z[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    size_t n = 5;
    h[0] = h0(a11, a12, a21, a22, b1, b2);
    h[1] = h1(a11, a12, a21, a22, b1, b2);
    h[2] = h2(a11, a12, a21, a22, b1, b2);
    h[3] = h3(a11, a12, a21, a22, b1, b2);
    h[4] = h4(a11, a12, a21, a22, b1, b2);



    // find roots
    size_t m;
    m = find_roots(&z[0], &h[0], n);



    double F;
    int i;

    if (m > 1)
    {
        if (f(0, a11, a12, a21, a22, b1, b2) > f(1, a11, a12, a21, a22, b1, b2))
        {

           F = f(0, a11, a12, a21, a22, b1, b2);

        }
        else
        {

           F = f(1, a11, a12, a21, a22, b1, b2);

        }


        int i;
        for(i=0; i<=m-1; i++)
        {

           if (( z[2 * i + 1] == 0 ) & (z[2 * i] > 0) & (z[2 * i] < 1))
           {

              F = f(z[2 * i], a11, a12, a21, a22, b1, b2);

           }

        }
    }
    else
    {
        F = 0.;
    }

    return F;



}