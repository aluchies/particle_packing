#include <math.h>
#include <stdio.h>
#include "../matrix_operations.c"
#include "parametric_function_roots.c"
#include "../polynomial.c"





void characteristic_ellipsoid_matrix(double *X, double *R, double phi, double *rotax, double exponent)
{


    // rotation matrix
    double I[3][3] = {{ 0.0 }};
    double v[3] = { 1.0, 1.0, 1.0 };
    set_diagonal(&I[0][0], &v[0], 3, 3);

    double s = cos(phi);

    double p[3][1] = {{ 0.0 }};
    p[0][0] = sin(phi) * rotax[0];
    p[1][0] = sin(phi) * rotax[1];
    p[2][0] = sin(phi) * rotax[2];

    double pT[1][3] = {{ 0 }};
    matrix_transpose(&pT[0][0], &p[0][0], 3, 1);

    double P[3][3] = {{ 0 }};
    matrix_multiply(&P[0][0], &p[0][0], &pT[0][0], 3, 1, 3);


    double s_temp1 = 1 - s;
    multiply_matrix_by_scalar(&P[3][3], s_temp1, 3, 3);


    double s_temp2 = pow(s, 2) - 0.5;
    multiply_matrix_by_scalar(&I[3][3], s_temp2, 3, 3);


    double Q[3][3] = {{ 0 }};
    double Qt[3][3] = {{ 0 }};
    add_matrices(&Q[0][0], &P[0][0], &I[0][0], 3, 3);
    matrix_transpose(&Qt[0][0], &Q[0][0], 3, 3);


    // radii matrix
    double O[3][3] = {{ 0 }};
    R[0] = pow(R[0], exponent * -2.0);
    R[1] = pow(R[1], exponent * -2.0);
    R[2] = pow(R[2], exponent * -2.0);
    set_diagonal(&O[0][0], R, 3, 3);


    // characteristic ellipse matrix
    double X_temp[3][3] = {{ 0 }};
    matrix_multiply(&X_temp[0][0], &O[0][0], &Q[0][0], 3, 3, 3);
    matrix_multiply(X, &Qt[0][0], &X_temp[0][0], 3, 3, 3);


}





double ellipsoid_overlap(double *rA, double *radiiA, double phiA, double *rotaxA, 
    double *rB, double *radiiB, double phiB, double *rotaxB)
{


    // find XA^(-1) and XB^(1/2)
    double XA[3][3] = {{ 0 }};
    double XB[3][3] = {{ 0 }};
    characteristic_ellipsoid_matrix(&XA[0][0], &radiiA[0], phiA, &rotaxA[0], -1.0);
    characteristic_ellipsoid_matrix(&XB[0][0], &radiiB[0], phiB, &rotaxB[0], 0.5);


    // find A_AB
    double A_AB[3][3] = {{ 0 }};
    double A_temp[3][3] = {{ 0 }};
    matrix_multiply(&A_temp[0][0], &XA[0][0], &XB[0][0], 3, 3, 3);
    matrix_multiply(&A_AB[0][0], &XB[0][0], &A_temp[0][0], 3, 3, 3);

    // find r_AB
    double rAB[3];
    rAB[0] = rB[0] - rA[0];
    rAB[1] = rB[1] - rA[1];
    rAB[2] = rB[2] - rA[2];

    // find a_AB
    double a_AB[3];
    matrix_multiply(&a_AB[0], &XB[0][0], &rAB[0], 3, 3, 1);

    // extract elements of the matrix A_AB and a_AB
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
    double b1, b2, b3;
    a11 = A_AB[0][0];
    a12 = A_AB[0][1];
    a13 = A_AB[0][2];
    a21 = A_AB[1][0];
    a22 = A_AB[1][1];
    a23 = A_AB[1][2];
    a31 = A_AB[2][0];
    a32 = A_AB[2][1];
    a33 = A_AB[2][2];
    b1 = a_AB[0];
    b2 = a_AB[1];
    b3 = a_AB[2];






    // find coefficients for the parametric polynomial derivative used to find max
    double h[7];
    double z[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    size_t n = 7;
    h[0] = h0(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);
    h[1] = h1(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);
    h[2] = h2(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);
    h[3] = h3(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);
    h[4] = h4(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);
    h[5] = h5(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);
    h[6] = h6(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);





    // find roots
    size_t m;
    m = find_roots(&z[0], &h[0], n);





    double F;
    int i;

    if (m > 1)
    {
        if (f(0, a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3) > f(1, a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3))
        {

           F = f(0, a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);

        }
        else
        {

           F = f(1, a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);

        }



        for(i=0; i<=m-2; i++)
        {

           if (( z[2 * i + 1] == 0 ) & (z[2 * i] > 0) & (z[2 * i] < 1))
           {

              F = f(z[2 * i], a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3);

           }

        }
    }
    else
    {
        F = 0.;
    }

    return F;



}