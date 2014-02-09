#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "matrix_operations.c"
#include "ellipsoid/parametric_function_roots.c"
#include "polynomial.c"





void characteristic_ellipsoid_matrix(double *X, double *R, double phi, double *rotax, double exponent)
{

    double w, x, y, z;
    w = cos(phi / 2.0);
    x = sin(phi / 2.0) * rotax[0];
    y = sin(phi / 2.0) * rotax[1];
    z = sin(phi / 2.0) * rotax[2];

    // Rotation matrix
    double Q[3][3] = {{ 0 }};

    Q[0][0] = 1.0 - 2.0 * (pow(y, 2) + pow(z, 2));
    Q[0][1] = 2.0 * (x * y + w * z);
    Q[0][2] = 2.0 * (x * z - w * y);

    Q[1][0] = 2.0 * (x * y - w * z);
    Q[1][1] = 1.0 - 2.0 * (pow(x, 2) + pow(z, 2));
    Q[1][2] = 2.0 * (y * z + w * x);

    Q[2][0] = 2.0 * (x * z + w * y);
    Q[2][1] = 2.0 * (y * z - w * x);
    Q[2][2] = 1.0 - 2.0 * (pow(x, 2) + pow(y, 2));


    //Rotation matrix transpose
    double Qt[3][3] = {{ 0 }};
    matrix_transpose(&Qt[0][0], &Q[0][0], 3, 3);


    // radii matrix
    double O[3][3] = {{ 0 }};
    double diag_vals[3];
    diag_vals[0] = pow(R[0], exponent * -2.0);
    diag_vals[1] = pow(R[1], exponent * -2.0);
    diag_vals[2] = pow(R[2], exponent * -2.0);
    set_diagonal(&O[0][0], diag_vals, 3, 3);


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










double container_cube_overlap_potential(double *rA, double *radiiA, double phiA, double *rotaxA)
{

    double rB[3];
    double radiiB[3];
    double phiB;
    double rotaxB[3];

    double top, bottom, left, right, front, back;

    // top
    rB[0] = 0.5;
    rB[1] = 0.5;
    rB[2] = 2;

    radiiB[0] = INFINITY;
    radiiB[1] = INFINITY;
    radiiB[2] = 1;

    rotaxB[0] = 1;
    rotaxB[1] = 0;
    rotaxB[2] = 0;

    phiB = 0.0;

    top = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0]);


    // bottom
    rB[0] = 0.5;
    rB[1] = 0.5;
    rB[2] = -1;

    radiiB[0] = INFINITY;
    radiiB[1] = INFINITY;
    radiiB[2] = 1.0;

    rotaxB[0] = 1;
    rotaxB[1] = 0;
    rotaxB[2] = 0;

    phiB = 0;

    bottom = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0]);



    // left
    rB[0] = 0.5;
    rB[1] = -1.0;
    rB[2] = 0.5;

    radiiB[0] = INFINITY;
    radiiB[1] = 1;
    radiiB[2] = INFINITY;

    rotaxB[0] = 1;
    rotaxB[1] = 0;
    rotaxB[2] = 0;

    phiB = 0;

    left = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0]);



    // right
    rB[0] = 0.5;
    rB[1] = 2;
    rB[2] = 0.5;

    radiiB[0] = INFINITY;
    radiiB[1] = 1;
    radiiB[2] = INFINITY;

    rotaxB[0] = 1;
    rotaxB[1] = 0;
    rotaxB[2] = 0;

    phiB = 0;

    right = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0]);


    // front
    rB[0] = -1.0;
    rB[1] = 0.5;
    rB[2] = 0.5;

    radiiB[0] = 1.0;
    radiiB[1] = INFINITY;
    radiiB[2] = INFINITY;

    rotaxB[0] = 1.0;
    rotaxB[1] = 0.0;
    rotaxB[2] = 0.0;

    phiB = 0.0;

    front = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0]);




    // back
    rB[0] = 2.0;
    rB[1] = 0.5;
    rB[2] = 0.5;

    radiiB[0] = 1;
    radiiB[1] = INFINITY;
    radiiB[2] = INFINITY;

    rotaxB[0] = 1.0;
    rotaxB[1] = 0.0;
    rotaxB[2] = 0.0;

    phiB = 0;

    back = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0]);

    return fminf(top, fminf(bottom, fminf(left, fminf(right, fminf(front, back)))));

}









void random_point_on_sphere(double *rotax, double theta, double z)
{

    double x, y;

    x = sqrt(1 - pow(z, 2)) * cos(theta);
    y = sqrt(1 - pow(z, 2)) * sin(theta);

    rotax[0] = x;
    rotax[1] = y;
    rotax[2] = z;

}






size_t rsa_align_cube(double *x, double *y, double *z,
    size_t npoints, double *radius, double *rotax, double phi,
    int step_limit, unsigned long randSeed)
{

    // Setup GSL random number generator
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // Set the seed
    // srand ( time(NULL) );
    // unsigned long randSeed = rand();
    gsl_rng_set(r, randSeed);

    double rA[3];
    double radiiA[3];
    double phiA;
    double rotaxA[3];

    double rB[3];
    double radiiB[3];
    double phiB;
    double rotaxB[3];

    size_t valid_pts;
    double F, G;
    int k, flag, step;

    step = 0;
    valid_pts = 1;



    // Get new ellipsoid position and orientation
    G = 0;
    while (G < 1)
    {

        rB[0] = gsl_rng_uniform (r);
        rB[1] = gsl_rng_uniform (r);
        rB[2] = gsl_rng_uniform (r);
        radiiB[0] = radius[0];
        radiiB[1] = radius[1];
        radiiB[2] = radius[2];
        phiB = phi;
        rotaxB[0] = rotax[0];
        rotaxB[1] = rotax[1];
        rotaxB[2] = rotax[2];

        G = container_cube_overlap_potential(&rB[0], &radiiB[0], phiB, &rotaxB[0]);

    }

    // store ellipsoid position
    x[0] = rB[0];
    y[0] = rB[1];
    z[0] = rB[2];





    // Generate ellipsoid positions
    while ((valid_pts < npoints) & (step < step_limit))
    {

        // Get new ellipsoid position and orientation
        G = 0;
        while (G < 1)
        {

            rB[0] = gsl_rng_uniform (r);
            rB[1] = gsl_rng_uniform (r);
            rB[2] = gsl_rng_uniform (r);
            radiiB[0] = radius[0];
            radiiB[1] = radius[1];
            radiiB[2] = radius[2];
            phiB = phi;
            rotaxB[0] = rotax[0];
            rotaxB[1] = rotax[1];
            rotaxB[2] = rotax[2];

            G = container_cube_overlap_potential(&rB[0], &radiiB[0], phiB, &rotaxB[0]);

        }



        flag = 1;
        for (k = 0; k < valid_pts; k++)
        {

            rA[0] = x[k];
            rA[1] = y[k];
            rA[2] = z[k];
            radiiA[0] = radius[0];
            radiiA[1] = radius[1];
            radiiA[2] = radius[2];
            phiA = phi;
            rotaxA[0] = rotax[0];
            rotaxA[1] = rotax[1];
            rotaxA[2] = rotax[2];


            F = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0],
                &rB[0], &radiiB[0], phiB, &rotaxB[0]);

            if (F < 1.0)
            {

                flag = 0;
                break;

            }
        }
        if (flag == 1)
        {

           x[valid_pts] = rB[0];
           y[valid_pts] = rB[1];
           z[valid_pts] = rB[2];
           valid_pts += 1;

        }

        step += 1;

    }


    gsl_rng_free (r);

    return valid_pts;

}