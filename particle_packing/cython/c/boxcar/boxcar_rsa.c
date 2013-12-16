#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <math.h>

// #include <stdlib.h>
// #include <time.h>

size_t gen_pts_rsa_1d(double *x,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)
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

    // Set the initial position
    double xn = gsl_rng_uniform (r) * (1 - 2 * radius) + radius;
    x[0] = xn;

    double diameter = 2 * radius;

    size_t valid_pts;
    double dist;
    int k, flag, step;

    step = 0;
    valid_pts = 1;

    while ((valid_pts < npoints) & (step < step_limit))
    {

        xn = gsl_rng_uniform (r) * (1 - 2 * radius) + radius;

        flag = 1;
        for (k = 0; k < valid_pts; k++)
        {

            dist = sqrt( pow(xn - x[k], 2) );
            if (dist < diameter)
            {

                flag = 0;
                break;

            }
        }
        if (flag == 1)
        {

           x[valid_pts] = xn;
           valid_pts += 1;

        }

        step += 1;
        
    }
    

    gsl_rng_free (r);

    return valid_pts;

}