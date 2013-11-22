#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <math.h>

// #include <stdlib.h>
// #include <time.h>


unsigned int metro_md_3d(double *x, double *y, double *z,
    double radius, size_t npoints, unsigned int step_limit,
    unsigned long randSeed)
{

    /* Setup GSL random number generator */
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    /* Set the seed */
    // srand ( time(NULL) );
    // unsigned long randSeed = rand();
    gsl_rng_set(r, randSeed);

    double diameter = 2 * radius;

    double dist, dx, dy, dz, xn, yn, zn;
    unsigned int step, i, k, flag, success_steps;

    step = 0;
    success_steps = 0;
    while (step < step_limit)
    {

        i = step % npoints;

        /* Generate new position */
        while (1)
        {

            dx = diameter * (gsl_rng_uniform (r) - 0.5);
            xn = x[i] + dx;

            dy = diameter * (gsl_rng_uniform (r) - 0.5);
            yn = y[i] + dy;

            dz = diameter * (gsl_rng_uniform (r) - 0.5);
            zn = z[i] + dz;

            if (((xn > radius) & (xn < 1 - radius)) & ((yn > radius) & (yn < 1 - radius)) & ((zn > radius) & (zn < 1 - radius)))
            {
               break;
            }

        }

        /* Determine if new position overlaps with other positions */
        flag = 1;
        for (k = 0; k < npoints; k++)
        {
            dist = sqrt( pow(xn - x[k], 2) + pow(yn - y[k], 2) + pow(zn - z[k], 2) );
            if ((dist < diameter) & (i != k))
            {
                flag = 0;
                break;
            }


        }

        if (flag == 1)
        {

            x[i] = xn;
            y[i] = yn;
            z[i] = zn;
            success_steps = success_steps + 1;

        }

        step = step + 1;


    }

    gsl_rng_free (r);

    return success_steps;

}




unsigned int metro_pd_3d(double *x, double *y, double *z,
    double *radius, size_t npoints, int step_limit,
    unsigned long randSeed)
{

    /* Setup GSL random number generator */
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    /* Set the seed */
    // srand ( time(NULL) );
    // unsigned long randSeed = rand();
    gsl_rng_set(r, randSeed);


    double dist, dx, dy, dz, xn, yn, zn, diameter;
    unsigned int step, i, k, flag, success_steps;

    step = 0;
    success_steps = 0;
    while (step < step_limit)
    {

        i = step % npoints;

        /* Generate new position */
        while (1)
        {

            diameter = 2 * radius[i];

            dx = diameter * (gsl_rng_uniform (r) - 0.5);
            xn = x[i] + dx;

            dy = diameter * (gsl_rng_uniform (r) - 0.5);
            yn = y[i] + dy;

            dz = diameter * (gsl_rng_uniform (r) - 0.5);
            zn = z[i] + dz;

            if (((xn > radius[i]) & (xn < 1 - radius[i])) & ((yn > radius[i]) & (yn < 1 - radius[i])) & ((zn > radius[i]) & (zn < 1 - radius[i])))
            {
               break;
            }

        }

        /* Determine if new position overlaps with other positions */
        flag = 1;
        for (k = 0; k < npoints; k++)
        {
            dist = sqrt( pow(xn - x[k], 2) + pow(yn - y[k], 2) + pow(zn - z[k], 2) );
            if ((dist < (radius[i] + radius[k])) & (i != k))
            {
                flag = 0;
                break;
            }


        }

        if (flag == 1)
        {

            x[i] = xn;
            y[i] = yn;
            z[i] = zn;
            success_steps += 1;

        }

        step = step + 1;


    }

    gsl_rng_free (r);

    return success_steps;

}