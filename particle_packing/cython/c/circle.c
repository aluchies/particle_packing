#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <math.h>

// #include <stdlib.h>
// #include <time.h>



double circle_overlap(double *rA, double radiiA, double *rB, double radiiB)
{

    double rAB[2];
    rAB[0] = rB[0] - rA[0];
    rAB[1] = rB[1] - rA[1];

    double F;
    F = (pow(rAB[0], 2) + pow(rAB[1], 2)) / pow(radiiA + radiiB, 2);

    return F;

}







size_t gen_pts_rsa_2d(double *x, double *y,
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
    double yn = gsl_rng_uniform (r) * (1 - 2 * radius) + radius;
    x[0] = xn;
    y[0] = yn;

    double diameter = 2 * radius;

    size_t valid_pts;
    double F;
    int k, flag, step;

    step = 0;
    valid_pts = 1;

    double rA[2];
    double rB[2];

    while ((valid_pts < npoints) & (step < step_limit))
    {

        xn = gsl_rng_uniform (r) * (1 - 2 * radius) + radius;
        yn = gsl_rng_uniform (r) * (1 - 2 * radius) + radius;

        flag = 1;
        for (k = 0; k < valid_pts; k++)
        {

            rA[0] = x[k];
            rA[1] = y[k];
            rB[0] = xn;
            rB[1] = yn;

            F = circle_overlap(&rA[0], radius, &rB[0], radius);

            if (F < 1.0)
            {

                flag = 0;
                break;

            }
        }
        if (flag == 1)
        {

           x[valid_pts] = xn;
           y[valid_pts] = yn;
           valid_pts += 1;

        }

        step += 1;
        
    }
    

    gsl_rng_free (r);

    return valid_pts;

}






unsigned int metro_md_2d(double *x, double *y,
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

    double F, dx, dy, xn, yn;
    unsigned int step, i, k, flag, success_steps;

    double rA[2];
    double rB[2];

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


            if (((xn > radius) & (xn < 1 - radius)) & ((yn > radius) & (yn < 1 - radius)) )
            {
               break;
            }

        }

        /* Determine if new position overlaps with other positions */
        flag = 1;
        for (k = 0; k < npoints; k++)
        {

            rA[0] = x[k];
            rA[1] = y[k];
            rB[0] = xn;
            rB[1] = yn;

            F = circle_overlap(&rA[0], radius, &rB[0], radius);

            if ((F < 1.0) & (i != k))
            {
                flag = 0;
                break;
            }


        }

        if (flag == 1)
        {

            x[i] = xn;
            y[i] = yn;
            success_steps = success_steps + 1;

        }

        step = step + 1;


    }

    gsl_rng_free (r);

    return success_steps;

}




unsigned int metro_pd_2d(double *x, double *y,
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


    double F, dx, dy, xn, yn, diameter;
    unsigned int step, i, k, flag, success_steps;

    double rA[2];
    double rB[2];

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


            if (((xn > radius[i]) & (xn < 1 - radius[i])) & ((yn > radius[i]) & (yn < 1 - radius[i])) )
            {
               break;
            }

        }

        /* Determine if new position overlaps with other positions */
        flag = 1;
        for (k = 0; k < npoints; k++)
        {

            rA[0] = x[k];
            rA[1] = y[k];
            rB[0] = xn;
            rB[1] = yn;

            F = circle_overlap(&rA[0], radius[i], &rB[0], radius[k]);

            if ((F < 1.0) & (i != k))
            {
                flag = 0;
                break;
            }


        }

        if (flag == 1)
        {

            x[i] = xn;
            y[i] = yn;
            success_steps += 1;

        }

        step = step + 1;


    }

    gsl_rng_free (r);

    return success_steps;

}