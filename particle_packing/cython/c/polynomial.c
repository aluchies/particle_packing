#include <gsl/gsl_poly.h>

size_t remove_leading_zeros(double *h, size_t n)
{


    size_t m = 0;
    double small = 1e-5;
    size_t i;
    for(i=n; i>0; i--)
    {

        if (fabs(h[i-1]) > small)
        {
            m = i;
            break;
        }
    }
    

    return m;

}


size_t find_roots(double *z, double *h, size_t n)
{

    size_t m;
    m = remove_leading_zeros(&h[0], n);

    if (m > 1)
    {

        gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(m);
        gsl_poly_complex_solve(h, m, w, z);
        gsl_poly_complex_workspace_free (w);

    }

    return m;

}