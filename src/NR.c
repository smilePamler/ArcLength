#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functiondef.h"

int NR(Sol yc, Sol *xout)
{
    double error = 1e-6;
    int converge = 1;

    double *jac = (double *)malloc(sizeof(double) * n * n);
    double *f = (double *)malloc(sizeof(double) * n);
    double *dx = (double *)malloc(sizeof(double) * n);

    memset(f, 0, n * sizeof(double));
    memset(jac, 0, n * n * sizeof(double));

    for (int i = 0; i < n; i++)
        xout->x[i] = yc.x[i];

    for (int i = 0; i < 100; i++)
    {
        // Compute jacobi matrix of F(xout)
        Fx(xout->x, f);
        Jacobian(xout->x, jac);

        // LU factorization for *jac
        double lower[n * n];
        double upper[n * n];
        memset(lower, 0, n * n * sizeof(double));
        memset(upper, 0, n * n * sizeof(double));
        LU(jac, lower, upper);

        for (int j = 0; j < n; j++)
        {
            f[j] *= -1;
        }

        // Linear solver
        memset(dx, 0, n * sizeof(double));
        Solver(lower, upper, f, dx);

        for (int j = 0; j < n; j++)
        {
            xout->x[j] += dx[j];
        }

        for (int k = 0; k < n; k++)
        {
            if (fabs(dx[k]) > error)
            {
                converge = 0;
                break;
            }
        }
        if (converge == 1)
        {
            printf("final NR iter.no = %d\n", i + 1);
            break;
        }
    }

    return converge;
}