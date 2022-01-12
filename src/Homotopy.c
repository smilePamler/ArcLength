#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functiondef.h"

void Homotopy(double *x0)
{
    int success;
    double sfinal = 30;

    Sol *xout = (Sol *)malloc(sizeof(Sol)); // final solution
    memset(xout, 0, sizeof(xout));

    // Homotopy function: H(x) = lambda*F(x) + (1-lambda)*G*(x-x0)

    HF *H = (HF *)malloc(sizeof(HF));

    // initialization
    memset(H->f, 0, n * sizeof(double));
    memset(H->Jac, 0, n * n * sizeof(double));
    memset(H->G, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        H->G[i * n + i] = 1e-3;
        H->s0.x[i] = x0[i];
    }
    H->s0.lambda = 0;

    // call Arclength interface
    success = ArcLength(H, 200, sfinal, xout);

    if (success == 1)
    {
        FILE *fp = fopen("final_sol.txt", "w");
        fprintf(fp, "The DC operating point is:\n");
        for (int i = 0; i < n; i++)
        {
            fprintf(fp, "x[%d] = %f ", i, xout->x[i]);
        }
        fclose(fp);
    }

    free(xout);
    free(H);
}