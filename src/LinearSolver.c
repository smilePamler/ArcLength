#include <stdio.h>
#include <stdlib.h>
#include "functiondef.h"

void Solver(double *L, double *U, double *b, double *solution)
{
    // Lv = b

    double *v = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        double rhs = b[i];
        for (int j = 0; j < i; j++)
        {
            rhs -= L[i * n + j] * v[j];
        }
        v[i] = rhs / L[i * n + i];
    }

    // Ux = v

    for (int i = n - 1; i >= 0; i--)
    {
        double rhs = v[i];
        for (int j = n - 1; j > i; j--)
        {
            rhs -= U[i * n + j] * solution[j];
        }
        solution[i] = rhs / U[i * n + i];
    }

    free(v);
}