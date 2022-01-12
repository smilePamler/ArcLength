#include <stdio.h>
#include "functiondef.h"

void LU(double *mat, double *lower, double *upper)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                lower[j*n+i] = 0;
            else
            {
                lower[j*n+i] = mat[j*n+i];
                for (k = 0; k < i; k++)
                {
                    lower[j*n+i] = lower[j*n+i] - lower[j*n+k] * upper[k*n+i];
                }
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                upper[i*n+j] = 0;
            else if (j == i)
                upper[i*n+j] = 1;
            else
            {
                upper[i*n+j] = mat[i*n+j] / lower[i*n+i];
                for (k = 0; k < i; k++)
                {
                    upper[i*n+j] = upper[i*n+j] - ((lower[i*n+k] * upper[k*n+j]) / lower[i*n+i]);
                }
            }
        }
    }
}