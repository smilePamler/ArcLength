#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "functiondef.h"

void ComputeVector(double *x, Vec *vold, Vec *v)
{
    // compute v_lambda
    for (int i = 0; i < n; i++)
    {
        v->v_lambda += vold->v_x[i] * x[i];
    }
    v->v_lambda = 1 / (vold->v_lambda - v->v_lambda);

    // compute v_x
    for (int i = 0; i < n; i++)
    {
        v->v_x[i] = -(v->v_lambda * x[i]);
    }

    // normalization
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += v->v_x[i] * v->v_x[i];
    }
    norm += v->v_lambda * v->v_lambda;
    norm = sqrt(norm);
    for (int i = 0; i < n; i++)
    {
        v->v_x[i] = v->v_x[i] / norm;
    }
    v->v_lambda = v->v_lambda / norm;

    // determine the direction
    SelectDirection(vold, v);
    
    return;
}