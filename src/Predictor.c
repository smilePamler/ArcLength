#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functiondef.h"

// LHS1 = lambda*J(x)+(1-lambda)*G
// LHS2 = F(x)-G*(x-a)

// |  LHS1     LHS2  |
// | v_x_p v_lambda_p|

void Predictor(double h, Vec *vold, Vec *v, Sol *yp, Sol *yc, HF *H)
{
    Fx(yc->x, H->f);
    Jacobian(yc->x, H->Jac);

    // Compute LHS1 and LHS2
    double LHS1[n * n];
    double LHS2[n] = { 0 };
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            LHS1[i * n + j] = yc->lambda * H->Jac[i * n + j] + (1 - yc->lambda) * H->G[i * n + j];
            LHS2[i] += H->G[i * n + j] * (yc->x[j] - H->s0.x[j]);
        }
        LHS2[i] = H->f[i] - LHS2[i];
    }

    // LU factorization for LHS1
    double lower[n * n];
    double upper[n * n];
    memset(lower, 0, n * n * sizeof(double));
    memset(upper, 0, n * n * sizeof(double));
    LU(LHS1, lower, upper);

    // Linear solver: solution = LHS1^(-1) * LHS2
    double solution[n];
    memset(solution, 0, n * sizeof(double));
    Solver(lower, upper, LHS2, solution);

    // Find tangent vector & normalization
    ComputeVector(solution, vold, v);

    // Compute predictor
    for (int i = 0; i < n; i++)
    {
        yp->x[i] = yc->x[i] + v->v_x[i] * h;
    }
    yp->lambda = yc->lambda + v->v_lambda * h;

    return;
}