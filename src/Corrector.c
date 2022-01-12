#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "functiondef.h"

// LHS1 = J(x_c)+(1-lambda_c)*G*lambda_c^(-1)
// LHS2 = [F(x_c)-G*(x_c-a)]*lambda_c^(-1)

// |LHS1     LHS2  |
// |v_x    v_lambda|

int Corrector(int stepnum, Vec *v, Sol yp, Sol *yc, HF *H)
{
    double dlambda;
    double dx[n];

    const int NR_iter = 15;
    const double tol = 1e-6;
    int converge;

    // Let the initial value of yc equal to yp
    yc->lambda = yp.lambda;
    memcpy(yc->x, yp.x, sizeof(yp.x));

    // FILE *fp;
    // fp = fopen("NR_sol.txt", "a");
    // if (fp == NULL)
    // {
    //     printf("File cannot open\n");
    //     exit(0);
    // }
    // fprintf(fp, "Arclength = %d:\n", stepnum+1);

    // Newton-Raphson
    for (int step = 0; step < NR_iter; step++)
    {
        // DebugInfo
        // fprintf(fp, "第%d步NR：\n", step + 1);

        dlambda = 0;
        memset(dx, 0, n * sizeof(double));
        converge = 1;

        // Compute jacobi matrix of F(x_c)
        Fx(yc->x, H->f);
        Jacobian(yc->x, H->Jac);

        // Compute LHS1 and LHS2
        double LHS1[n * n];
        double LHS2[n] = {0};
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                LHS1[i * n + j] = yc->lambda * H->Jac[i * n + j] + (1 - yc->lambda) * H->G[i * n + j];
                LHS2[i] += H->G[i * n + j] * (yc->x[j] - H->s0.x[j]);
            }
            LHS2[i] = H->f[i] - LHS2[i];
        }

        // Compute RHS and hyperplane
        double RHS[n] = {0};
        double hyperplane = 0;
        double temp1;
        double temp2 = 0;
        for (int i = 0; i < n; i++)
        {
            temp1 = 0;
            for (int j = 0; j < n; j++)
            {
                temp1 += H->G[i * n + j] * (yc->x[j] - H->s0.x[j]);
            }
            RHS[i] = -yc->lambda * H->f[i] - (1 - yc->lambda) * temp1;
            temp2 += v->v_x[i] * (yc->x[i] - yp.x[i]);
        }
        hyperplane = -temp2 - v->v_lambda * (yc->lambda - yp.lambda);

        // LU factorization for LHS1
        double lower[n * n];
        double upper[n * n];
        memset(lower, 0, n * n * sizeof(double));
        memset(upper, 0, n * n * sizeof(double));
        LU(LHS1, lower, upper);

        // Linear solver: solution2 = LHS1^(-1) * LHS2  solution1 = LHS1^(-1) * RHS
        double solution1[n];
        double solution2[n];
        memset(solution1, 0, n * sizeof(double));
        memset(solution2, 0, n * sizeof(double));
        Solver(lower, upper, RHS, solution1);
        Solver(lower, upper, LHS2, solution2);

        // Compute dx and dlambda
        double temp3 = 0;
        double temp4 = 0;
        for (int k = 0; k < n; k++)
        {
            temp3 += v->v_x[k] * solution1[k];
            temp4 += v->v_x[k] * solution2[k];
        }
        dlambda = (hyperplane - temp3) /  (v->v_lambda - temp4);

        for (int i = 0; i < n; i++)
        {
            dx[i] = solution1[i] - dlambda * solution2[i];
        }

        // Update the corrector solution and compute delta
        for (int i = 0; i < n; i++)
        {
            yc->x[i] += dx[i];
        }
        yc->lambda += dlambda;

        // DEBUG INFO
        // fprintf(fp, "第%d步NR解：\n", step+1);
        // for (int i = 0; i < n; i++)
        // {
        //     fprintf(fp, "dx[%d] = %f ", i, dx[i]);
        // }
        // fprintf(fp, "dlambda = %f\n", dlambda);

        // Check NR's convergence
        for (int i = 0; i < n + 1; i++)
        {
            if (i < n)
            {
                if (fabs(dx[i]) > tol)
                {
                    converge = 0;
                    break;
                }
            }
            else
            {
                if (fabs(dlambda) > tol)
                {
                    converge = 0;
                    break;
                }
            }
        }
        if (converge == 1)
        {
            printf("NR_iter = %d\n", step + 1);
            break;
        }
    }

    // fprintf(fp, "\n");
    // fclose(fp);

    if (converge == 0)
    {
        printf("NR iteration not converged\n");
    }

    return converge;
}