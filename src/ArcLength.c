#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "functiondef.h"

int ArcLength(HF *H, int arclength_step, double sfinal, Sol *xout)
{
    const double hmax = sfinal;
    const double hmin = 1e-12;
    const double hstart = 0.1;
    const double tol = 1e-2;
    double h = hstart;
    int converge;
    int endgame = 0;

    double delta = 0;

    Sol *yp; // predict solution
    Sol *yc; // correct solution
    Sol *ycold;

    Vec *vold; // previous tangent vector
    Vec *v;    // tangent vector

    yp = (Sol *)malloc(sizeof(Sol));
    yc = (Sol *)malloc(sizeof(Sol));
    ycold = (Sol *)malloc(sizeof(Sol));
    vold = (Vec *)malloc(sizeof(Vec));
    v = (Vec *)malloc(sizeof(Vec));

    // initialzation
    yp->lambda = 0;
    memset(yp->x, 0, sizeof(yp->x));
    v->v_lambda = 0;
    memset(v->v_x, 0, sizeof(v->v_x));

    FILE *fp = NULL;
    fp = fopen("solution.txt", "w");
    if (fp == NULL)
    {
        printf("File cannot open!\n");
        exit(0);
    }

    // Start the arclength method
    for (int i = 0; i < arclength_step; i++)
    {
        delta = 0;
        converge = 0;

        // fprintf(fp, "Arclength = %d\n", i + 1);

        if (i == 0)
        {
            yc->lambda = H->s0.lambda;
            memcpy(yc->x, H->s0.x, sizeof(yc->x));
            memset(vold->v_x, 0, sizeof(vold->v_x));
            vold->v_lambda = 1;
        }

        memcpy(ycold, yc, sizeof(yc));

        while (converge == 0)
        {
            // fprintf(fp, "step size = %f\n", h);
            Predictor(h, vold, v, yp, yc, H); // get yp

            converge = Corrector(i, v, *yp, yc, H); // get yc

            if (converge == 0)
            {
                if (0.5 * h > hmin)
                {
                    h = 0.5 * h;
                }
                else
                {
                    h = hmin;
                }
                memcpy(yc, ycold, sizeof(ycold));
            }
        }

        if (converge == 1)
        {
            // record the corrected solution
            for (int k = 0; k < n; k++)
            {
                fprintf(fp, "%f ", yc->x[k]);
            }
            fprintf(fp, "%f\n", yc->lambda);

            // check arclength convergence
            if (fabs(1 - yc->lambda) < tol)
            {
                endgame = NR(*yc, xout); // get xout
                break;
            }
            else
            {
                // if (h >= hstart)
                // {
                //     h = pow(0.8 * h * (tol / (delta + 1e-10)), 0.2);
                //     h = h < hmax ? h : hmax;
                // }
                // else
                //     h = hstart;
                h = hstart;
                memcpy(vold, v, sizeof(v));
            }
        }
    }
    fprintf(fp, "\n");

    if (endgame = 1)
    {
        printf("ArcLength converged\n");
    }
    if (endgame == 0)
    {
        printf("Arclength not converged\n");
    }

    fclose(fp);

    free(yp);
    free(yc);
    free(vold);
    free(v);

    return endgame;
}