#ifndef __FUNCTIONDEF_H__
#define __FUNCTIONDEF_H__

#define n 8

typedef struct Solution{
    double lambda;
    double x[n];
}Sol;

typedef struct Vector{
    double v_x[n];
    double v_lambda;
}Vec;

typedef struct HomoFun
{
    double G[n*n];
    double f[n];
    double Jac[n*n];
    Sol s0;
}HF;

void Homotopy(double *x0);
int ArcLength(HF *H, int arclength_step, double sfinal, Sol *xout);
void Predictor(double h, Vec *vold, Vec *v, Sol *yp, Sol *yc, HF *H);
void ComputeVector(double *x, Vec *vold, Vec *v);
void Fx(double *x, double *f);
void Jacobian(double *x, double *jac);
void LU(double *mat, double *lower, double *upper);
void Solver(double *L, double *U, double *b, double *solution);
void SelectDirection(Vec *vold, Vec *v);
int Corrector(int stepnum, Vec *v, Sol yp, Sol *yc, HF *H);
int NR(Sol yc, Sol *xout);

#endif