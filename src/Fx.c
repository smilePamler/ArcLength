#include <stdio.h>
#include <math.h>
#include "functiondef.h"

void Fx(double *x, double *f)
{
    //Circuit Parameters

    double R3 = 10000.0;
    double Rc1 = 2000.0;
    double Rc2 = 1000.0;
    double Re = 100.0;

    double Vcc = 10.0;
    double Vin = 1.5;

    //BJT Parameters

    double Is = 1 * pow(10, -16);
    double af = 0.99;
    double ar = 0.5;
    double N = 38.78;
    
    //Ebers-Moll Transistor Model

    double fe1 = -Is / af * (exp(-N * (x[1] - x[4])) - 1);
    double fc1 = -Is / ar * (exp(-N * (x[0] - x[4])) - 1);
    double fe2 = -Is / af * (exp(-N * (x[1] - x[3])) - 1);
    double fc2 = -Is / ar * (exp(-N * (x[2] - x[3])) - 1);

    double Ie1 = fe1 - ar * fc1;
    double Ic1 = fc1 - af * fe1;
    double Ie2 = fe2 - ar * fc2;
    double Ic2 = fc2 - af * fe2;

    //System of Equations for Schmitt Trigger Circuit 1 (Modified Nodal Equations)

    f[0] = (x[0] - x[3]) / R3 + (x[0] - x[5]) / Rc1 + Ic1; 
    f[1] = x[1] / Re + Ie1 + Ie2;
    f[2] = (x[2] - x[5]) / Rc2 + Ic2;
    f[3] = (x[3] - x[0]) / R3 - Ic2 - Ie2;
    f[4] = x[4] - Vin;
    f[5] = x[5] - Vcc;
    f[6] = (x[5] - x[0]) / Rc1 + (x[5] - x[2]) / Rc2 + x[6];
    f[7] = x[7] - Ic1 - Ie1;
}