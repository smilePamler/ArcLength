#include <stdio.h>
#include <math.h>
#include <string.h>
#include "functiondef.h"

void Jacobian(double *x, double *jac)
{
    //Circuit Parameters

    double R3 = 10000;
    double Rc1 = 2000;
    double Rc2 = 1000;
    double Re = 100;

    //Vcc = 10;
    //Vin = 1.5;

    //BJT Parameters

    double Is = 1 * pow(10, -16);
    double af = 0.99;
    double ar = 0.5;
    double N = 38.78;

    //Ebers-Moll Transistor Model

    // double fe1 = -Is / af * (exp(-N * (x[2] - x[5])) - 1);
    // double fc1 = -Is / ar * (exp(-N * (x[1] - x[5])) - 1);
    // double fe2 = -Is / af * (exp(-N * (x[2] - x[4])) - 1);
    // double fc2 = -Is / ar * (exp(-N * (x[3] - x[4])) - 1);

    // double Ie1 = fe1 - ar * fc1;
    // double Ic1 = fc1 - af * fe1;
    // double Ie2 = fe2 - ar * fc2;
    // double Ic2 = fc2 - af * fe2;

    // f[0] = (x[0] - x[3]) / R3 + (x[0] - x[5]) / Rc1 + Ic1; 

    jac[0 * 8 + 0] = 1 / R3 + 1 / Rc1 + Is * N / ar * (exp(-N * (x[0] - x[4])));
    jac[0 * 8 + 1] = 0 + 0 - Is * N * (exp(-N * (x[1] - x[4]))); 
    jac[0 * 8 + 2] = 0 + 0 + 0;
    jac[0 * 8 + 3] = - 1 / R3 + 0 + 0; 
    jac[0 * 8 + 4] = 0 + 0 + Is * N * (exp(-N * (x[1] - x[4]))) - Is * N / ar * (exp(-N * (x[0] - x[4])));
    jac[0 * 8 + 5] = 0 - 1 / Rc1 + 0;
    jac[0 * 8 + 6] = 0 + 0 + 0;
    jac[0 * 8 + 7] = 0 + 0 + 0;

    // f[1] = x[1] / Re + Ie1 + Ie2;

    jac[1 * 8 + 0] = 0 - Is * N * (exp(-N * (x[0] - x[4]))) + 0;
    jac[1 * 8 + 1] = 1 / Re + Is * N / af * (exp(-N * (x[1] - x[4]))) + Is * N / af * (exp(-N * (x[1] - x[3]))); 
    jac[1 * 8 + 2] = 0 + 0 - Is * N * (exp(-N * (x[2] - x[3])));
    jac[1 * 8 + 3] = 0 + 0 - Is * N / af * (exp(-N * (x[1] - x[3]))) + Is * N * (exp(-N * (x[2] - x[3]))); 
    jac[1 * 8 + 4] = 0 - Is * N / af * (exp(-N * (x[1] - x[4]))) + Is * N * (exp(-N * (x[0] - x[4]))) + 0;
    jac[1 * 8 + 5] = 0 + 0 + 0;
    jac[1 * 8 + 6] = 0 + 0 + 0;
    jac[1 * 8 + 7] = 0 + 0 + 0;

    //f[2] = (x[2] - x[5]) / Rc2 + Ic2;

    jac[2 * 8 + 0] = 0 + 0;
    jac[2 * 8 + 1] = 0 - Is * N * (exp(-N * (x[1] - x[3]))); 
    jac[2 * 8 + 2] = 1 / Rc2 + Is * N / ar * (exp(-N * (x[2] - x[3])));
    jac[2 * 8 + 3] = 0 + Is * N * (exp(-N * (x[1] - x[3]))) - Is * N / ar * (exp(-N * (x[2] - x[3]))); 
    jac[2 * 8 + 4] = 0 + 0;
    jac[2 * 8 + 5] = - 1 / Rc2 + 0;
    jac[2 * 8 + 6] = 0 + 0 + 0;
    jac[2 * 8 + 7] = 0 + 0 + 0;

    //f[3] = (x[3] - x[0]) / R3 - Ic2 - Ie2;

    jac[3 * 8 + 0] = -1 / R3 + 0 + 0;
    jac[3 * 8 + 1] = 0 - Is * N / af * (exp(-N * (x[1] - x[3]))) + Is * N * (exp(-N * (x[1] - x[3]))); 
    jac[3 * 8 + 2] = 0 - Is * N / ar * (exp(-N * (x[2] - x[3]))) + Is * N * (exp(-N * (x[2] - x[3])));
    jac[3 * 8 + 3] = 1 / R3 + Is * N / af * (exp(-N * (x[1] - x[3]))) - Is * N * (exp(-N * (x[2] - x[3]))) - Is * N *(exp(-N * (x[1] - x[3]))) + Is * N / ar * (exp(-N * (x[2] - x[3]))); 
    jac[3 * 8 + 4] = 0 + 0 + 0 + 0;
    jac[3 * 8 + 5] = 0 + 0 + 0 + 0;
    jac[3 * 8 + 6] = 0 + 0 + 0;
    jac[3 * 8 + 7] = 0 + 0 + 0;

    //f[4] = x[4] - Vin;

    jac[4 * 8 + 0] = 0 + 0;
    jac[4 * 8 + 1] = 0 + 0; 
    jac[4 * 8 + 2] = 0 + 0;
    jac[4 * 8 + 3] = 0 + 0; 
    jac[4 * 8 + 4] = 1 + 0;
    jac[4 * 8 + 5] = 0 + 0;
    jac[4 * 8 + 6] = 0 + 0;
    jac[4 * 8 + 7] = 0 + 0;

    //f[5] = x[5] - Vcc;

    jac[5 * 8 + 0] = 0 + 0;
    jac[5 * 8 + 1] = 0 + 0; 
    jac[5 * 8 + 2] = 0 + 0;
    jac[5 * 8 + 3] = 0 + 0; 
    jac[5 * 8 + 4] = 0 + 0;
    jac[5 * 8 + 5] = 1 + 0;
    jac[5 * 8 + 6] = 0 + 0;
    jac[5 * 8 + 7] = 0 + 0;

    //f[6] = (x[5] - x[0]) / Rc1 + (x[5] - x[2]) / Rc2 + x[6];

    jac[6 * 8 + 0] = - 1 / Rc1 + 0 + 0;
    jac[6 * 8 + 1] = 0 + 0 + 0;
    jac[6 * 8 + 2] = 0 - 1 / Rc2 + 0;
    jac[6 * 8 + 3] = 0 + 0 + 0;
    jac[6 * 8 + 4] = 0 + 0 + 0;
    jac[6 * 8 + 5] = 1 / Rc1 + 1 / Rc2 + 0;
    jac[6 * 8 + 6] = 0 + 0 + 1;
    jac[6 * 8 + 7] = 0 + 0 + 0;

    //f[7] = x[7] - Ic1 - Ie1;

    jac[7 * 8 + 0] = 0 + Is * N * (exp(-N * (x[0] - x[4]))) - Is * N / ar * (exp(-N * (x[0] - x[4])));
    jac[7 * 8 + 1] = 0 - Is * N / af * (exp(-N * (x[1] - x[4]))) + Is * N * (exp(-N * (x[1] - x[4])));
    jac[7 * 8 + 2] = 0 + 0 + 0;
    jac[7 * 8 + 3] = 0 + 0 + 0;
    jac[7 * 8 + 4] = 0 + Is * N / af * (exp(-N * (x[1] - x[4]))) - Is * N *(exp(-N * (x[0] - x[4]))) - Is * N * (exp(-N * (x[1] - x[4]))) + Is * N / ar * (exp(-N * (x[0] - x[4])));
    jac[7 * 8 + 5] = 0 + 0 + 0;
    jac[7 * 8 + 6] = 0 + 0 + 0;
    jac[7 * 8 + 7] = 1 + 0 + 0;
}