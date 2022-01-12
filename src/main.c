/*
 * @Arclength method for Schmitt Trigger circuit.
 * @author: feng tian
 * @date: July 1, 2021
 * @H(x, lambda) = lambda*F(x) + (1-lambda)G*(x-a)
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "functiondef.h"

int main()
{
    // define a
    double a[n] = {0.5000, 0.4799, 0.9047, 0.6099, 0.6177, 0.8594, 0.8055, 0.5767};
    double *x0 = (double *)malloc(n * sizeof(double));
    memcpy(x0, a, n * sizeof(double));

    Homotopy(x0);

    free(x0);

    system("pause");
    return 0;
}
