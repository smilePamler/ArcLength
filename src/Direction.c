#include <stdio.h>
#include <stdlib.h>
#include "functiondef.h"

int sgn(double x)
{
    if (x > 0)
    {
        return 1;
    }
    else if (x == 0)
    {
        return 0;
    }
    else
    {
        return -1;
    }
}

void SelectDirection(Vec *vold, Vec *v)
{
    int sign;
    double ans = 0;
    for (int i = 0; i < n; i++)
    {
        ans += vold->v_x[i] * v->v_x[i];
    }
    ans += vold->v_lambda * v->v_lambda;
    sign = sgn(ans);
    for (int i = 0; i < n; i++)
    {
        v->v_x[i] *= sign;
    }
    v->v_lambda *= sign;

    return;
}