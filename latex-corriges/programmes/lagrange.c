/* Chapitre III : interpolation de Lagrange (algorithme du cours, III-2-1) */
#include <stdio.h>

#define NMAX 50

int main(void)
{
    double x[NMAX], y[NMAX], xp, yp, ys;
    int n1, i, k;

    printf("Nombre de points n+1 : ");
    scanf("%d", &n1);
    for (k = 0; k < n1; k++) {
        printf("x[%d] y[%d] : ", k + 1, k + 1);
        scanf("%lf %lf", &x[k], &y[k]);
    }
    printf("xp : ");
    scanf("%lf", &xp);

    yp = 0.0;
    for (k = 0; k < n1; k++) {
        ys = 1.0;
        for (i = 0; i < n1; i++)
            if (i != k)
                ys = ys * (xp - x[i]) / (x[k] - x[i]);
        yp = yp + y[k] * ys;
    }
    printf("\nxp = %g   yp = %.6f\n", xp, yp);
    return 0;
}
