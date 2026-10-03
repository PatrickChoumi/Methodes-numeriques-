/* Exercice 16 : u_tt + a u_t = c^2 u_xx + d u^3 sur ]0,l[,
   u(0,t) = u(l,t) = 0, u(x,0) = f(x), u_t(x,0) = g(x).
   Differences centrees en temps et en espace (schema explicite a 3 niveaux). */
#include <stdio.h>
#include <math.h>

#define MMAX 2001
#define PI 3.14159265358979324

double f(double x) { return sin(PI * x); }   /* deplacement initial (l = 1) */
double g(double x) { (void)x; return 0.0; }  /* vitesse initiale            */

int main(void)
{
    double a, c, d, l, k, tmax, h, lam2, x;
    double uold[MMAX], u[MMAX], unew[MMAX];
    int m, i, j, nt;

    printf("a, c, d, l, m, k, tmax : ");
    if (scanf("%lf %lf %lf %lf %d %lf %lf", &a, &c, &d, &l, &m, &k, &tmax) != 7)
        return 1;
    h = l / m;
    lam2 = (c * k / h) * (c * k / h);
    printf("lambda = c k / h = %g\n", c * k / h);
    if (c * k / h > 1.0) printf("Attention : lambda > 1, schema instable\n");

    /* niveau 0 : u(x,0) = f(x) */
    for (i = 0; i <= m; i++) uold[i] = f(i * h);
    uold[0] = uold[m] = 0.0;
    /* niveau 1 : developpement de Taylor a l'ordre 2 en k,
       u_tt(x,0) = c^2 f'' + d f^3 - a g */
    for (i = 1; i < m; i++) {
        x = i * h;
        u[i] = uold[i] + k * g(x)
             + 0.5 * k * k * (c * c * (uold[i + 1] - 2 * uold[i] + uold[i - 1]) / (h * h)
                              + d * pow(uold[i], 3) - a * g(x));
    }
    u[0] = u[m] = 0.0;
    /* niveaux 2, 3, ... */
    nt = (int)floor(tmax / k + 0.5);
    for (j = 1; j < nt; j++) {
        for (i = 1; i < m; i++)
            unew[i] = (2 * u[i] - (1 - a * k / 2) * uold[i]
                       + lam2 * (u[i + 1] - 2 * u[i] + u[i - 1])
                       + k * k * d * pow(u[i], 3)) / (1 + a * k / 2);
        unew[0] = unew[m] = 0.0;
        for (i = 0; i <= m; i++) { uold[i] = u[i]; u[i] = unew[i]; }
    }
    printf("t = %g\n", nt * k);
    for (i = 0; i <= m; i += (m / 10 > 0 ? m / 10 : 1))
        printf("%8.4f %14.8f\n", i * h, u[i]);
    return 0;
}
