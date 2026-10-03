/* Chapitre V : predicteur Adams-Bashforth d'ordre 2, correcteur RK2
   (formule des trapezes). Exemple : y' = -2 x y, y(0) = 1,
   solution exacte y = exp(-x^2). */
#include <stdio.h>
#include <math.h>

double f(double x, double y) { return -2.0 * x * y; }

int main(void)
{
    double x0 = 0.0, y0 = 1.0, xmax = 1.0;
    int n = 10, k;
    double h = (xmax - x0) / n;
    double x = x0, y = y0, fprec, fk, p;

    /* Demarrage : y1 par Runge-Kutta d'ordre 2 */
    fprec = f(x, y);
    y = y + h / 2.0 * (fprec + f(x + h, y + h * fprec));
    x = x + h;
    printf("%5.2f %12.8f %12.8f\n", x0, y0, exp(-x0 * x0));
    printf("%5.2f %12.8f %12.8f\n", x, y, exp(-x * x));

    for (k = 1; k < n; k++) {
        fk = f(x, y);
        p = y + h / 2.0 * (3.0 * fk - fprec);         /* predicteur AB2 */
        y = y + h / 2.0 * (fk + f(x + h, p));          /* correcteur RK2 */
        fprec = fk;
        x = x + h;
        printf("%5.2f %12.8f %12.8f\n", x, y, exp(-x * x));
    }
    return 0;
}
