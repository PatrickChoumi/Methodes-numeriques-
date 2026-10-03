"""Exercice 8 : y'' = p(x) y' + q(x) y + g(x), y(a) = c, y(b) = d.
Différences finies centrées (ordre 2), puis méthode itérative de Jacobi.
Exemple (Burden-Faires, exemple 1 de la section 11.3, p. 703) :
p = -2/x, q = 2/x^2, g = sin(ln x)/x^2, 1 <= x <= 2, y(1) = 1, y(2) = 2."""
import math


def differences_finies_jacobi(p, q, g, a, b, c, d, n, eps, kmax):
    h = (b - a) / n
    # équation i : B_i y_(i-1) + A_i y_i + C_i y_(i+1) = D_i
    A, B, C, D = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    for i in range(1, n):
        x = a + i * h
        B[i] = 1 + h / 2 * p(x)
        A[i] = -(2 + h * h * q(x))
        C[i] = 1 - h / 2 * p(x)
        D[i] = h * h * g(x)
    y = [0.0] * (n + 1)          # valeurs initiales nulles à l'intérieur
    y[0], y[n] = c, d
    for k in range(1, kmax + 1):
        z = y[:]                 # Jacobi : nouvelles valeurs à part
        num = den = 0.0
        for i in range(1, n):
            z[i] = (D[i] - B[i] * y[i - 1] - C[i] * y[i + 1]) / A[i]
            num += abs(z[i] - y[i])
            den += abs(z[i])
        y = z
        if num / den < eps:
            return y, k
    return y, kmax


if __name__ == "__main__":
    y, k = differences_finies_jacobi(
        lambda x: -2 / x, lambda x: 2 / x ** 2,
        lambda x: math.sin(math.log(x)) / x ** 2,
        a=1.0, b=2.0, c=1.0, d=2.0, n=10, eps=1e-12, kmax=10000)
    print(f"Itérations : {k}")
    for i, yi in enumerate(y):
        print(f"{1 + i / 10:6.2f} {yi:14.8f}")
