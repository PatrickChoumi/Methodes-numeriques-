"""Exercice 3 : facteur de compressibilité
    z = (1 + y + y^2 - y^3) / (1 - y)^3,   y = b / (4 v).
On cherche y (puis b) par Newton-Raphson sur
    f(y) = z (1 - y)^3 - (1 + y + y^2 - y^3) = 0."""


def desantis(z, v, y0, eps=1e-12, kmax=50):
    """Retourne (y, b, nombre d'itérations)."""
    y = y0
    for k in range(1, kmax + 1):
        f = z * (1 - y) ** 3 - (1 + y + y ** 2 - y ** 3)
        fp = -3 * z * (1 - y) ** 2 - (1 + 2 * y - 3 * y ** 2)
        dy = f / fp
        y = y - dy
        if abs(dy) < eps:
            return y, 4 * v * y, k
    raise RuntimeError(f"pas de convergence en {kmax} itérations")


if __name__ == "__main__":
    z, v = 3.9737609, 1.0                 # valeur de z pour y = 0,3
    for y0 in (0.5, (z - 1) / 4):
        y, b, k = desantis(z, v, y0)
        print(f"y0 = {y0:.4f} : y = {y:.10f}   b = 4 v y = {b:.10f}"
              f"   ({k} itérations)")
