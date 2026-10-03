"""Exercice 11, questions 2 et 3.
T(Y0)/T0 = (2/pi) * intégrale de 0 à pi/2 de dx / sqrt(1 - A sin^2 x),
avec A = sin^2(Y0/2), par la formule composite des trapèzes (sans
tableaux), puis moindres carrés T/T0 = a + b Y0^2."""
import math


def trapezes(f, a, b, n):
    h = (b - a) / n
    s = (f(a) + f(b)) / 2
    for k in range(1, n):
        s = s + f(a + k * h)
    return h * s


def periode_relative(y0, n):
    A = math.sin(y0 / 2) ** 2
    f = lambda x: 1 / math.sqrt(1 - A * math.sin(x) ** 2)
    return 2 / math.pi * trapezes(f, 0.0, math.pi / 2, n)


def droite_moindres_carres(X, Y):
    """(a, b) tels que Y = a + b X : équations normales, Cramer."""
    m = len(X)
    s1, s2 = sum(X), sum(x * x for x in X)
    sy, sxy = sum(Y), sum(x * y for x, y in zip(X, Y))
    det = m * s2 - s1 * s1
    return (s2 * sy - s1 * sxy) / det, (m * sxy - s1 * sy) / det


if __name__ == "__main__":
    Y0 = [math.pi / 18, math.pi / 9, math.pi / 6, 2 * math.pi / 9]
    for y0 in Y0:
        print(f"Y0 = {y0:.4f} : T/T0 = {periode_relative(y0, 4):.8f}"
              f" (n = 4)   {periode_relative(y0, 8):.8f} (n = 8)")
    X = [y * y for y in Y0]
    a, b = droite_moindres_carres(X, [1.002, 1.008, 1.020, 1.030])
    print(f"tableau de l'énoncé : a = {a:.6f}   b = {b:.6f}")
    exactes = [periode_relative(y, 8) for y in Y0]
    a, b = droite_moindres_carres(X, exactes)
    print(f"valeurs exactes     : a = {a:.6f}   b = {b:.6f}")
