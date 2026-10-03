"""Chapitre V : prédicteur Adams-Bashforth d'ordre 2,
correcteur RK2 (trapèzes).
Exemple : y' = -2 x y, y(0) = 1, solution exacte y = exp(-x^2)."""
import math


def predicteur_correcteur(f, x0, y0, xmax, n):
    """Retourne la liste des couples (x_k, y_k), k = 0, ..., n."""
    h = (xmax - x0) / n
    x, y = x0, y0
    points = [(x, y)]
    # démarrage : y_1 par Runge-Kutta d'ordre 2
    fprec = f(x, y)
    y = y + h / 2 * (fprec + f(x + h, y + h * fprec))
    x = x0 + h
    points.append((x, y))
    for k in range(1, n):
        fk = f(x, y)
        p = y + h / 2 * (3 * fk - fprec)        # prédicteur AB2
        y = y + h / 2 * (fk + f(x + h, p))      # correcteur RK2
        fprec = fk
        x = x0 + (k + 1) * h
        points.append((x, y))
    return points


if __name__ == "__main__":
    f = lambda x, y: -2.0 * x * y
    for x, y in predicteur_correcteur(f, 0.0, 1.0, 1.0, 10):
        print(f"{x:5.2f} {y:12.8f} {math.exp(-x * x):12.8f}")
    for n in (10, 20, 40):
        y1 = predicteur_correcteur(f, 0.0, 1.0, 1.0, n)[-1][1]
        print(f"n = {n:3d}   erreur en x = 1 : {y1 - math.exp(-1):.3e}")
