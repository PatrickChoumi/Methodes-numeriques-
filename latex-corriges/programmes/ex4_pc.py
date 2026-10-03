"""Exercice 4 : y' = -a y - b y^4 + c, y(0) = y0.
Prédicteur : Adams-Bashforth d'ordre 2 ; correcteur : RK2 (trapèzes)."""


def resoudre(a, b, c, y0, xmax, n, affichage=10):
    f = lambda y: -a * y - b * y ** 4 + c      # f ne dépend pas de x
    h = xmax / n
    y = y0
    print(f"{0.0:8.4f} {y:16.10f}")
    # démarrage : y_1 par RK2
    fprec = f(y)
    y = y + h / 2 * (fprec + f(y + h * fprec))
    for k in range(1, n):
        fk = f(y)
        p = y + h / 2 * (3 * fk - fprec)          # prédicteur AB2
        y = y + h / 2 * (fk + f(p))               # correcteur RK2
        fprec = fk
        if (k + 1) % affichage == 0:
            print(f"{(k + 1) * h:8.4f} {y:16.10f}")
    return y


if __name__ == "__main__":
    y2 = resoudre(a=1.0, b=0.5, c=2.0, y0=0.0, xmax=2.0, n=40)
    print(f"y(2) = {y2:.10f}")
