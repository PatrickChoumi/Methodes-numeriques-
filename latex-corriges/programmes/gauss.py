"""Exercices 10 et 13 (et chapitre I) : méthode de Gauss avec pivot
partiel (triangularisation), puis résolution du système triangulaire
(remontée). La matrice a est une liste de lignes, modifiée sur place."""
import random


def triangularise(a, b):
    n = len(b)
    for k in range(n - 1):
        # pivot maximal : ligne ip telle que |a[ip][k]| soit maximal
        ip = k
        for i in range(k + 1, n):
            if abs(a[i][k]) > abs(a[ip][k]):
                ip = i
        if a[ip][k] == 0:
            raise ValueError("matrice singulière")
        if ip != k:                       # échange des lignes k et ip
            a[k], a[ip] = a[ip], a[k]
            b[k], b[ip] = b[ip], b[k]
        for i in range(k + 1, n):         # élimination de x_k
            m = a[i][k] / a[k][k]
            for j in range(k + 1, n):
                a[i][j] = a[i][j] - m * a[k][j]
            a[i][k] = 0.0
            b[i] = b[i] - m * b[k]


def remontee(a, b):
    n = len(b)
    x = [0.0] * n
    x[n - 1] = b[n - 1] / a[n - 1][n - 1]
    for k in range(n - 2, -1, -1):
        s = b[k]
        for l in range(k + 1, n):
            s = s - a[k][l] * x[l]
        x[k] = s / a[k][k]
    return x


def gauss(a, b):
    a = [ligne[:] for ligne in a]         # copies des données
    b = b[:]
    triangularise(a, b)
    return remontee(a, b)


if __name__ == "__main__":
    # système de l'exercice I.3 (a11 = 0 : échange de lignes obligatoire)
    x = gauss([[0, 2, 1], [1, 1, 1], [2, 1, 3]], [5, 6, 13])
    for i, xi in enumerate(x, start=1):
        print(f"x{i} = {xi:16.10f}")
    # système 6 x 6 aléatoire de solution connue, avec a11 = 0
    random.seed(1)
    n = 6
    a = [[random.uniform(-1, 1) for j in range(n)] for i in range(n)]
    a[0][0] = 0.0
    xe = [random.uniform(-1, 1) for i in range(n)]
    b = [sum(a[i][j] * xe[j] for j in range(n)) for i in range(n)]
    x = gauss(a, b)
    erreur = max(abs(x[i] - xe[i]) for i in range(n))
    print(f"test 6 x 6 : erreur maximale = {erreur:.2e}")
