"""Exercice 11, question 1 : Y'' + w0^2 sin Y = 0, Y(0) = Y0, Y'(0) = Z0.
Système Y' = Z, Z' = -w0^2 sin Y, résolu par RK2 (trapèzes)."""
import math


def pendule_rk2(w0, y0, z0, tmax, n, npas):
    """Affiche (t, Y, Z) tous les npas pas."""
    h = tmax / n
    y, z = y0, z0
    print(f"{0.0:16.8f}{y:16.8f}{z:16.8f}")
    for i in range(1, n + 1):
        k1y = h * z
        k1z = -h * w0 ** 2 * math.sin(y)
        k2y = h * (z + k1z)
        k2z = -h * w0 ** 2 * math.sin(y + k1y)
        y = y + (k1y + k2y) / 2
        z = z + (k1z + k2z) / 2
        if i % npas == 0:
            print(f"{i * h:16.8f}{y:16.8f}{z:16.8f}")
    return y, z


if __name__ == "__main__":
    pendule_rk2(w0=1.0, y0=2 * math.pi / 9, z0=0.0,
                tmax=4.0, n=4000, npas=1000)
