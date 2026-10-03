"""Exercice 14, question 1 : mur d'épaisseur L, dT/dt = alpha d2T/dx2,
T(x,0) = T0, T(0,t) = T0 + T1 sin(pi t / 2), T(L,t) = TL.
Schéma explicite (différences progressives en temps)."""
import math


def mur_explicite(L, alpha, T0, T1, TL, m, dt, tmax):
    dx = L / m
    lam = alpha * dt / dx ** 2
    print(f"lambda = {lam:.4f}")
    if lam > 0.5:
        print("Attention : lambda > 1/2, schéma instable")
    U = [T0] * m + [TL]            # condition initiale (TL dès t > 0)
    nt = round(tmax / dt)
    for j in range(1, nt + 1):
        V = U[:]
        for i in range(1, m):
            V[i] = (1 - 2 * lam) * U[i] + lam * (U[i + 1] + U[i - 1])
        V[0] = T0 + T1 * math.sin(math.pi * j * dt / 2)
        V[m] = TL
        U = V
    return [(i * dx, U[i]) for i in range(m + 1)]


if __name__ == "__main__":
    profil = mur_explicite(L=1, alpha=1, T0=20, T1=0, TL=100,
                           m=10, dt=0.004, tmax=2)
    for x, T in profil:
        print(f"{x:6.2f} {T:12.6f}")
    profil = mur_explicite(L=1, alpha=0.1, T0=20, T1=10, TL=30,
                           m=10, dt=0.04, tmax=3)
    for x, T in profil:
        print(f"{x:6.2f} {T:12.6f}")
    profil = mur_explicite(L=1, alpha=0.1, T0=20, T1=10, TL=30,
                           m=10, dt=0.06, tmax=3)
    print("dt = 0,06 :", [f"{T:.3e}" for x, T in profil[1:3]])
