"""Exercice 14, question 2 : schéma de Richardson (saute-mouton)
  T(i,j+1) = T(i,j-1) + 2 lambda (T(i+1,j) - 2 T(i,j) + T(i-1,j)),
démarré avec T(i,1) donné par le schéma explicite de la question 1."""
import math


def mur_richardson(L, alpha, T0, T1, TL, m, dt, tmax):
    dx = L / m
    lam = alpha * dt / dx ** 2
    print(f"lambda = {lam:.4f}")
    bord0 = lambda t: T0 + T1 * math.sin(math.pi * t / 2)
    ancien = [T0] * m + [TL]           # niveau j = 0
    courant = ancien[:]                # niveau j = 1 : schéma explicite
    for i in range(1, m):
        courant[i] = ((1 - 2 * lam) * ancien[i]
                      + lam * (ancien[i + 1] + ancien[i - 1]))
    courant[0], courant[m] = bord0(dt), TL
    nt = round(tmax / dt)
    for j in range(1, nt):             # niveaux j + 1 = 2, 3, ...
        nouveau = courant[:]
        for i in range(1, m):
            nouveau[i] = ancien[i] + 2 * lam * (
                courant[i + 1] - 2 * courant[i] + courant[i - 1])
        nouveau[0], nouveau[m] = bord0((j + 1) * dt), TL
        ancien, courant = courant, nouveau
        if (j + 1) % (nt // 10) == 0:
            print(f"t = {(j + 1) * dt:6.3f}"
                  f"   T(L/2) = {courant[m // 2]:.4e}")


if __name__ == "__main__":
    mur_richardson(L=1, alpha=0.1, T0=20, T1=10, TL=30,
                   m=10, dt=0.01, tmax=3)
