"""Exercice 16 : u_tt + a u_t = c^2 u_xx + d u^3 sur ]0, l[,
u(0,t) = u(l,t) = 0, u(x,0) = f(x), u_t(x,0) = g(x).
Différences centrées en temps et en espace (schéma à trois niveaux)."""
import math


def onde(a, c, d, l, m, k, tmax, f, g):
    h = l / m
    lam2 = (c * k / h) ** 2
    beta = a * k / 2
    if c * k / h > 1:
        print("Attention : lambda = c k / h > 1, schéma instable")
    # niveau 0 : u(x,0) = f(x)
    ancien = [f(i * h) for i in range(m + 1)]
    ancien[0] = ancien[m] = 0.0
    # niveau 1 : Taylor à l'ordre 2, u_tt(x,0) = c^2 f'' + d f^3 - a g
    courant = [0.0] * (m + 1)
    for i in range(1, m):
        gi = g(i * h)
        d2f = (ancien[i + 1] - 2 * ancien[i] + ancien[i - 1]) / h ** 2
        courant[i] = (ancien[i] + k * gi
                      + k * k / 2 * (c * c * d2f + d * ancien[i] ** 3
                                     - a * gi))
    # niveaux 2, 3, ...
    nt = round(tmax / k)
    for j in range(1, nt):
        nouveau = [0.0] * (m + 1)
        for i in range(1, m):
            d2u = courant[i + 1] - 2 * courant[i] + courant[i - 1]
            nouveau[i] = (2 * courant[i] - (1 - beta) * ancien[i]
                          + lam2 * d2u
                          + k * k * d * courant[i] ** 3) / (1 + beta)
        ancien, courant = courant, nouveau
    return courant


if __name__ == "__main__":
    f = lambda x: math.sin(math.pi * x)    # déplacement initial (l = 1)
    g = lambda x: 0.0                      # vitesse initiale
    u = onde(a=0, c=1, d=0, l=1, m=50, k=0.01, tmax=1.5, f=f, g=g)
    print(f"a = d = 0       : u(0,5 ; 1,5) = {u[25]:.8f}   (exact 0)")
    u = onde(a=0.4, c=1, d=0, l=1, m=50, k=0.01, tmax=1.5, f=f, g=g)
    w = math.sqrt(math.pi ** 2 - 0.04)
    exact = math.exp(-0.3) * (math.cos(1.5 * w)
                              + 0.2 / w * math.sin(1.5 * w))
    print(f"a = 0,4, d = 0  : u(0,5 ; 1,5) = {u[25]:.8f}"
          f"   (exact {exact:.8f})")
    for m in (50, 100, 200):
        u = onde(a=0.4, c=1, d=-2, l=1, m=m, k=0.5 / m, tmax=1.5,
                 f=f, g=g)
        print(f"a = 0,4, d = -2, m = {m:3d} : u(0,5 ; 1,5) ="
              f" {u[m // 2]:.8f}")
