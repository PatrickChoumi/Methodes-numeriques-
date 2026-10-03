"""Exercice 6 : V'' = a V^(-1/2), écrit comme le système
    V' = W,   W' = a / sqrt(V),
résolu par RK2 (trapèzes). Le second membre n'est pas défini en V = 0 :
on part de x = eps.
Question 1 : V(0) = V'(0) = 0 ; départ avec la solution exacte
    V = C x^(4/3),  C = (9a/4)^(2/3).
Question 2 : V(0) = 0, V(l) = V0 ; tir sur gamma = V'(0) (sécante)."""
import math


def rk2(a, x, v, w, xfin, n):
    """Intègre de x à xfin en n pas ; retourne V(xfin), W(xfin)."""
    h = (xfin - x) / n
    for i in range(n):
        k1v = h * w
        k1w = h * a / math.sqrt(v)
        k2v = h * (w + k1w)
        k2w = h * a / math.sqrt(v + k1v)
        v = v + (k1v + k2v) / 2
        w = w + (k1w + k2w) / 2
    return v, w


def question1(a, l, eps, n):
    c = (9 * a / 4) ** (2 / 3)
    v0 = c * eps ** (4 / 3)                 # V(eps) exact
    w0 = 4 / 3 * c * eps ** (1 / 3)         # V'(eps) exact
    v, w = rk2(a, eps, v0, w0, l, n)
    return v, c * l ** (4 / 3)


def tir(a, l, V0, eps, n, g0, g1, tol=1e-10, kmax=50):
    """Cherche gamma = V'(0) > 0 tel que V(l) = V0."""
    def V_l(g):
        # départ approché : V ~ g x + (4a/3) x^(3/2) / sqrt(g)
        v = g * eps + 4 * a / 3 * eps ** 1.5 / math.sqrt(g)
        w = g + 2 * a * math.sqrt(eps / g)
        return rk2(a, eps, v, w, l, n)[0]
    r0 = V_l(g0) - V0
    for k in range(kmax):
        r1 = V_l(g1) - V0
        if abs(r1) < tol:
            return g1, k + 1
        g0, g1, r0 = g1, g1 - (g1 - g0) * r1 / (r1 - r0), r1
    raise RuntimeError("pas de convergence")


if __name__ == "__main__":
    for n in (100, 200, 400):
        v, exact = question1(a=1.0, l=1.0, eps=1e-2, n=n)
        print(f"n = {n:3d} : V(1) = {v:.8f}   exact = {exact:.8f}"
              f"   erreur = {v - exact:.2e}")
    v, exact = question1(a=1.0, l=1.0, eps=1e-3, n=400)
    print(f"eps = 1e-3, n = 400 : V(1) = {v:.8f}")
    # près de la singularité x = 0, eps et le pas doivent être petits
    for eps, n in ((1e-2, 4000), (1e-3, 16000), (1e-4, 64000)):
        g, k = tir(a=1.0, l=1.0, V0=2.5, eps=eps, n=n, g0=0.1, g1=1.0)
        print(f"tir (eps = {eps:g}, n = {n}) : gamma = {g:.6f}"
              f" en {k} itérations")
