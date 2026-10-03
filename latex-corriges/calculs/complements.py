# Calculs intermediaires detailles cites dans les corriges.
import math
import numpy as np
from scipy.integrate import solve_ivp, quad

def titre(s):
    print("=== " + s)

titre("III, ex. 2 : sommes intermediaires")
x = np.array([1, 2, 3, 4, 5.0]); y = np.array([1.6, 2, 2.3, 2.4, 2.5])
print(" S(1/x) =", (1 / x).sum(), " S(1/x^2) =", (1 / x ** 2).sum(),
      " S(y) =", y.sum(), " S(y/x) =", (y / x).sum())
A = np.array([[5, (1 / x).sum()], [(1 / x).sum(), (1 / x ** 2).sum()]])
det = np.linalg.det(A)
print(" det =", det, " numerateur a0 =", (1 / x ** 2).sum() * y.sum() - (1 / x).sum() * (y / x).sum(),
      " numerateur a1 =", 5 * (y / x).sum() - (1 / x).sum() * y.sum())
a0, a1 = np.linalg.solve(A, [y.sum(), (y / x).sum()])
print(" residus a) :", (y - a0 - a1 / x).round(6).tolist())
print(" S(e^x) =", np.exp(x).sum(), " S(e^x/x) =", (np.exp(x) / x).sum(),
      " S(e^2x) =", np.exp(2 * x).sum(), " S(y e^x) =", (y * np.exp(x)).sum())
G = np.column_stack([np.ones(5), 1 / x, np.exp(x)])
print(" conditionnement de la matrice normale b) :", np.linalg.cond(G.T @ G))
c = np.linalg.solve(G.T @ G, G.T @ y)
print(" residus b) :", (y - G @ c).round(6).tolist())

titre("III, ex. 3 : poids de Lagrange")
X = np.array([93.0, 96.2, 100.0, 104.2, 108.7]); Y = np.array([11.38, 12.80, 14.70, 17.07, 19.91])
def poids(nodes, xp):
    w = []
    for k in range(len(nodes)):
        p = 1.0
        for i in range(len(nodes)):
            if i != k:
                p *= (xp - nodes[i]) / (nodes[k] - nodes[i])
        w.append(p)
    return np.array(w)
w = poids(X, 102.0); print(" L_k(102) =", w.round(6).tolist(), " somme =", w.sum(), " y_p =", w @ Y)
w = poids(Y, 13.5); print(" L_k(13,5) (interpolation inverse) =", w.round(6).tolist(), " x_p =", w @ X)

titre("IV, ex. 3 : valeurs et Romberg")
for h in [0.2, 0.1]:
    n = round(1 / h); xs = 1 + h * np.arange(n + 1); f = 1 / xs
    print(f" h={h}: f =", f.round(6).tolist(), " somme interieure =", f[1:-1].sum())
T2 = 0.1 * (1.5 + 2 * sum(1 / (1 + 0.2 * k) for k in range(1, 5)))
T1 = 0.05 * (1.5 + 2 * sum(1 / (1 + 0.1 * k) for k in range(1, 10)))
print(" Richardson (4 T_0.1 - T_0.2)/3 =", (4 * T1 - T2) / 3, " erreur =", (4 * T1 - T2) / 3 - math.log(2))

titre("IV, ex. 6 : test de la formule (b-a)/18 [f(a) + 8 f(c) + 8 f(d) + f(b)] sur e^x, [0,1]")
F = lambda a, b, f: (b - a) / 18 * (f(a) + 8 * f(a + (b - a) / 4) + 8 * f(a + 3 * (b - a) / 4) + f(b))
S38 = lambda a, b, f: (b - a) / 8 * (f(a) + 3 * f(a + (b - a) / 3) + 3 * f(a + 2 * (b - a) / 3) + f(b))
ex = math.e - 1
print(" formule :", F(0, 1, math.exp), " erreur exact-formule =", ex - F(0, 1, math.exp),
      " prevision (b-a)^5 f''''/11520 ~", 1 / 11520 * math.exp(0.5))
print(" Simpson 3/8 : erreur =", ex - S38(0, 1, math.exp))

titre("IV, ex. 7 : sommes et borne")
xs = np.linspace(0, 1, 11); f = 4 / (1 + xs ** 2)
print(" f =", f.round(6).tolist())
print(" somme impairs =", f[1::2].sum(), " somme pairs interieurs =", f[2:-1:2].sum())
print(" crochet =", f[0] + 4 * f[1::2].sum() + 2 * f[2:-1:2].sum() + f[-1])
t = np.linspace(0, 1, 100001)
f4 = 4 * 24 * (5 * t ** 4 - 10 * t ** 2 + 1) / (1 + t ** 2) ** 5
print(" max |f''''| sur [0,1] =", np.abs(f4).max(), " borne h^4/180 max =", 1e-4 / 180 * np.abs(f4).max())

titre("IV, ex. 2 : differences")
xd = [1.6487, 1.7333, 1.8221, 1.9155, 2.0138, 2.1170]
print(" differences =", [round(xd[i + 1] - xd[i], 4) for i in range(5)])

titre("V, ex. 6 : version trapezes (ordre 2) du schema integro-differentiel")
def integro(N, schema):
    h = 1.0 / N; xg = h * np.arange(N + 1); y = np.zeros(N + 1); y[0] = 1.0
    f = lambda s: 1 - 2 * math.exp(s)
    def I(i):  # trapezes sur [0, x_i] avec y_0..y_i, noyau k = 1
        if i == 0:
            return 0.0
        return h * (0.5 * y[0] + y[1:i].sum() + 0.5 * y[i])
    for i in range(N):
        Fi = 2 * y[i] + f(xg[i]) + I(i)
        if schema == "euler":
            y[i + 1] = y[i] + h * Fi
        else:
            # trapezes en x : y_{i+1} = y_i + h/2 (F_i + F_{i+1}), F_{i+1} lineaire en y_{i+1}
            # F_{i+1} = 2 y_{i+1} + f_{i+1} + h (0.5 y_0 + sum y_1..y_i + 0.5 y_{i+1})
            Sk = h * (0.5 * y[0] + y[1:i + 1].sum()) if i + 1 > 0 else 0.0
            num = y[i] + h / 2 * (Fi + f(xg[i + 1]) + Sk)
            y[i + 1] = num / (1 - h / 2 * (2 + h / 2))
    return y[-1] - math.e
for N in [10, 20, 40, 80]:
    print(f" N={N}: erreur Euler = {integro(N, 'euler'):.3e}   erreur trapezes = {integro(N, 'trap'):.3e}")

titre("Ex. 1 : constante de Wien")
hP, cL, kB = 6.62607015e-34, 299792458.0, 1.380649e-23
print(" hc/(k x_m) =", hP * cL / (kB * 4.965114231744276), "m.K ; g(0) =", (0 - 5) * 1 + 5)

titre("Ex. 2 : racines du polynome de Van der Waals (CO2, T = 300 K)")
R = 0.08206; aw = 3.592; bw = 0.04267
for P in [10.0, 50.0]:
    K = R * 300
    print(f" P={P}: racines =", np.roots([P, -(P * bw + K), aw, -aw * bw]).round(6).tolist())

titre("Ex. 3 : depart y0 = (z-1)/4")
z = 3.9737609; yy = (z - 1) / 4
for k in range(30):
    f = z * (1 - yy) ** 3 - (1 + yy + yy * yy - yy ** 3)
    fp = -3 * z * (1 - yy) ** 2 - (1 + 2 * yy - 3 * yy * yy)
    d = f / fp; yy -= d
    if abs(d) < 1e-12:
        break
print(f" y0 = {(z - 1) / 4:.6f} -> y = {yy:.10f} en {k + 1} iterations")
print(" z(0.3) =", (1 + 0.3 + 0.09 - 0.027) / 0.7 ** 3)

titre("Ex. 4 : point d'equilibre de y' = -y - 0.5 y^4 + 2")
r = [v.real for v in np.roots([-0.5, 0, 0, -1, 2]) if abs(v.imag) < 1e-12 and v.real > 0]
print(" y* =", r)

titre("Ex. 5 : y(1, gamma) en fonction de gamma (a=b=c=1) et seconde solution (a=1, b=1, c=0.5)")
def y1(a, b, c, g):
    def F(x, Y):
        if x == 0:
            return [Y[1], -b * math.exp(c * Y[0]) / (1 + a)]
        return [Y[1], -a / x * Y[1] - b * math.exp(c * Y[0])]
    s = solve_ivp(F, [0, 1], [g, 0], rtol=1e-11, atol=1e-12)
    return s.y[0, -1]
gs = np.linspace(-5, 5, 101)
v = [y1(1, 1, 1, g) for g in gs]
print(" a=b=c=1 : max_gamma y(1,gamma) =", max(v), "atteint pour gamma ~", gs[int(np.argmax(v))])
print(" lambda = b c e^c (forme de Gelfand) : a=b=c=1 ->", math.e, " ; a=b=1, c=0.5 ->", 0.5 * math.exp(0.5))
from scipy.optimize import brentq
G = lambda g: y1(1, 1, 0.5, g) - 1
gg = np.linspace(0, 10, 201); vals = [G(g) for g in gg]
racines = [brentq(G, gg[i], gg[i + 1]) for i in range(200) if vals[i] * vals[i + 1] < 0]
print(" a=b=1, c=0.5 : racines gamma =", [round(r_, 8) for r_ in racines])

titre("Ex. 6, 2) : monotonie de V(l, gamma)")
def ex6_V(gam, a=1.0, l=1.0, eps=1e-6):
    V = gam * eps + 4 * a / 3 * eps ** 1.5 / math.sqrt(gam)
    W = gam + 2 * a * math.sqrt(eps / gam)
    s = solve_ivp(lambda x, Y: [Y[1], a / math.sqrt(Y[0])], [eps, l], [V, W], rtol=1e-11, atol=1e-13)
    return s.y[0, -1]
print(" V(1,gamma) pour gamma = 0.01, 0.1, 0.5, 1, 2 :", [round(ex6_V(g), 6) for g in [0.01, 0.1, 0.5, 1, 2]])

titre("Ex. 7 : tir + predicteur-correcteur sur y''' = 2x - y', y(0)=0, y'(0)=1, y(1)=sin 1 + 1 (exact y''(0)=2)")
def pc_sys(F, x, Y, h, n):
    Y = np.array(Y, float); Fp = F(x, Y)
    Y = Y + h / 2 * (Fp + F(x + h, Y + h * Fp)); x += h          # demarrage RK2
    for _ in range(n - 1):
        Fk = F(x, Y); P = Y + h / 2 * (3 * Fk - Fp)                 # predicteur AB2
        Y = Y + h / 2 * (Fk + F(x + h, P)); Fp = Fk; x += h         # correcteur trapezes
    return Y
F7 = lambda x, Y: np.array([Y[1], Y[2], 2 * x - Y[1]])
target = math.sin(1) + 1
for n in [20, 40, 80]:
    g0, g1 = 0.0, 1.0
    y0 = pc_sys(F7, 0, [0, 1, g0], 1 / n, n)[0] - target
    for it in range(20):
        y1_ = pc_sys(F7, 0, [0, 1, g1], 1 / n, n)[0] - target
        if abs(y1_) < 1e-12:
            break
        g0, g1, y0 = g1, g1 - (g1 - g0) * y1_ / (y1_ - y0), y1_
    print(f" n={n}: gamma = y''(0) = {g1:.8f} (exact 2), erreur = {g1 - 2:.2e}")

titre("Ex. 13, 5) : bissection sur phi obtenue par Nystrom (exact phi = x - 1/2, zero en 0,5)")
def nystrom(n, f, k, a=0.0, b=1.0):
    h = (b - a) / (n - 1); t = a + h * np.arange(n)
    c = np.full(n, 1 / (n - 1)); c[0] = c[-1] = 1 / (2 * (n - 1))
    A = np.eye(n) - (b - a) * np.array([[c[i] * k(t[j], t[i]) for i in range(n)] for j in range(n)])
    ph = np.linalg.solve(A, f(t))
    return lambda x: f(x) + (b - a) * sum(c[i] * k(x, t[i]) * ph[i] for i in range(n))
fz = lambda x: x - 0.5 - 1 / 36
kz = lambda x, t: (x + t) / 3
for n in [5, 9, 17]:
    phi = nystrom(n, fz, kz)
    a_, b_ = 0.0, 1.0
    while b_ - a_ > 1e-10:
        c_ = (a_ + b_) / 2
        if phi(a_) * phi(c_) <= 0:
            b_ = c_
        else:
            a_ = c_
    print(f" n={n}: zero = {(a_ + b_) / 2:.8f} (exact 0,5)")
