# Tests numeriques des algorithmes rediges dans les corriges.
# Chaque fonction applique l'algorithme tel qu'il est ecrit dans le corrige
# et le compare a une solution connue (exacte ou de reference).
import math
import numpy as np
from scipy.integrate import solve_ivp

def titre(s):
    print("=== " + s)

# ---------------------------------------------------------------- Chapitre I
titre("I : remontee, pivot partiel, LU (Crout), L^-1 b, inverse, Cholesky")

def gauss_pivot(A, b):
    A = A.astype(float).copy(); b = b.astype(float).copy(); n = len(b)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(A[k:, k]))
        if p != k:
            A[[k, p]] = A[[p, k]]; b[[k, p]] = b[[p, k]]
        for l in range(k + 1, n):
            m = A[l, k] / A[k, k]
            b[l] -= m * b[k]
            A[l, k + 1:] -= m * A[k, k + 1:]
            A[l, k] = 0.0
    return remontee(A, b)

def remontee(A, b):
    n = len(b); x = np.zeros(n)
    x[n - 1] = b[n - 1] / A[n - 1, n - 1]
    for k in range(n - 2, -1, -1):
        s = b[k]
        for l in range(k + 1, n):
            s -= A[k, l] * x[l]
        x[k] = s / A[k, k]
    return x

def crout(A):
    n = A.shape[0]; L = np.zeros((n, n)); S = np.eye(n)
    for p in range(n):
        for j in range(p, n):              # colonne p de L
            L[j, p] = A[j, p] - sum(L[j, k] * S[k, p] for k in range(p))
        for j in range(p + 1, n):          # ligne p de S
            S[p, j] = (A[p, j] - sum(L[p, k] * S[k, j] for k in range(p))) / L[p, p]
    return L, S

def descente(L, b):
    n = len(b); y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - sum(L[i, j] * y[j] for j in range(i))) / L[i, i]
    return y

def inverse(A):
    n = A.shape[0]; L, S = crout(A); X = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n); e[j] = 1.0
        X[:, j] = remontee(S, descente(L, e))
    return X

def cholesky_M(A):
    n = A.shape[0]; M = np.zeros((n, n))
    for p in range(n):
        s = A[p, p] - sum(M[k, p] ** 2 for k in range(p))
        if s <= 0:
            raise ValueError("matrice non definie positive")
        M[p, p] = math.sqrt(s)
        for j in range(p + 1, n):
            M[p, j] = (A[p, j] - sum(M[k, p] * M[k, j] for k in range(p))) / M[p, p]
    return M

rng = np.random.default_rng(1)
A = rng.uniform(-1, 1, (5, 5)); A[0, 0] = 0.0; xe = rng.uniform(-1, 1, 5); b = A @ xe
print(" Gauss + pivot (a11 = 0) : erreur max =", np.abs(gauss_pivot(A, b) - xe).max())
A = rng.uniform(-1, 1, (5, 5)) + 5 * np.eye(5); b = A @ xe
L, S = crout(A)
print(" Crout : max|LS - A| =", np.abs(L @ S - A).max(),
      " ; resolution LSx=b : erreur =", np.abs(remontee(S, descente(L, b)) - xe).max())
print(" inverse : max|A A^-1 - I| =", np.abs(A @ inverse(A) - np.eye(5)).max())
Ach = np.array([[4, -1, 1], [-1, 4.25, 2.75], [1, 2.75, 3.5]])   # B&F, ex. 4 p. 424
M = cholesky_M(Ach)
print(" Cholesky (exemple B&F p. 424) : M =", M.round(6).tolist(), " max|M^t M - A| =",
      np.abs(M.T @ M - Ach).max())

titre("I : Jacobi et Gauss-Seidel avec test d'arret (systeme de B&F, ex. 1 p. 456)")
Aj = np.array([[10, -1, 2, 0], [-1, 11, -1, 3], [2, -1, 10, -1], [0, 3, -1, 8]], float)
bj = np.array([6, 25, -11, 15], float)

def jacobi(A, b, eps, kmax):
    n = len(b); x = np.zeros(n)
    for k in range(1, kmax + 1):
        xn = np.array([(b[i] - sum(A[i, j] * x[j] for j in range(n) if j != i)) / A[i, i]
                       for i in range(n)])
        d = np.abs(xn - x).sum() / np.abs(xn).sum()
        x = xn
        if d < eps:
            return x, k
    return x, kmax

def gauss_seidel(A, b, eps, kmax):
    n = len(b); x = np.zeros(n)
    for k in range(1, kmax + 1):
        num = 0.0
        for i in range(n):
            xi = (b[i] - sum(A[i, j] * x[j] for j in range(n) if j != i)) / A[i, i]
            num += abs(xi - x[i]); x[i] = xi
        if num / np.abs(x).sum() < eps:
            return x, k
    return x, kmax

for nom, meth in [("Jacobi", jacobi), ("Gauss-Seidel", gauss_seidel)]:
    x, k = meth(Aj, bj, 1e-3, 100)
    print(f" {nom:12s}: {k} iterations, x = {x.round(4).tolist()}")

# --------------------------------------------------------------- Chapitre II
titre("II : Newton avec test |f(x_n)|/m < E, bissection, Lin")

def newton(f, fp, x0, m, E, kmax):
    x = x0
    for k in range(1, kmax + 1):
        x = x - f(x) / fp(x)
        if abs(f(x)) / m < E:
            return x, k
    return x, kmax

def bissection(f, a, b, eps):
    fa = f(a); k = 0
    while b - a > eps:
        c = (a + b) / 2; fc = f(c); k += 1
        if fc == 0:
            return c, k
        if fa * fc < 0:
            b = c
        else:
            a, fa = c, fc
    return (a + b) / 2, k

g = lambda x: (x - 5) * math.exp(x) + 5
gp = lambda x: (x - 4) * math.exp(x)
m = min(abs(gp(4.9)), abs(gp(5.0)))
x, k = newton(g, gp, 5.0, m, 1e-12, 50)
print(f" Newton (Planck) : x = {x:.12f} en {k} iterations, m = {m:.3f}")
x, k = bissection(g, 4.9, 5.0, 1e-10)
print(f" Bissection (Planck) sur [4.9 ; 5] : x = {x:.10f} en {k} iterations")

def lin(a, x0, eps, kmax):
    n = len(a) - 1; x = x0
    for k in range(1, kmax + 1):
        bk = a[0]
        for j in range(1, n):            # b_1 ... b_{n-1}
            bk = a[j] + bk * x
        xn = -a[n] / bk
        if abs(xn - x) < eps:
            return xn, k
        x = xn
    return x, None

P = [1, -6, 11, -6]                      # (x-1)(x-2)(x-3)
for x0 in [0.0, 1.5, 2.5, 3.5]:
    r, k = lin(P, x0, 1e-10, 500)
    print(f" Lin sur x^3-6x^2+11x-6, x0 = {x0} :", f"x = {r:.10f} ({k} it.)" if k else "pas de convergence")

# --------------------------------------------------------------- Chapitre IV
titre("IV, ex. 6 : noyau de Peano de la formule a 4 points (0, 1/4, 3/4, 1)")
w = np.array([1 / 18, 4 / 9, 4 / 9, 1 / 18]); nds = np.array([0, 0.25, 0.75, 1.0])
ss = np.linspace(0, 1, 2001)
K = [((1 - s) ** 4 / 4 - np.sum(w * np.maximum(nds - s, 0) ** 3)) / 6 for s in ss]
print(f" min K = {min(K):.3e}, max K = {max(K):.3e}, integrale |K| = {np.trapezoid(np.abs(K), ss):.6e},"
      f" integrale K = {np.trapezoid(K, ss):.6e}"
      f" (1/11520 = {1 / 11520:.6e})")

# ---------------------------------------------------------------- Chapitre V
titre("V : algorithmes RK2, AB2 sur y' = -2xy, y(0) = 1 (exact e^{-x^2})")

def rk2(f, x, y, h, n):
    for _ in range(n):
        k1 = f(x, y); y = y + h / 2 * (k1 + f(x + h, y + h * k1)); x += h
    return y

def ab2(f, x, y, h, n):
    fprec = f(x, y); y = y + h / 2 * (fprec + f(x + h, y + h * fprec)); x += h
    for _ in range(n - 1):
        fk = f(x, y); y = y + h / 2 * (3 * fk - fprec); fprec = fk; x += h
    return y

f = lambda x, y: -2 * x * y
for n in [10, 20, 40]:
    print(f" n={n}: err RK2 = {rk2(f, 0, 1, 1 / n, n) - math.exp(-1):.3e}"
          f"   err AB2 = {ab2(f, 0, 1, 1 / n, n) - math.exp(-1):.3e}")

titre("V : ordre de la methode 'de Runge d'ordre 4' du cours (L1, L2, L3, L4) et RK4 classique")

def runge_cours(f, x, y, h, n):
    for _ in range(n):
        L1 = h * f(x, y); L2 = h * f(x + h, y + L1); L3 = h * f(x + h, y + L2)
        L4 = h * f(x + h / 2, y + L1 / 2)
        y = y + (L1 + 4 * L4 + L3) / 6; x += h
    return y

def rk4(f, x, y, h, n):
    for _ in range(n):
        L1 = h * f(x, y); L2 = h * f(x + h / 2, y + L1 / 2)
        L3 = h * f(x + h / 2, y + L2 / 2); L4 = h * f(x + h, y + L3)
        y = y + (L1 + 2 * L2 + 2 * L3 + L4) / 6; x += h
    return y

for nom, f_, ex in [("y' = y", lambda x, y: y, math.e), ("y' = -2xy", lambda x, y: -2 * x * y, math.exp(-1))]:
    e = [runge_cours(f_, 0, 1, 1 / n, n) - ex for n in [10, 20, 40, 80]]
    e4 = [rk4(f_, 0, 1, 1 / n, n) - ex for n in [10, 20, 40, 80]]
    print(f" {nom} : Runge (cours) erreurs", [f"{v:.3e}" for v in e], "rapports", [round(e[i] / e[i + 1], 2) for i in range(3)])
    print(f" {nom} : RK4 classique erreurs", [f"{v:.3e}" for v in e4], "rapports", [round(e4[i] / e4[i + 1], 2) for i in range(3)])

titre("V : methode de tir (secante du cours) sur y'' = -y, y(0)=0, y(pi/2)=1 (exact sin x)")

def rk4_sys(F, x, Y, h, n):
    Y = np.array(Y, float)
    for _ in range(n):
        K1 = h * F(x, Y); K2 = h * F(x + h / 2, Y + K1 / 2)
        K3 = h * F(x + h / 2, Y + K2 / 2); K4 = h * F(x + h, Y + K3)
        Y = Y + (K1 + 2 * K2 + 2 * K3 + K4) / 6; x += h
    return Y

def tir_secante(F, a, b, eta1, eta2, g0, g1, n, eps, kmax):
    yb0 = rk4_sys(F, a, [eta1, g0], (b - a) / n, n)[0]
    for k in range(kmax):
        yb1 = rk4_sys(F, a, [eta1, g1], (b - a) / n, n)[0]
        if abs(yb1 - eta2) < eps:
            return g1, k + 1
        g0, g1, yb0 = g1, g1 - (g1 - g0) * (yb1 - eta2) / (yb1 - yb0), yb1
    return g1, None

gam, k = tir_secante(lambda x, Y: np.array([Y[1], -Y[0]]), 0, math.pi / 2, 0, 1, 0.5, 2.0, 50, 1e-10, 20)
print(f" gamma = y'(0) = {gam:.10f} (exact 1) en {k} iterations")

# ------------------------------------------------------- Exercices de fin
titre("Ex. 2 : Van der Waals, CO2 (a = 3.592 L^2 atm/mol^2, b = 0.04267 L/mol)")
R = 0.08206; aw = 3.592; bw = 0.04267
for T, Pr in [(300.0, 10.0), (300.0, 50.0)]:
    K = R * T                      # une mole : second membre R T
    fv = lambda V: Pr * V ** 3 - (Pr * bw + K) * V ** 2 + aw * V - aw * bw
    fpv = lambda V: 3 * Pr * V ** 2 - 2 * (Pr * bw + K) * V + aw
    V = K / Pr
    for it in range(50):
        dV = fv(V) / fpv(V); V -= dV
        if abs(dV) < 1e-12 * abs(V):
            break
    print(f" T = {T} K, P = {Pr} atm : V = {V:.6f} L (gaz parfait {K / Pr:.6f} L), {it + 1} iterations,"
          f" verif (P+a/V^2)(V-b)-RT = {(Pr + aw / V ** 2) * (V - bw) - K:.2e}")

titre("Ex. 5 : tir + Newton (equation variationnelle) sur y'' + (a/x)y' + b e^{cy} = 0")

def ex5(a, bb, c, gam, n=200):
    # inconnues : y, y', z = dy/dgamma, z' ; depart en x = 0 avec les limites
    def F(x, Y):
        y, yp, z, zp = Y
        if x == 0.0:
            ypp = -bb * math.exp(c * y) / (1 + a)
            zpp = -bb * c * math.exp(c * y) * z / (1 + a)
        else:
            ypp = -a / x * yp - bb * math.exp(c * y)
            zpp = -a / x * zp - bb * c * math.exp(c * y) * z
        return np.array([yp, ypp, zp, zpp])
    return rk4_sys(F, 0.0, [gam, 0.0, 1.0, 0.0], 1.0 / n, n)

for (a, bb, c) in [(1.0, 1.0, 0.0), (1.0, 1.0, 0.5), (2.0, 0.5, 1.0)]:
    gam = 1.0
    for it in range(30):
        Y = ex5(a, bb, c, gam)
        d = (Y[0] - 1.0) / Y[2]
        gam -= d
        if abs(d) < 1e-12:
            break
    print(f" a={a}, b={bb}, c={c} : gamma = y(0) = {gam:.10f} ({it + 1} it.), y(1) = {ex5(a, bb, c, gam)[0]:.12f}")
print(" (cas c = 0 : exact gamma = 1 + b/(2(1+a)) =", 1 + 1.0 / 4, ")")

titre("Ex. 6, question 2 : tir sur gamma = V'(0) pour V(l) = V0 (a = 1, l = 1)")

def ex6_V(gam, a=1.0, l=1.0, eps=1e-6, n=2000):
    # depart en x = eps (V ~ gam x + (4a/3) x^{3/2}/sqrt(gam))
    V = gam * eps + 4 * a / 3 * eps ** 1.5 / math.sqrt(gam)
    W = gam + 2 * a * math.sqrt(eps / gam)
    s = solve_ivp(lambda x, Y: [Y[1], a / math.sqrt(Y[0])], [eps, l], [V, W], rtol=1e-11, atol=1e-13)
    return s.y[0, -1]

C = (9 / 4) ** (2 / 3)
print(f" valeur de Child-Langmuir (gamma = 0) : V(1) = {C:.6f}")
V0 = 2.5; g0, g1 = 0.1, 1.0
for it in range(30):
    v0, v1 = ex6_V(g0) - V0, ex6_V(g1) - V0
    g0, g1 = g1, g1 - (g1 - g0) * v1 / (v1 - v0)
    if abs(g1 - g0) < 1e-10:
        break
print(f" V0 = {V0} : gamma = {g1:.8f}, V(1) = {ex6_V(g1):.8f}")

titre("Ex. 9 : Duffing adimensionne, RK4 puis moindres carres de Fourier")
mu, beta, Om = 0.2, 0.1, 1.2           # Y'' + mu Y' + Y + beta Y^3 = cos(Om tau)
F9 = lambda t, Y: np.array([Y[1], -mu * Y[1] - Y[0] - beta * Y[0] ** 3 + math.cos(Om * t)])
Tp = 2 * math.pi / Om; Mpts = 64; Nh = 5
h = Tp / Mpts
Y = rk4_sys(F9, 0.0, [0.0, 0.0], h, 100 * Mpts)     # 100 periodes de transitoire
ts, ys = [], []
t = 100 * Mpts * h
for i in range(Mpts):
    ts.append(t); ys.append(Y[0])
    Y = rk4_sys(F9, t, Y, h, 1); t += h
ts = np.array(ts); ys = np.array(ys)
G = np.column_stack([np.cos(n * Om * ts) for n in range(1, Nh + 1)] +
                    [np.sin(n * Om * ts) for n in range(1, Nh + 1)])
coef_normales = np.linalg.solve(G.T @ G, G.T @ ys)
coef_ortho = np.array([2 / Mpts * (np.cos(n * Om * ts) @ ys) for n in range(1, Nh + 1)] +
                      [2 / Mpts * (np.sin(n * Om * ts) @ ys) for n in range(1, Nh + 1)])
print(" a_n =", coef_normales[:Nh].round(6).tolist())
print(" b_n =", coef_normales[Nh:].round(6).tolist())
print(" ecart equations normales / formules d'orthogonalite :", np.abs(coef_normales - coef_ortho).max())
print(" residu R =", float(np.sum((ys - G @ coef_normales) ** 2)), " ; moyenne des y_i =", ys.mean())

titre("Ex. 17 : Gauss-Seidel sur le cadre (verification avec u = x^2 - y^2, harmonique)")

def gs_cadre(n, m, inside_outer, bc, eps=1e-12, kmax=100000):
    # grille (n+1) x (m+1) de pas h = k ; inside_outer(i,j) True si noeud inconnu
    U = np.zeros((n + 1, m + 1))
    for i in range(n + 1):
        for j in range(m + 1):
            if not inside_outer(i, j):
                U[i, j] = bc(i, j)
    for it in range(1, kmax + 1):
        num = den = 0.0
        for i in range(1, n):
            for j in range(1, m):
                if inside_outer(i, j):
                    new = 0.25 * (U[i + 1, j] + U[i - 1, j] + U[i, j + 1] + U[i, j - 1])
                    num += abs(new - U[i, j]); den += abs(new); U[i, j] = new
        if num / den < eps:
            return U, it
    return U, kmax

# L = 1, l = 0.8, L1 = 0.5, l1 = 0.4 (cadre centre), h = 0.05 : n = 20, m = 16
inc = lambda i, j: 0 < i < 20 and 0 < j < 16 and not (5 <= i <= 15 and 4 <= j <= 12)
U, it = gs_cadre(20, 16, inc, lambda i, j: (0.05 * i) ** 2 - (0.05 * j) ** 2)
err = max(abs(U[i, j] - ((0.05 * i) ** 2 - (0.05 * j) ** 2)) for i in range(21) for j in range(17))
print(f" test harmonique : {it} iterations, erreur max = {err:.2e}")
U, it = gs_cadre(20, 16, inc, lambda i, j: 1.0 if (5 <= i <= 15 and 4 <= j <= 12) else 0.0, eps=1e-10)
print(f" probleme physique (u = 1 dedans, 0 dehors) : {it} iterations,"
      f" min = {U.min():.3f}, max = {U.max():.3f}, u au milieu du bord gauche du cadre (i=2,j=8) = {U[2, 8]:.6f}")

titre("Ex. 17, 4) : u'' + u'/r = 0, u(R1) = 1, u(R2) = 0 (exact ln(r/R2)/ln(R1/R2))")
R1, R2, N = 0.4, 0.8, 40
dr = (R2 - R1) / N; r = R1 + dr * np.arange(N + 1)
Am = np.zeros((N - 1, N - 1)); bm = np.zeros(N - 1)
for i in range(1, N):
    lo, di, up = 1 - dr / (2 * r[i]), -2.0, 1 + dr / (2 * r[i])
    Am[i - 1, i - 1] = di
    if i > 1:
        Am[i - 1, i - 2] = lo
    else:
        bm[0] -= lo * 1.0
    if i < N - 1:
        Am[i - 1, i] = up
ui = np.linalg.solve(Am, bm)
exact = np.log(r[1:-1] / R2) / math.log(R1 / R2)
print(f" N = {N} : erreur max = {np.abs(ui - exact).max():.2e}")

titre("Ex. 18 : Jacobi sur la plaque entaillee (test u = x^2 + y^2, f = 4)")
p, q = 4, 3                               # h = L/(3p), k = l/(2q)
Lx, ly = 1.0, 0.8
n, m = 3 * p, 2 * q; hx, ky = Lx / n, ly / m
def dans_plaque(i, j):                    # noeud strictement interieur au domaine
    if not (0 < i < n and 0 < j < m):
        return False
    return not (p <= i <= 2 * p and j <= q)
ue = lambda i, j: (i * hx) ** 2 + (j * ky) ** 2
U = np.zeros((n + 1, m + 1))
for i in range(n + 1):
    for j in range(m + 1):
        if not dans_plaque(i, j):
            U[i, j] = ue(i, j)
r2 = (hx / ky) ** 2
for it in range(1, 100000):
    V = U.copy(); num = den = 0.0
    for i in range(n + 1):
        for j in range(m + 1):
            if dans_plaque(i, j):
                V[i, j] = (U[i + 1, j] + U[i - 1, j] + r2 * (U[i, j + 1] + U[i, j - 1])
                           - hx ** 2 * 4.0) / (2 * (1 + r2))
                num += abs(V[i, j] - U[i, j]); den += abs(V[i, j])
    U = V
    if num / den < 1e-12:
        break
err = max(abs(U[i, j] - ue(i, j)) for i in range(n + 1) for j in range(m + 1) if dans_plaque(i, j))
print(f" {it} iterations, erreur max = {err:.2e}")

titre("Ex. 19 : triangle, Neumann par points fictifs (test l = h, exact u = a(l - x - y))")
aN, l19, N = 2.0, 1.0, 10
d = l19 / N
U = np.zeros((N + 1, N + 1))
def ghost(U, i, j):
    # valeurs voisines avec points fictifs : u_{-1,j} = u_{1,j} + 2 a dx, u_{i,-1} = u_{i,1} + 2 a dy
    ue_ = U[i + 1, j]
    uw = U[i - 1, j] if i > 0 else U[1, j] + 2 * aN * d
    un = U[i, j + 1]
    us = U[i, j - 1] if j > 0 else U[i, 1] + 2 * aN * d
    return ue_, uw, un, us
for it in range(1, 200000):
    num = den = 0.0
    for i in range(N):
        for j in range(N - i):
            e_, w_, n_, s_ = ghost(U, i, j)
            new = 0.25 * (e_ + w_ + n_ + s_)
            num += abs(new - U[i, j]); den += abs(new); U[i, j] = new
    if num / den < 1e-13:
        break
err = max(abs(U[i, j] - aN * (l19 - i * d - j * d)) for i in range(N + 1) for j in range(N + 1 - i))
print(f" Gauss-Seidel : {it} iterations, erreur max = {err:.2e}")

titre("VI : regressif, Crank-Nicolson, Dufort-Frankel sur u_t = u_xx, u = e^{-pi^2 t} sin(pi x)")

def crout_tridiag(lo, di, up, r):
    n = len(di); l = np.zeros(n); u = np.zeros(n); z = np.zeros(n); x = np.zeros(n)
    l[0] = di[0]; u[0] = up[0] / l[0]; z[0] = r[0] / l[0]
    for i in range(1, n):
        l[i] = di[i] - lo[i] * u[i - 1]
        if i < n - 1:
            u[i] = up[i] / l[i]
        z[i] = (r[i] - lo[i] * z[i - 1]) / l[i]
    x[n - 1] = z[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = z[i] - u[i] * x[i + 1]
    return x

def chaleur(schema, m, k, T):
    h = 1.0 / m; lam = k / h ** 2; x = h * np.arange(m + 1)
    U = np.sin(np.pi * x); nt = int(round(T / k)); Uold = None
    for j in range(nt):
        Ui = U[1:-1]
        if schema == "regressif":
            new = crout_tridiag(np.full(m - 1, -lam), np.full(m - 1, 1 + 2 * lam), np.full(m - 1, -lam), Ui.copy())
        elif schema == "CN":
            rhs = (1 - lam) * Ui + lam / 2 * (U[2:] + U[:-2])
            new = crout_tridiag(np.full(m - 1, -lam / 2), np.full(m - 1, 1 + lam), np.full(m - 1, -lam / 2), rhs)
        elif schema == "DF":
            if Uold is None:                     # demarrage par Crank-Nicolson
                rhs = (1 - lam) * Ui + lam / 2 * (U[2:] + U[:-2])
                new = crout_tridiag(np.full(m - 1, -lam / 2), np.full(m - 1, 1 + lam), np.full(m - 1, -lam / 2), rhs)
            else:
                new = ((1 - 2 * lam) * Uold[1:-1] + 2 * lam * (U[2:] + U[:-2])) / (1 + 2 * lam)
        Uold = U.copy(); U = U.copy(); U[1:-1] = new
    return np.abs(U - np.exp(-np.pi ** 2 * T) * np.sin(np.pi * x)).max()

for sch in ["regressif", "CN", "DF"]:
    print(" " + sch + " :", [f"{chaleur(sch, m_, k_, 0.5):.2e}" for m_, k_ in [(10, 0.01), (20, 0.005), (40, 0.0025)]],
          " (lambda = 1, 2, 4)")
print(" DF avec lambda = 1 fixe (k/h -> 0) :", [f"{chaleur('DF', m_, 1.0 / m_ ** 2, 0.5):.2e}" for m_ in [10, 20, 40]])
