"""Calculs numériques des corrigés (vérification de toutes les valeurs citées)."""
import numpy as np
from fractions import Fraction as Fr
from scipy.special import ellipk
from scipy.integrate import solve_ivp

np.set_printoptions(precision=6, suppress=True)
print("=== III, exercice 2 : moindres carrés")
x = np.array([1, 2, 3, 4, 5.]); y = np.array([1.6, 2, 2.3, 2.4, 2.5])
for nom, G in [("a) {1, 1/x}", np.c_[np.ones(5), 1/x]),
               ("b) {1, 1/x, e^x}", np.c_[np.ones(5), 1/x, np.exp(x)])]:
    A = G.T @ G; b = G.T @ y
    c = np.linalg.solve(A, b)
    r = y - G @ c
    print(nom); print(" A =", A.tolist()); print(" b =", b.tolist())
    print(" coefficients =", c, " R =", (r**2).sum())
    print(" valeurs ajustées =", G @ c)

print("=== III, exercice 3 : Lagrange")
X = [93.0, 96.2, 100.0, 104.2, 108.7]; Y = [11.38, 12.80, 14.70, 17.07, 19.91]
def lag(xs, ys, xp):
    s = 0
    for k in range(len(xs)):
        t = ys[k]
        for i in range(len(xs)):
            if i != k: t *= (xp - xs[i]) / (xs[k] - xs[i])
        s += t
    return s
print(" y(102) =", lag(X, Y, 102)); print(" x(13,5) =", lag(Y, X, 13.5))

print("=== IV, exercice 2 : vitesses")
t = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5]); xt = np.array([1.6487, 1.7333, 1.8221, 1.9155, 2.0138, 2.1170])
print(" écart max |x - e^(t/2)| =", np.abs(xt - np.exp(t/2)).max())
v = np.diff(xt) / 0.1; tm = (t[:-1] + t[1:]) / 2
for a, b_ in zip(tm, v): print(f"  t={a:.2f} v={b_:.3f}  0.5 e^(t/2)={0.5*np.exp(a/2):.4f}")

print("=== IV, exercice 3 : trapèzes 1/x")
for h in (0.2, 0.1):
    n = round(1/h); xs = np.linspace(1, 2, n+1); f = 1/xs
    T = h*(f[0]/2 + f[1:-1].sum() + f[-1]/2)
    print(f" h={h}: T={T:.6f}  erreur={T-np.log(2):.6f}  borne h^2/12*(b-a)*max|f''|={h*h/12*2:.6f}")
print(" ln 2 =", np.log(2))

print("=== IV, exercice 6 : poids, noeuds 0, 1/4, 3/4, 1 sur [0,1]")
nodes = [Fr(0), Fr(1, 4), Fr(3, 4), Fr(1)]
def integ_poly(coeffs):  # coeffs ascending
    return sum(c / (i+1) for i, c in enumerate(coeffs))
def polymul(p, q):
    r = [Fr(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b_ in enumerate(q): r[i+j] += a*b_
    return r
W = []
for k, xk in enumerate(nodes):
    p = [Fr(1)]
    for i, xi in enumerate(nodes):
        if i != k: p = polymul(p, [-xi/(xk-xi), 1/(xk-xi)])
    W.append(integ_poly(p))
print(" poids =", W, " somme =", sum(W))
for d in range(6):
    exact = Fr(1, d+1); q = sum(w * xk**d for w, xk in zip(W, nodes))
    print(f"  x^{d}: formule={q} exact={exact} {'OK' if q == exact else 'ECART '+str(q-exact)}")

print("=== IV, exercice 7 : Simpson 4/(1+x^2)")
h = 0.1; xs = np.linspace(0, 1, 11); f = 4/(1+xs**2)
S = h/3*(f[0] + 4*f[1:-1:2].sum() + 2*f[2:-1:2].sum() + f[-1])
print(" valeurs f =", f); print(f" S={S:.10f} pi={np.pi:.10f} erreur={S-np.pi:.3e}")

print("=== Fin 1 : Planck / Wien, Newton sur g(x)=(x-5)e^x+5")
g = lambda x: (x-5)*np.exp(x) + 5; dg = lambda x: (x-4)*np.exp(x)
xk = 5.0
for k in range(6):
    print(f"  x{k}={xk:.12f}  g={g(xk):.3e}"); xk = xk - g(xk)/dg(xk)

print("=== Fin 11 : pendule, moindres carrés T/T0 = a + b Y0^2")
Y0 = np.array([np.pi/18, np.pi/9, np.pi/6, 2*np.pi/9]); TT = np.array([1.002, 1.008, 1.020, 1.030])
G = np.c_[np.ones(4), Y0**2]; A = G.T@G; bb = G.T@TT; ab = np.linalg.solve(A, bb)
print(" A =", A.tolist(), " b =", bb.tolist(), " (a, b) =", ab)
ex = 2/np.pi*ellipk(np.sin(Y0/2)**2)
print(" valeurs exactes (2/pi)K(sin^2(Y0/2)) =", ex)
G2 = G; ab2 = np.linalg.solve(G2.T@G2, G2.T@ex); print(" (a,b) avec valeurs exactes =", ab2)
# trapèzes composites pour T(Y0)/T0
for y0 in Y0:
    A_ = np.sin(y0/2)**2
    for n in (4, 8):
        xs = np.linspace(0, np.pi/2, n+1); f = 1/np.sqrt(1 - A_*np.sin(xs)**2)
        T_ = (np.pi/2/n)*(f[0]/2 + f[1:-1].sum() + f[-1]/2)
        print(f"  Y0={y0:.4f} n={n}: T/T0 = {2/np.pi*T_:.8f}")

print("=== Fin 12 : capacité calorifique")
T = np.array([0.1, 0.3, 0.5]); Cv = np.array([0.211, 0.694, 1.361])
Q = 0.2/3*(Cv[0] + 4*Cv[1] + Cv[2]); print(f" Q = {Q:.6f} x 1e-4 cal/mol")
G = np.c_[np.ones(3), T]; A = G.T@G; bb = G.T@Cv; ab = np.linalg.solve(A, bb)
print(" système :", A.tolist(), bb.tolist(), " (a,b) =", ab, "(x1e-4)")
print(" Q avec la droite =", ab[0]*0.4 + ab[1]*(0.5**2-0.1**2)/2)

print("=== Fin 5 : test du tir avec c = 0 (solution exacte)")
def tir(a, b, c, gam, n=200):
    # y'' = -(a/x) y' - b e^{c y}, y(0)=gam, y'(0)=0 ; en x=0 : y'' = -b e^{c gam}/(1+a)
    def F(x, Y):
        y, yp = Y
        ypp = -b*np.exp(c*y)/(1+a) if x == 0 else -(a/x)*yp - b*np.exp(c*y)
        return np.array([yp, ypp])
    h = 1/n; x = 0.0; Y = np.array([gam, 0.0])
    for i in range(n):
        k1 = h*F(x, Y); k2 = h*F(x+h/2, Y+k1/2); k3 = h*F(x+h/2, Y+k2/2); k4 = h*F(x+h, Y+k3)
        Y = Y + (k1 + 2*k2 + 2*k3 + k4)/6; x += h
    return Y[0]
a_, b_ = 1.0, 1.0
g0, g1 = 1.0, 2.0
for it in range(6):
    y0_, y1_ = tir(a_, b_, 0.0, g0), tir(a_, b_, 0.0, g1)
    g0, g1 = g1, g1 - (g1-g0)*(y1_-1)/(y1_-y0_)
    if abs(g1-g0) < 1e-12: break
print(f" gamma trouvé = {g1:.10f}  exact 1 + b/(2(1+a)) = {1 + b_/(2*(1+a_)):.10f}")
a_, b_, c_ = 1.0, 1.0, 1.0
g0, g1 = 0.0, 0.5
for it in range(30):
    y0_, y1_ = tir(a_, b_, c_, g0), tir(a_, b_, c_, g1)
    g0, g1 = g1, g1 - (g1-g0)*(y1_-1)/(y1_-y0_)
    if abs(g1-g0) < 1e-12: break
print(f" exemple a=b=c=1 : y(0) = {g1:.10f}, y(1) = {tir(a_, b_, c_, g1):.12f}")

print("=== Fin 6 : Child-Langmuir V = (9a/4)^(2/3) x^(4/3)")
a_ = 1.0; C = (9*a_/4)**(2/3)
xx = 0.7; V = C*xx**(4/3); Vpp = C*(4/3)*(1/3)*xx**(-2/3)
print(" V'' =", Vpp, " a V^-1/2 =", a_*V**-0.5)

print("=== V, ex. 6 : équation intégro-différentielle, test k=1, f=1-2e^x, y=e^x")
def integro(N, xmax=1.0):
    h = xmax/N; xs = np.linspace(0, xmax, N+1); y = np.zeros(N+1); y[0] = 1.0
    k = lambda x, t: 1.0; f = lambda x: 1 - 2*np.exp(x)
    def F(n):
        I = 0.0
        if n > 0:
            I = h*(k(xs[n], xs[0])*y[0]/2 + sum(k(xs[n], xs[j])*y[j] for j in range(1, n)) + k(xs[n], xs[n])*y[n]/2)
        return 2*y[n] + f(xs[n]) + I
    for n in range(N): y[n+1] = y[n] + h*F(n)
    return y[-1]
for N in (10, 20, 40, 80): print(f"  N={N}: y(1)={integro(N):.6f} erreur={integro(N)-np.e:.2e}")

print("=== Fin 13 : équation intégrale du cours par trapèzes")
for n in (3, 5, 9, 17):
    tt = np.linspace(0, 1, n); hh = 1/(n-1); c = np.full(n, 1/(n-1)); c[0] = c[-1] = 1/(2*(n-1))
    K = lambda x, t: (x + t)/3
    A = np.eye(n) - np.array([[c[j]*K(tt[i], tt[j]) for j in range(n)] for i in range(n)])
    f = 5*tt/6 - 1/9; phi = np.linalg.solve(A, f)
    print(f"  n={n}: max|phi - x| = {np.abs(phi - tt).max():.3e}")
