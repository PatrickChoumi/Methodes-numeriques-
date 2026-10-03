"""Exercice 12 : moindres carrés Cv = a + b T (équations normales),
puis chaleur Q par Simpson composite sur les points expérimentaux
(abscisses équidistantes, nombre de points impair)."""


def moindres_carres(T, C):
    m = len(T)
    s1 = s2 = sy = sty = 0.0
    for t, c in zip(T, C):         # 1) coefficients du système normal
        s1 += t
        s2 += t * t
        sy += c
        sty += t * c
    print(f"Système : {m:8.4f} {s1:10.4f} | {sy:10.4f}")
    print(f"          {s1:8.4f} {s2:10.4f} | {sty:10.4f}")
    det = m * s2 - s1 * s1         # 2) résolution 2 x 2 (Cramer)
    return (s2 * sy - s1 * sty) / det, (m * sty - s1 * sy) / det


def simpson_tableau(T, C):
    m = len(C)
    if m % 2 == 0:
        raise ValueError("il faut un nombre impair de points")
    h = T[1] - T[0]
    q = C[0] + C[m - 1]
    for i in range(1, m - 1):      # indices impairs : poids 4
        q += 4 * C[i] if i % 2 == 1 else 2 * C[i]
    return h / 3 * q


if __name__ == "__main__":
    T = [0.1, 0.3, 0.5]
    C = [0.211, 0.694, 1.361]      # Cv en 1e-4 cal/(K.mol)
    a, b = moindres_carres(T, C)
    print(f"a = {a:.6f}   b = {b:.6f}")
    print(f"Q (Simpson) = {simpson_tableau(T, C):.6f}  (x 1e-4 cal/mol)")
    q = a * (T[-1] - T[0]) + b * (T[-1] ** 2 - T[0] ** 2) / 2
    print(f"Q (droite)  = {q:.6f}")
