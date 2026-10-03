"""Chapitre III : interpolation de Lagrange (algorithme du cours).
Exercice III.4 : y_p pour x_p = 102, puis x_p pour y_p = 13,5."""


def lagrange(x, y, xp):
    """Valeur en xp du polynôme qui passe par les points (x[k], y[k])."""
    yp = 0.0
    for k in range(len(x)):
        ys = 1.0                  # ys = f_k(xp)
        for i in range(len(x)):
            if i != k:
                ys = ys * (xp - x[i]) / (x[k] - x[i])
        yp = yp + y[k] * ys
    return yp


if __name__ == "__main__":
    x = [93.0, 96.2, 100.0, 104.2, 108.7]
    y = [11.38, 12.80, 14.70, 17.07, 19.91]
    print(f"x_p = 102    y_p = {lagrange(x, y, 102.0):.6f}")
    # interpolation inverse : on échange les rôles de x et de y
    print(f"y_p = 13,5   x_p = {lagrange(y, x, 13.5):.6f}")
