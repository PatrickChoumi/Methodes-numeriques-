{ Exercice 6, question 1 : V'' = a V^(-1/2), V(0) = V'(0) = 0.
  Systeme : V' = W, W' = a / sqrt(V), resolu par RK2 (trapezes).
  Le second membre n'est pas defini en V = 0 : on part de x0 = eps
  avec la solution exacte V = C x^(4/3), C = (9a/4)^(2/3). }
program ex6rk2;
var
  a, l, eps, h, x, v, w, c, k1v, k1w, k2v, k2w : real;
  n, i : integer;

function g(v : real) : real;           { W' = g(V) }
begin
  g := a / sqrt(v)
end;

begin
  write('a, l, eps, n : ');
  readln(a, l, eps, n);
  c := exp(2.0 / 3.0 * ln(9.0 * a / 4.0));
  x := eps;
  v := c * exp(4.0 / 3.0 * ln(eps));
  w := 4.0 / 3.0 * c * exp(1.0 / 3.0 * ln(eps));
  h := (l - eps) / n;
  for i := 1 to n do
  begin
    k1v := h * w;          k1w := h * g(v);
    k2v := h * (w + k1w);  k2w := h * g(v + k1v);
    v := v + (k1v + k2v) / 2;
    w := w + (k1w + k2w) / 2;
    x := x + h
  end;
  writeln('x = ', x:8:4, '  V = ', v:14:8,
          '  exact = ', c * exp(4.0 / 3.0 * ln(x)):14:8)
end.
