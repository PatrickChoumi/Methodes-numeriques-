{ Exercice 8 : y'' = p(x) y' + q(x) y + g(x), y(a) = c, y(b) = d.
  Differences finies centrees (ordre 2) puis methode iterative de Jacobi.
  Exemple (Burden-Faires, ex. 1 du section 11.3) :
  p = -2/x, q = 2/x^2, g = sin(ln x)/x^2, 1 <= x <= 2, y(1) = 1, y(2) = 2. }
program ex8jacobi;
const
  NMAX = 1000;
var
  y, ynew, aa, bb, cc, dd : array[0..NMAX] of real;
  a, b, ya, yb, h, x, eps, num, den : real;
  n, i, k, kmax : integer;

function p(x : real) : real; begin p := -2 / x end;
function q(x : real) : real; begin q := 2 / sqr(x) end;
function g(x : real) : real; begin g := sin(ln(x)) / sqr(x) end;

begin
  write('a, b, y(a), y(b), n (sous-intervalles), eps, kmax : ');
  readln(a, b, ya, yb, n, eps, kmax);
  h := (b - a) / n;
  { Coefficients : bb_i y_(i-1) + aa_i y_i + cc_i y_(i+1) = dd_i }
  for i := 1 to n - 1 do
  begin
    x := a + i * h;
    bb[i] := 1 + h / 2 * p(x);
    aa[i] := -(2 + sqr(h) * q(x));
    cc[i] := 1 - h / 2 * p(x);
    dd[i] := sqr(h) * g(x)
  end;
  { Valeurs initiales : conditions aux limites et 0 a l'interieur }
  y[0] := ya; y[n] := yb; ynew[0] := ya; ynew[n] := yb;
  for i := 1 to n - 1 do y[i] := 0;
  k := 0;
  repeat
    k := k + 1;
    num := 0; den := 0;
    for i := 1 to n - 1 do
    begin
      ynew[i] := (dd[i] - bb[i] * y[i - 1] - cc[i] * y[i + 1]) / aa[i];
      num := num + abs(ynew[i] - y[i]);
      den := den + abs(ynew[i])
    end;
    for i := 1 to n - 1 do y[i] := ynew[i]
  until (num / den < eps) or (k >= kmax);
  writeln('Iterations : ', k);
  for i := 0 to n do
    writeln(a + i * h:6:2, y[i]:14:8)
end.
