{ Exercice 4 : y' = -a y - b y^4 + c, y(0) = y0.
  Predicteur : Adams-Bashforth d'ordre 2 ; correcteur : RK2 (trapezes). }
program ex4pc;
var
  a, b, c, y0, xmax, h, x, y, fk, fprec, p : real;
  n, k : integer;

function f(x, y : real) : real;
begin
  f := -a * y - b * sqr(sqr(y)) + c
end;

begin
  write('a, b, c, y0, xmax, n : ');
  readln(a, b, c, y0, xmax, n);
  h := xmax / n;
  x := 0.0;
  y := y0;
  writeln(x:8:4, y:16:10);
  { Demarrage : y1 par RK2 }
  fprec := f(x, y);
  y := y + h / 2 * (fprec + f(x + h, y + h * fprec));
  x := x + h;
  writeln(x:8:4, y:16:10);
  for k := 1 to n - 1 do
  begin
    fk := f(x, y);
    p := y + h / 2 * (3 * fk - fprec);      { predicteur AB2 }
    y := y + h / 2 * (fk + f(x + h, p));    { correcteur RK2 }
    fprec := fk;
    x := x + h;
    writeln(x:8:4, y:16:10)
  end
end.
