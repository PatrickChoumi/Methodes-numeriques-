{ Exercice 11, question 2 : T(Y0)/T0 = (2/pi) * integrale de 0 a pi/2
  de dx / sqrt(1 - A sin^2 x), A = sin^2(Y0/2), par les trapezes composites. }
program trapeze;
var
  y0, A, h, x, S : real;
  n, k : integer;

function f(x : real) : real;
begin
  f := 1 / sqrt(1 - A * sqr(sin(x)))
end;

begin
  write('Y0 (radians), n : ');
  readln(y0, n);
  A := sqr(sin(y0 / 2));
  h := (pi / 2) / n;
  S := (f(0) + f(pi / 2)) / 2;
  x := 0;
  for k := 1 to n - 1 do
  begin
    x := x + h;
    S := S + f(x)
  end;
  S := h * S;
  writeln('T/T0 = ', 2 / pi * S : 12 : 8)
end.
