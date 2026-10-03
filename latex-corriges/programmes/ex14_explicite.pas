{ Exercice 14, question 1 : mur d'epaisseur L, dT/dt = alpha d2T/dx2,
  T(x,0) = T0, T(0,t) = T0 + T1 sin(pi t / 2), T(L,t) = TL.
  Schema explicite (differences progressives en temps). }
program chaleurexplicite;
const
  MMAX = 1000;
var
  U, V : array[0..MMAX] of real;
  L, alpha, T0, T1, TL, dt, tmax, dx, lambda, t : real;
  m, i, j, nt : integer;

begin
  write('L, alpha, T0, T1, TL, m, dt, tmax : ');
  readln(L, alpha, T0, T1, TL, m, dt, tmax);
  dx := L / m;
  lambda := alpha * dt / sqr(dx);
  writeln('lambda = ', lambda:8:4);
  if lambda > 0.5 then writeln('Attention : lambda > 1/2, schema instable');
  { condition initiale (la face x = L est a TL des t > 0) }
  for i := 0 to m - 1 do U[i] := T0;
  U[m] := TL;
  nt := round(tmax / dt);
  t := 0;
  for j := 1 to nt do
  begin
    for i := 1 to m - 1 do
      V[i] := (1 - 2 * lambda) * U[i] + lambda * (U[i + 1] + U[i - 1]);
    t := j * dt;
    V[0] := T0 + T1 * sin(pi * t / 2);
    V[m] := TL;
    for i := 0 to m do U[i] := V[i]
  end;
  writeln('Temperatures a t = ', t:8:4);
  for i := 0 to m do writeln(i * dx:10:4, U[i]:14:6)
end.
