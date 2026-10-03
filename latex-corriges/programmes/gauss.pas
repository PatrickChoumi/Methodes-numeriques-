{ Exercices 10 et 13 : methode de Gauss avec pivot partiel
  (triangularisation) puis resolution du systeme triangulaire (remontee). }
program gausspivot;
const
  NMAX = 50;
type
  matrice = array[1..NMAX, 1..NMAX] of real;
  vecteur = array[1..NMAX] of real;
var
  a : matrice;
  b, x : vecteur;
  n, i, j : integer;

procedure triangularise(n : integer; var a : matrice; var b : vecteur);
var
  i, j, k, ip : integer;
  m, t : real;
begin
  for k := 1 to n - 1 do
  begin
    { pivot maximal : ligne ip telle que |a[ip,k]| soit maximal }
    ip := k;
    for i := k + 1 to n do
      if abs(a[i, k]) > abs(a[ip, k]) then ip := i;
    if a[ip, k] = 0 then
    begin
      writeln('Matrice singuliere'); halt
    end;
    if ip <> k then
    begin
      for j := k to n do
      begin t := a[k, j]; a[k, j] := a[ip, j]; a[ip, j] := t end;
      t := b[k]; b[k] := b[ip]; b[ip] := t
    end;
    for i := k + 1 to n do
    begin
      m := a[i, k] / a[k, k];
      for j := k + 1 to n do a[i, j] := a[i, j] - m * a[k, j];
      a[i, k] := 0;
      b[i] := b[i] - m * b[k]
    end
  end
end;

procedure remontee(n : integer; var a : matrice; var b, x : vecteur);
var
  k, l : integer;
  s : real;
begin
  x[n] := b[n] / a[n, n];
  for k := n - 1 downto 1 do
  begin
    s := b[k];
    for l := k + 1 to n do s := s - a[k, l] * x[l];
    x[k] := s / a[k, k]
  end
end;

begin
  readln(n);
  for i := 1 to n do
  begin
    for j := 1 to n do read(a[i, j]);
    readln(b[i])
  end;
  triangularise(n, a, b);
  remontee(n, a, b, x);
  for i := 1 to n do writeln('x', i, ' = ', x[i]:16:10)
end.
