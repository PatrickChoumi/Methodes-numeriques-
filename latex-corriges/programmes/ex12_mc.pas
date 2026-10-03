{ Exercice 12 : moindres carres Cv = a + b T (equations normales)
  puis chaleur Q par Simpson composite sur les points experimentaux
  (abscisses equidistantes, nombre de points impair). }
program moindrescarres;
const
  MMAX = 1000;
var
  T, C : array[1..MMAX] of real;
  m, i : integer;
  s0, s1, s2, sy, sty, det, a, b, h, Q : real;

begin
  write('Nombre de mesures m : ');
  readln(m);
  for i := 1 to m do readln(T[i], C[i]);
  { 1) coefficients du systeme normal }
  s0 := m; s1 := 0; s2 := 0; sy := 0; sty := 0;
  for i := 1 to m do
  begin
    s1 := s1 + T[i];
    s2 := s2 + sqr(T[i]);
    sy := sy + C[i];
    sty := sty + T[i] * C[i]
  end;
  { 2) resolution du systeme 2 x 2 (Cramer) }
  det := s0 * s2 - sqr(s1);
  a := (s2 * sy - s1 * sty) / det;
  b := (s0 * sty - s1 * sy) / det;
  writeln('Systeme : ', s0:8:4, s1:10:4, ' | ', sy:10:4);
  writeln('          ', s1:8:4, s2:10:4, ' | ', sty:10:4);
  writeln('a = ', a:12:6, '   b = ', b:12:6);
  { 3) Simpson composite : m impair, pas constant }
  if odd(m) then
  begin
    h := T[2] - T[1];
    Q := C[1] + C[m];
    for i := 2 to m - 1 do
      if odd(i) then Q := Q + 2 * C[i]
      else Q := Q + 4 * C[i];
    Q := h / 3 * Q;
    writeln('Q (Simpson) = ', Q:12:6)
  end
  else
    writeln('Simpson composite : il faut un nombre impair de points')
end.
