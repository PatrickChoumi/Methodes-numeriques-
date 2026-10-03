#!/bin/sh
# Compile et execute tous les programmes des corriges avec les donnees de test
# utilisees dans le document (gcc, gfortran et fpc requis).
set -e
cd "$(dirname "$0")"
mkdir -p bin
gcc -Wall -o bin/lagrange lagrange.c
gcc -Wall -o bin/pc_ab2 pc_ab2.c -lm
gcc -Wall -o bin/ex16_onde ex16_onde.c -lm
for f in desantis ex4_pc ex6_rk2 gauss ex11_pendule ex11_trapeze ex14_richardson; do
  gfortran -Wall -o bin/${f}_f $f.f90
done
for f in ex4_pc ex6_rk2 ex8_jacobi gauss ex11_trapeze ex12_mc ex14_explicite; do
  fpc -v0 -FUbin -obin/${f}_p $f.pas > /dev/null
done

echo "--- III, ex. 3 : Lagrange"
printf "5\n93.0 11.38\n96.2 12.80\n100.0 14.70\n104.2 17.07\n108.7 19.91\n102\n" | bin/lagrange | tail -1
printf "5\n11.38 93.0\n12.80 96.2\n14.70 100.0\n17.07 104.2\n19.91 108.7\n13.5\n" | bin/lagrange | tail -1
echo "--- V : predicteur-correcteur AB2/RK2 sur y' = -2xy"
bin/pc_ab2 | tail -1
echo "--- Ex. 3 : De Santis (z = 3.9737609, v = 1, y0 = 0.5)"
echo "3.9737609 1 0.5 1e-12 50" | bin/desantis_f | tail -1
echo "--- Ex. 4 : y' = -y - 0.5 y^4 + 2, y(0) = 0, x = 2, n = 40"
echo "1 0.5 2 0 2 40" | bin/ex4_pc_p | tail -1
echo "1 0.5 2 0 2 40" | bin/ex4_pc_f | tail -1
echo "--- Ex. 6 : V'' = a V^(-1/2), a = 1, l = 1, eps = 1e-2"
for n in 100 200 400; do echo "1 1 1e-2 $n" | bin/ex6_rk2_p | tail -1; done
echo "1 1 1e-2 400" | bin/ex6_rk2_f | tail -1
echo "--- Ex. 8 : differences finies + Jacobi (exemple de Burden-Faires)"
echo "1 2 1 2 10 1e-12 10000" | bin/ex8_jacobi_p
echo "--- Ex. 10 et 13 : Gauss avec pivot partiel (a11 = 0)"
printf "3\n0 2 1 5\n1 1 1 6\n2 1 3 13\n" | bin/gauss_p
printf "3\n0 2 1 5\n1 1 1 6\n2 1 3 13\n" | bin/gauss_f
echo "--- Ex. 11 : pendule RK2 (w0 = 1, Y0 = 2 pi/9)"
echo "1 0.6981317007977318 0 4.0 4000 1000" | bin/ex11_pendule_f | tail -4
echo "--- Ex. 11 : periode par les trapezes (n = 4)"
for y in 0.17453292519943295 0.3490658503988659 0.5235987755982988 0.6981317007977318; do
  echo "$y 4" | bin/ex11_trapeze_p; echo
  echo "$y 4" | bin/ex11_trapeze_f | tail -1
done
echo "--- Ex. 12 : moindres carres et Simpson"
printf "3\n0.1 0.211\n0.3 0.694\n0.5 1.361\n" | bin/ex12_mc_p
echo "--- Ex. 14 : schema explicite (L = 1, alpha = 0.1, T0 = 20, T1 = 10, TL = 30)"
echo "1 0.1 20 10 30 10 0.04 3" | bin/ex14_explicite_p | sed -n '1p;3,14p'
echo "--- Ex. 14 : Richardson (lambda = 0.1) : divergence"
echo "1 0.1 20 10 30 10 0.01 3" | bin/ex14_richardson_f | sed -n '3p;6p;12p'
echo "--- Ex. 16 : onde amortie non lineaire"
echo "0 1 0 1 50 0.01 1.5" | bin/ex16_onde | grep '^  0.5000'
echo "0.4 1 0 1 50 0.01 1.5" | bin/ex16_onde | grep '^  0.5000'
for m in 50 100 200; do k=$(echo "scale=8; 0.5/$m" | bc); echo "0.4 1 -2 1 $m $k 1.5" | bin/ex16_onde | grep '^  0.5000'; done
