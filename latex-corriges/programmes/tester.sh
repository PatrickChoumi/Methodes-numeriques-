#!/bin/sh
# Exécute tous les programmes Python des corrigés avec les données utilisées
# dans le document (Python 3 requis, aucune bibliothèque externe).
set -e
cd "$(dirname "$0")"
for f in lagrange pc_ab2 desantis ex4_pc ex6_rk2 ex8_jacobi gauss ex11_pendule \
         ex11_trapeze ex12_mc ex14_explicite ex14_richardson ex16_onde; do
  echo "=== $f.py"
  python3 "$f.py"
done
