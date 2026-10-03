! Exercice 6, question 1 : V'' = a V^(-1/2), V(0) = V'(0) = 0.
! Systeme : V' = W, W' = a / sqrt(V), resolu par RK2 (trapezes).
! Le second membre n'est pas defini en V = 0 : on part de x0 = eps
! avec la solution exacte V = C x^(4/3), C = (9a/4)^(2/3).
program ex6rk2
  implicit none
  double precision :: a, l, eps, h, x, v, w, c, k1v, k1w, k2v, k2w
  integer :: n, i

  write(*,*) 'a, l, eps, n :'
  read(*,*) a, l, eps, n
  c = (9d0*a/4d0)**(2d0/3d0)
  x = eps
  v = c * eps**(4d0/3d0)
  w = 4d0/3d0 * c * eps**(1d0/3d0)
  h = (l - eps) / n
  do i = 1, n
     k1v = h*w;          k1w = h*a/sqrt(v)
     k2v = h*(w + k1w);  k2w = h*a/sqrt(v + k1v)
     v = v + (k1v + k2v)/2d0
     w = w + (k1w + k2w)/2d0
     x = x + h
  end do
  write(*,'(a, f8.4, a, f14.8, a, f14.8)') ' x = ', x, '  V = ', v, &
       '  exact = ', c*x**(4d0/3d0)
end program ex6rk2
