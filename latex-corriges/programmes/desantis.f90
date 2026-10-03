! Exercice 3 : facteur de compressibilite de De Santis
!   z = (1 + y + y^2 - y^3) / (1 - y)^3 ,  y = b / (4 v)
! On cherche y (puis b) par Newton-Raphson sur
!   f(y) = z (1 - y)^3 - (1 + y + y^2 - y^3) = 0
program desantis
  implicit none
  double precision :: z, v, y, dy, eps, f, fp
  integer :: k, kmax

  write(*,*) 'z, v, y0 (approximation initiale), precision, nb max iterations :'
  read(*,*) z, v, y, eps, kmax

  do k = 1, kmax
     f  = z*(1d0 - y)**3 - (1d0 + y + y**2 - y**3)
     fp = -3d0*z*(1d0 - y)**2 - (1d0 + 2d0*y - 3d0*y**2)
     dy = f / fp
     y  = y - dy
     write(*,'(i4, f16.12)') k, y
     if (abs(dy) < eps) exit
  end do

  if (abs(dy) >= eps) then
     write(*,*) 'Pas de convergence en ', kmax, ' iterations'
  else
     write(*,'(a, f14.10, a, f14.10)') ' y = ', y, '   b = 4 v y = ', 4d0*v*y
  end if
end program desantis
