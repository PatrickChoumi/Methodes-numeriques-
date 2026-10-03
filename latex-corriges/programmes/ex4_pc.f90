! Exercice 4 : y' = -a y - b y^4 + c, y(0) = y0.
! Predicteur : Adams-Bashforth d'ordre 2 ; correcteur : RK2 (trapezes).
program ex4pc
  implicit none
  double precision :: a, b, c, y0, xmax, h, x, y, fk, fprec, p
  integer :: n, k

  write(*,*) 'a, b, c, y0, xmax, n :'
  read(*,*) a, b, c, y0, xmax, n
  h = xmax / n
  x = 0d0
  y = y0
  write(*,'(f8.4, f16.10)') x, y
  ! Demarrage : y1 par RK2
  fprec = f(x, y)
  y = y + h/2d0 * (fprec + f(x + h, y + h*fprec))
  x = x + h
  write(*,'(f8.4, f16.10)') x, y
  do k = 1, n - 1
     fk = f(x, y)
     p  = y + h/2d0 * (3d0*fk - fprec)      ! predicteur AB2
     y  = y + h/2d0 * (fk + f(x + h, p))    ! correcteur RK2
     fprec = fk
     x = x + h
     write(*,'(f8.4, f16.10)') x, y
  end do

contains

  double precision function f(x, y)
    double precision, intent(in) :: x, y
    f = -a*y - b*y**4 + c
  end function f

end program ex4pc
