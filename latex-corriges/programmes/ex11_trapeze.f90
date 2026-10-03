! Exercice 11, question 2 : T(Y0)/T0 = (2/pi) * integrale de 0 a pi/2
! de dx / sqrt(1 - A sin^2 x), A = sin^2(Y0/2), par les trapezes composites.
program trapeze
  implicit none
  double precision, parameter :: pi = 3.14159265358979324d0
  double precision :: y0, a, h, x, s
  integer :: n, k

  write(*,*) 'Y0 (radians), n :'
  read(*,*) y0, n
  a = sin(y0/2d0)**2
  h = (pi/2d0) / n
  s = (f(0d0) + f(pi/2d0)) / 2d0
  x = 0d0
  do k = 1, n - 1
     x = x + h
     s = s + f(x)
  end do
  s = h * s
  write(*,'(a, f12.8)') ' T/T0 = ', 2d0/pi*s
contains
  double precision function f(x)
    double precision, intent(in) :: x
    f = 1d0 / sqrt(1d0 - a*sin(x)**2)
  end function f
end program trapeze
