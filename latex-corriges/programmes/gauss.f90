! Exercice 13 : methode de Gauss avec pivot partiel
! (triangularisation) puis resolution du systeme triangulaire (remontee).
program gausspivot
  implicit none
  integer, parameter :: nmax = 50
  double precision :: a(nmax, nmax), b(nmax), x(nmax)
  integer :: n, i

  read(*,*) n
  do i = 1, n
     read(*,*) a(i, 1:n), b(i)
  end do
  call triangularise(n, a, b)
  call remontee(n, a, b, x)
  do i = 1, n
     write(*,'(a, i0, a, f16.10)') 'x', i, ' = ', x(i)
  end do

contains

  subroutine triangularise(n, a, b)
    integer, intent(in) :: n
    double precision, intent(inout) :: a(nmax, nmax), b(nmax)
    integer :: i, k, ip
    double precision :: m, t(nmax), tb
    do k = 1, n - 1
       ip = k - 1 + maxloc(abs(a(k:n, k)), dim=1)   ! pivot maximal
       if (a(ip, k) == 0d0) stop 'Matrice singuliere'
       if (ip /= k) then
          t(k:n) = a(k, k:n); a(k, k:n) = a(ip, k:n); a(ip, k:n) = t(k:n)
          tb = b(k); b(k) = b(ip); b(ip) = tb
       end if
       do i = k + 1, n
          m = a(i, k) / a(k, k)
          a(i, k+1:n) = a(i, k+1:n) - m * a(k, k+1:n)
          a(i, k) = 0d0
          b(i) = b(i) - m * b(k)
       end do
    end do
  end subroutine triangularise

  subroutine remontee(n, a, b, x)
    integer, intent(in) :: n
    double precision, intent(in) :: a(nmax, nmax), b(nmax)
    double precision, intent(out) :: x(nmax)
    integer :: k, l
    double precision :: s
    x(n) = b(n) / a(n, n)
    do k = n - 1, 1, -1
       s = b(k)
       do l = k + 1, n
          s = s - a(k, l) * x(l)
       end do
       x(k) = s / a(k, k)
    end do
  end subroutine remontee

end program gausspivot
