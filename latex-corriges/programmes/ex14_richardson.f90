! Exercice 14, question 2 : schema de Richardson (saute-mouton)
!   T(i,j+1) = T(i,j-1) + 2 lambda (T(i+1,j) - 2 T(i,j) + T(i-1,j)),
! demarre avec T(i,1) donne par le schema explicite de la question 1.
program richardson
  implicit none
  integer, parameter :: mmax = 1000
  double precision, parameter :: pi = 3.14159265358979324d0
  double precision :: uold(0:mmax), u(0:mmax), unew(0:mmax)
  double precision :: l, alpha, t0, t1, tl, dt, tmax, dx, lambda, t
  integer :: m, i, j, nt

  write(*,*) 'L, alpha, T0, T1, TL, m, dt, tmax :'
  read(*,*) l, alpha, t0, t1, tl, m, dt, tmax
  dx = l / m
  lambda = alpha*dt/dx**2
  write(*,'(a, f8.4)') ' lambda = ', lambda
  ! niveau j = 0 : condition initiale
  uold(0:m-1) = t0
  uold(m) = tl
  ! niveau j = 1 : schema explicite (question 1)
  do i = 1, m - 1
     u(i) = (1d0 - 2d0*lambda)*uold(i) + lambda*(uold(i+1) + uold(i-1))
  end do
  u(0) = t0 + t1*sin(pi*dt/2d0)
  u(m) = tl
  nt = nint(tmax/dt)
  ! niveaux j = 2, 3, ... : Richardson
  do j = 1, nt - 1
     t = (j + 1)*dt
     do i = 1, m - 1
        unew(i) = uold(i) + 2d0*lambda*(u(i+1) - 2d0*u(i) + u(i-1))
     end do
     unew(0) = t0 + t1*sin(pi*t/2d0)
     unew(m) = tl
     uold(0:m) = u(0:m)
     u(0:m) = unew(0:m)
     if (mod(j + 1, nt/10) == 0) write(*,'(a, f8.3, a, es12.4)') &
          ' t = ', t, '   T(L/2) = ', u(m/2)
  end do
end program richardson
