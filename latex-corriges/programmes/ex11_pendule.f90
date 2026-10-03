! Exercice 11, question 1 : Y'' + w0^2 sin Y = 0, Y(0) = Y0, Y'(0) = Z0.
! Systeme : Y' = Z, Z' = -w0^2 sin Y, resolu par RK2 (trapezes).
program pendule
  implicit none
  double precision :: w0, y0, z0, tmax, h, t, y, z, k1y, k1z, k2y, k2z
  integer :: n, i, npas

  write(*,*) 'w0, Y0, Z0, tmax, n, npas (affichage tous les npas pas) :'
  read(*,*) w0, y0, z0, tmax, n, npas
  h = tmax / n
  t = 0d0
  y = y0
  z = z0
  write(*,'(3f16.8)') t, y, z
  do i = 1, n
     k1y = h*z;            k1z = -h*w0**2*sin(y)
     k2y = h*(z + k1z);    k2z = -h*w0**2*sin(y + k1y)
     y = y + (k1y + k2y)/2d0
     z = z + (k1z + k2z)/2d0
     t = t + h
     if (mod(i, npas) == 0) write(*,'(3f16.8)') t, y, z
  end do
end program pendule
