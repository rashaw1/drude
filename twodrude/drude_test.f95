program drude_test
  !
  ! PURPOSE: To test the module drude_wf
  !

  use drude_wf
  implicit none

  real(dbl), dimension(3) :: r1, r2, R ! Positions 
  real(dbl) :: wf, laplacian, param, pmu, pq ! Wavefunction and Laplacian, and parameters
  real(dbl), dimension(3) :: grad1_wf, grad2_wf ! Gradients
  
  ! Set parameters
  param = 1d0
  pmu = 2d0
  pq = 1.5d0
  call init_drude(param, param, param, param, pmu, pq)

  ! Initialise positions
  r1(1) = 0.5d0
  r1(2) = -0.1d0
  r1(3) = 1.0d0
  r2(1) = 0.5d0
  r2(2) = 0.5d0
  r2(3) = -0.5d0
  R(1) = 1.0d0
  R(2) = 1.0d0
  R(3) = 0.0d0

  ! Calculate values and print out
  call calc_wf_derivatives(r1, r2, R, wf, grad1_wf, grad2_wf, laplacian)

  write(*, *) 'wf = ', wf
  write(*, *) 'grad1 = ', grad1_wf
  write(*, *) 'grad2 = ', grad2_wf
  write(*, *) 'laplacian = ', laplacian

end program drude_test
