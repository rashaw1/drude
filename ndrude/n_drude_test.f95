program n_drude_test
  !
  ! PURPOSE: To test the module n_drude_wf
  !

  use n_drude_wf
  implicit none
 
  integer, parameter :: No = 4 ! Number of oscillators
  real(dbl), dimension(:, :), allocatable :: rk, R ! Positions 
  real(dbl) :: wf, laplacian, palpha, pmu, pq ! Wavefunction and Laplacian, and parameters
  real(dbl), dimension(:, :), allocatable :: grads ! Gradients
  real(dbl), dimension(:, :), allocatable :: params ! gamma, eta, epsilon 

  integer :: i, j ! Loop indices
  
  ! Allocate arrays
  allocate( rk(No, 3) )
  allocate( R(No, 3) )
  allocate( grads(No, 3) )
  allocate( params(No, No) )
  
  ! Set parameters
  do i = 1, No
     do j = 1, No
        params(i, j) = 0.6d0
     end do
  end do
  pmu = 1d0
  pq = 0.5d0
  palpha = 1.0d0

  ! Initialise module
  call init_drude(No, palpha, pmu, pq, params, params, params)

  ! Initialise positions
  rk(1, 1) = 0.5d0
  rk(1, 2) = 0.5d0
  rk(1, 3) = 0.5d0
  R(1, 1) = -0.2d0
  R(1, 2) = 1.0d0
  R(1, 3) = 0.8d0
  rk(2, 1) = -0.1d0
  rk(2, 2) = 1.0d0
  rk(2, 3) = 0.5d0
  R(2, 1) = 0.7d0
  R(2, 2) = -0.3d0
  R(2, 3) = -0.6d0
  rk(3, 1) = 0.4d0
  rk(3, 2) = 0.5d0
  rk(3, 3) = -0.5d0
  R(3, 1) = 1.0d0
  R(3, 2) = 0.0d0
  R(3, 3) = 1.0d0
  rk(4, 1) = 0.5d0
  rk(4, 2) = -0.1d0
  rk(4, 3) = 0.6d0

  ! Calculate values and print out
  call calc_wf_derivatives(rk, R, wf, grads, laplacian)

  ! Output
  write(*, *) 'wf = ', wf
  write(*, *) 'grads = '
  do i = 1, No
     write(*, *) grads(i, :)
  end do
  write(*, *) 'laplacian = ', laplacian

  ! Deallocate arrays
  call fin_drude()
  deallocate ( params )
  deallocate ( rk )
  deallocate ( R )
  deallocate ( grads )
end program n_drude_test
