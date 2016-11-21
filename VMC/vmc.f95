!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           !
!  VARIATIONAL MONTE CARLO FOR DRUDE OSCILLATORS            !
!  Robert Shaw, 2016                                        !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program vmc

  use constants
  use random
  use vmc_input
  use n_drude_wf
  use omp_lib

  implicit none
  
  integer                                    :: nwalk ! Number of walkers
  integer                                    :: ndrude ! Number of oscillators
  integer                                    :: nsteps, nequil
  real(dbl), allocatable, dimension(:, :)    :: rk, rknew, R ! Positions
  real(dbl), allocatable, dimension(:, :, :) :: walkers ! Walker positions
  real(dbl), allocatable, dimension(:)       :: p_omega, p_q, p_mu, p_alpha ! Frew, charge, mass, parameter
  real(dbl), allocatable, dimension(:, :)    :: p_gamma, p_eta, p_epsilon ! Parameters
  real(dbl)                                  :: dt ! Step size
  real(dbl)                                  :: e_sum, e_sq_sum ! Accumulators
  real(dbl)                                  :: rand, psiold, psinew
  real(dbl)                                  :: eold, enew, ekinloc, rho, enew_sum, enewsq_sum
  real(dbl), allocatable, dimension(:, :)    :: grads
  real(dbl)                                  :: start_time, end_time

  integer                                    :: iwalker, j, k, istep, icount
  character(len=10)                          :: time, date
  
  ! Read in parameters
  call vmc_read_input(nwalk, ndrude, nsteps, nequil, dt)

  ! Allocate position arrays
  allocate(rk(ndrude, 3))
  allocate(rknew(ndrude, 3))
  allocate(R(ndrude, 3))
  allocate(walkers(nwalk, ndrude, 3))
  allocate(grads(ndrude, 3))
  allocate(p_gamma(ndrude, ndrude))
  allocate(p_eta(ndrude, ndrude))
  allocate(p_epsilon(ndrude, ndrude))
  allocate(p_alpha(ndrude))
  allocate(p_omega(ndrude))
  allocate(p_q(ndrude))
  allocate(p_mu(ndrude))

  ! Read geometry and parameters
  call vmc_read_geom(ndrude, R, p_omega, p_q, p_mu)
  call vmc_read_params(ndrude, p_alpha, p_gamma, p_eta, p_epsilon)

  ! Initialise n_drude_wf
  call init_drude(ndrude, p_alpha, p_mu, p_q, p_gamma, p_eta, p_epsilon)

  ! Initialise random number generator
  call random_init()
  
  ! Zero accumulators
  e_sum = 0d0
  e_sq_sum = 0d0

  ! Initialise walker positions
  do iwalker = 1, nwalk
     do j = 1, ndrude
        do k = 1, 3
           rand = random_uniform(0d0, 1d0)
           walkers(iwalker, j, k) = R(j, k) + rand
        enddo
     enddo
  enddo

  call date_and_time(time=time, date=date)

  write(*, *) 'VARIATIONAL MONTE CARLO FOR DRUDE OSCILLATORS'
  write(*, *) 'Robert Shaw 2016'
  write(*, *) 'Date: ', date(7:8), '/', date(5:6), '/', date(1:4)
  write(*, *) 'Time: ', time(1:2), ':', time(3:4)
  write(*, *) ''
  write(*, *) 'Beginning equilibration:'
  write(*, *) 'nequil =', nequil
  write(*, *) ''

10 format(1X, I8, F15.8, F15.8)
  write(*, '(1X, A8, A15, A15)') 'Step', 'E Step', 'E Avg'
  icount = 1
  start_time = OMP_GET_WTIME()
  do istep = 1, nequil + nsteps
     enew_sum = 0d0
     enewsq_sum = 0d0
     !$OMP parallel default(none), shared(R, ndrude, p_mu, p_q, p_omega, walkers), &
     !$OMP private(rk, psiold, grads, ekinloc, eold, rknew, psinew, rand, rho, enew, iwalker), &
     !$OMP reduction(+: enewsq_sum, e_sq_sum)
     
     !$OMP do
     do iwalker = 1, nwalk
        rk(:, :) = walkers(iwalker, :, :)
        call calc_wf_derivatives(rk, R, psiold, grads, ekinloc)
        call calc_energy(ndrude, rk, R, p_mu, p_q, p_omega, ekinloc, psiold, eold)
        
        call move(ndrude, rk, rknew, dt)
        call calc_wf_derivatives(rknew, R, psinew, grads, ekinloc)

        rand = random_uniform(0d0, 1d0)
        rho = psinew**2 / psiold**2
        enew = eold
        if (rho > rand) then
           call calc_energy(ndrude, rknew, R, p_mu, p_q, p_omega, ekinloc, psinew, enew)
           walkers(iwalker, :, :) = rknew(:, :)
           psiold = psinew
        end if

        ! Accumulate
        enew_sum = enew_sum + enew
        enewsq_sum = enewsq_sum + enew**2

        eold = enew

     enddo
     !$OMP end do
     !$OMP end parallel
     e_sum = e_sum + enew_sum
     e_sq_sum = e_sq_sum + enewsq_sum
     write(*, 10) istep, enew_sum/dble(nwalk), e_sum/(dble(icount*nwalk))
     
     icount = icount + 1

     if ( istep == nequil ) then
        ! Move on to main run
        write(*, *) 'Finished equilibration, beginning main run:'
        write(*, *) 'nsteps =', nsteps
        write(*, *) ''
        e_sum = 0d0
        e_sq_sum = 0d0
        icount = 1
     endif
        
  enddo
  end_time = OMP_GET_WTIME()

  call fin_drude()
  
  ! Write output
  write(*, *) ''
  write(*, *) 'Run finished.'
  e_sum = e_sum/dble(nsteps*nwalk)
  e_sq_sum = e_sq_sum/dble(nsteps*nwalk)
  write(*, '(1X, A15, F15.8)') 'Average energy: ', e_sum
  write(*, '(1X, A15, F15.8)') 'Variance: ', e_sq_sum - e_sum**2
  write(*, '(1X, A15, F15.8)') 'Std. dev.: ', sqrt(e_sq_sum-e_sum**2)
  write(*, '(1X, A15, F15.8)') 'Time elapsed: ', end_time - start_time

contains

  subroutine calc_energy(nd, rk, R, mu, q, omega, ekinloc, psi, e)
    ! Calculate the local energy
    implicit none

    integer,                     intent(in)  :: nd
    real(dbl), dimension(nd, 3), intent(in)  :: rk, R
    real(dbl), dimension(nd),    intent(in)  :: mu, q, omega
    real(dbl),                   intent(in)  :: ekinloc, psi
    real(dbl),                   intent(out) :: e

    integer                                  :: i, j
    real(dbl), dimension(3)                  :: rij
    real(dbl)                                :: drij
    
    e = 0d0
    ! Drude potential
    do i = 1, nd
       rij = rk(i, :) - R(i, :)
       e = e + 0.5d0 * mu(i) * omega(i)**2 * dot_product(rij, rij)
    enddo

    ! Coulomb attraction
    do i = 1, nd
       do j = 1, nd
          if (i == j) cycle
          rij = rk(j, :) - R(i, :)
          drij = sqrt(dot_product(rij, rij))
          e = e - q(i) * q(j) / drij
       enddo
    enddo

    ! Coulomb repulsion
    do i = 1, nd
       do j = 1, i-1
          rij = rk(j, :) - rk(i, :)
          drij = sqrt(dot_product(rij, rij))
          e = e + q(i) * q(j) / drij

          rij = R(j, :) - R(i, :)
          drij = sqrt(dot_product(rij, rij))
          e = e + q(i) * q(j) / drij
       enddo
    enddo      
    
    e = e - 0.5d0 * ekinloc
  end subroutine calc_energy

  subroutine move(nd, rk, rknew, dt)
    ! Make a metropolis move
    implicit none
    
    integer,                     intent(in)  :: nd
    real(dbl),                   intent(in)  :: dt
    real(dbl), dimension(nd, 3), intent(in)  :: rk
    real(dbl), dimension(nd, 3), intent(out) :: rknew

    real(dbl)                                :: rand
    integer                                  :: i, j
    
    rknew(:, :) = rk(:, :)
    do i = 1, nd
       do j = 1, 3
          rand = random_normal(0d0, 1d0)
          rknew(i, j) = rknew(i, j) + dt * rand
       enddo
    enddo
  end subroutine move

end program vmc

    
          
