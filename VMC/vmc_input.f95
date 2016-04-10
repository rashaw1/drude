  ! Input module for VMC

module vmc_input

  use constants
  
  implicit none
  
contains
  subroutine vmc_read_input(N, ndrude, nsteps, nequil, dt)
    ! Read the input parameters from control
    integer,   intent(out) :: N, ndrude, nsteps, nequil
    real(dbl), intent(out) :: dt

    integer                :: filestream, ioerr
    character(20)          :: keyword, str1
    
    open(newunit = filestream, file='control', action="read", iostat = ioerr)
    if (ioerr > 0) write(*, *) 'Failed to open control.'

    str1 = ""
    do
       read(filestream, *, iostat = ioerr) keyword
       if ( ioerr < 0 ) then
          exit ! EOF
       else if (ioerr > 0) then
          write(*, *) 'Could not read control file.'
       end if

       backspace(filestream)
       select case (keyword)
       case ("nwalkers")
          read(filestream, *) keyword, N
       case ("ndrude")
          read(filestream, *) keyword, ndrude
       case ("nsteps")
          read(filestream, *) keyword, nsteps
       case ("nequil")
          read(filestream, *) keyword, nequil
       case ("timestep")
          read(filestream, *) keyword, dt
       case default
          read(filestream, '(A)') str1
          write(*, *) 'Could not handle control line ', str1
       end select
    end do

    close(filestream)
  end subroutine vmc_read_input

  subroutine vmc_read_geom(ndrude, R, omega, q, mu)
    ! Read in the geometry of the oscillators
    integer,                         intent(in)  :: ndrude
    real(dbl), dimension(ndrude, 3), intent(out) :: R
    real(dbl), dimension(ndrude),    intent(out) :: omega, q, mu

    integer                                      :: filestream, ioerr, i
       
    open(newunit = filestream, file='coord', action="read", iostat = ioerr)
    if (ioerr > 0) write(*, *) 'Failed to open geometry file.'

    do i = 1, ndrude
       read(filestream, *, iostat = ioerr) R(i, 1), R(i, 2), R(i, 3), omega(i), q(i), mu(i)
       if (ioerr < 0) then
          exit ! EOF
       else if (ioerr > 0) then
          write(*, *) 'Could not read geometry file.'
       end if
    enddo

    close(filestream)
  end subroutine vmc_read_geom

  subroutine vmc_read_params(ndrude, alpha, gamma, eta, epsilon)
    ! Read the trial wavefunction parameters

    integer, intent(in) :: ndrude
    real(dbl), dimension(ndrude), intent(out) :: alpha
    real(dbl), dimension(ndrude, ndrude), intent(out) :: gamma, eta, epsilon

    integer                                      :: filestream, ioerr, i, j

    open(newunit = filestream, file='psidat', action="read", iostat = ioerr)
    if (ioerr > 0) write(*, *) 'Failed to open psi data file.'

    do i = 1, ndrude
       do j = 1, ndrude
          read(filestream, *, iostat = ioerr) alpha(i), gamma(i, j), eta(i, j), epsilon(i, j)
          if (ioerr < 0) then
             exit ! EOF
          else if (ioerr > 0) then
             write(*, *) 'Could not read psidat file.'
          end if
       enddo
    enddo

    close (filestream)
  end subroutine vmc_read_params
  
end module vmc_input
