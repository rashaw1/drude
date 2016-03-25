module n_drude_wf
  !
  ! PURPOSE: To calculate the wavefunction, gradient, and Laplacian
  !          for a system of N drude oscillators.
  !

  implicit none
  
  ! Declare parameters
  integer, parameter :: dbl = selected_real_kind(15, 307)
  integer :: N ! Number of oscillators
  real(dbl) :: alpha, kappa ! Parameters
  real(dbl),  dimension(:, :), allocatable :: gamma, eta, epsilon, delta 
  
contains

  subroutine init_drude(p_N, p_alpha, p_mu, p_q, p_gamma, p_eta, p_epsilon)
    ! Intialise the parameters N, alpha, gamma, eta, epsilon, mu, q

    integer, intent(in) :: p_N
    real(dbl), intent(in) :: p_alpha, p_mu, p_q
    real(dbl), dimension(p_N, p_N), intent(in) :: p_gamma, p_eta, p_epsilon
    integer :: i, j ! Indices
    
    N = p_N
    alpha = p_alpha
    kappa = p_mu*p_q*p_q ! kappa = mu q^2

    ! Allocate arrays
    allocate( gamma(N, N) )
    allocate( eta(N, N) )
    allocate( epsilon(N, N) )
    allocate( delta(N, N) )

    ! Copy arrays
    gamma(:, :) = p_gamma(:, :)
    eta(:, :) = p_eta(:, :)
    epsilon(:, :) = p_epsilon(:, :)

    ! Populate the Kronecker delta
    outer: do i = 1, N
       inner: do j = 1, N
          if (i == j) then
             delta(i, j) = 1.0d0
          else
             delta(i, j) = 0.0d0
          end if
       end do inner
    end do outer

  end subroutine init_drude

  subroutine fin_drude()
    ! Deallocate the arrays
    if (allocated(gamma)) deallocate(gamma)
    if (allocated(epsilon)) deallocate(epsilon)
    if (allocated(eta)) deallocate(eta)
    if (allocated(delta)) deallocate(delta)
  end subroutine fin_drude
  
  subroutine calc_abc_vectors(rk, R, avecs, bvecs, cvecs)
    ! Calculates the vectors aij, bji, cij

    real(dbl), dimension(N, 3), intent(in) :: rk, R
    real(dbl), dimension(N, N, 3), intent(out) :: avecs, bvecs, cvecs

    integer :: i, j
    real(dbl), dimension(3) :: Rij
    ! Loop over all NxN-1 vectors
    outer: do i = 1, N-1
       Rij = R(i, :)
       inner: do j = i+1, N
          avecs(i, j, :) = Rij - rk(i, :)
          bvecs(j, i, :) = -Rij - rk(j, :)
          cvecs(i, j, :) = Rij - rk(i, :) + rk(j, :)
          Rij = Rij + R(j, :)
       end do inner
    end do outer

  end subroutine calc_abc_vectors

  subroutine calc_abc(rk, R, avals, bvals, cvals, avecs, bvecs, cvecs)
    ! Calculates the vectors a, b, c and their magnitudes

    real(dbl), dimension(N, 3), intent(in) :: rk, R
    real(dbl), dimension(N, N), intent(out) :: avals, bvals, cvals
    real(dbl), dimension(N, N, 3), intent(out) :: avecs, bvecs, cvecs
  
    integer :: i, j
    
    ! Get the vectors for a, b, c
    call calc_abc_vectors(rk, R, avecs, bvecs, cvecs)

    ! Compute magnitudes
    outer: do i = 1, N-1
       inner: do j = i+1, N
          avals(i, j) = sqrt(dot_product(avecs(i, j, :), avecs(i, j, :)))
          bvals(j, i) = sqrt(dot_product(bvecs(j, i, :), bvecs(j, i, :)))
          cvals(i, j) = sqrt(dot_product(cvecs(i, j, :), cvecs(i, j, :)))
       end do inner
    end do outer
          
  end subroutine calc_abc

  real(dbl) pure function tau(z, theta)
    ! Utility function for grad and laplacian, see notes

    real(dbl), intent(in) :: z, theta

    tau = theta*theta/(z*(theta+z)*(theta+z))
  end function tau

  real(dbl) pure function zeta(z, theta, tau)
    ! Utility function for laplacian, see notes
    ! Depends on value of tau for same z, theta

    real(dbl), intent(in) :: z, theta, tau

    zeta = (theta + 3.0d0*z)*tau/(theta + z)
  end function zeta

  real(dbl) pure function xi(tau, zeta)
    ! Utility for laplacian, see notes
    real(dbl), intent(in) :: tau, zeta

    xi = 3.0d0*tau - zeta
  end function xi

  real(dbl) pure function Jfunc(z, theta)
    ! Utility for wavefunction, see notes
    real(dbl), intent(in) :: z, theta

    Jfunc = theta*z/(theta + z)
  end function Jfunc
    
  subroutine calc_wf_derivatives(rk, R, wf, grads, laplacian)
    ! Computes the value of the wavefunction, wf,
    ! the gradients of the wavefunction, grads
    ! and the value of the laplacian, given the coordinates, rk and R

    real(dbl), dimension(N, 3), intent(in) :: rk, R
    real(dbl) :: wf, laplacian
    real(dbl), dimension(N, 3), intent(out) :: grads

    ! a, b, c vectors and magnitudes
    real(dbl), dimension(N, N) :: avals, bvals, cvals
    real(dbl), dimension(N, N, 3) :: avecs, bvecs, cvecs
    ! tau, zeta and xi values
    real(dbl), dimension(N, N) :: tau_a, tau_b, tau_c
    real(dbl), dimension(N, N) :: zeta_a, zeta_b, zeta_c
    real(dbl), dimension(N, N) :: xi_a, xi_b, xi_c

    integer :: i, j, k ! Loop indices
    real(dbl) :: Tijk
    
    ! Get the values of a, b, c
    call calc_abc(rk, R, avals, bvals, cvals, avecs, bvecs, cvecs)

    ! Calculate tau, zeta, and xi values
    outer: do i = 1, N-1
       inner: do j = i+1, N
          tau_a(i, j) = tau(avals(i, j), eta(i, j))
          tau_b(j, i) = tau(bvals(j, i), epsilon(j, i))
          tau_c(i, j) = tau(cvals(i, j), gamma(i, j))
          zeta_a(i, j) = zeta(avals(i, j), eta(i, j), tau_a(i, j))
          zeta_b(j, i) = zeta(bvals(j, i), epsilon(j, i), tau_b(j, i))
          zeta_c(i, j) = zeta(cvals(i, j), gamma(i, j), tau_c(i, j))
          xi_a(i, j) = xi(tau_a(i, j), zeta_a(i, j))
          xi_b(j, i) = xi(tau_b(j, i), zeta_b(j, i))
          xi_c(i, j) = xi(tau_c(i, j), zeta_c(i, j))
       end do inner
    end do outer
    
    ! Compute log wavefunction
    wf = 0.0d0
    outer2: do i = 1, N-1
       inner2: do j = i+1, N
          wf = wf + 0.5d0*Jfunc(cvals(i, j), gamma(i, j))
          wf = wf - Jfunc(avals(i, j), eta(i, j))
          wf = wf - Jfunc(bvals(j, i), epsilon(j, i))
       end do inner2
    end do outer2
    wf = kappa*wf
    do i = 1, N
       wf = wf - alpha*dot_product(rk(i, :), rk(i, :))
    end do
    ! Exponentiate
    wf = exp(wf)
    
    ! Compute the gradient prefactors
    kloop: do k = 1, N
       grads(k, :) = 0.0d0 ! Initialise
       
       outer3: do i = 1, N-1
          inner3: do j = i+1, N
             Tijk = delta(j, k) - delta(i, k)
             grads(k, :) = grads(k, :) + 0.5*Tijk*tau_c(i, j)*cvecs(i, j, :)
             grads(k, :) = grads(k, :) + delta(i, k)*tau_a(i, j)*avecs(i, j, :)
             grads(k, :) = grads(k, :) + delta(j, k)*tau_b(j, i)*bvecs(j, i, :)
          end do inner3
       end do outer3

       grads(k, :) = kappa*grads(k, :) - 2.0d0*alpha*rk(k, :)
    end do kloop

    ! Compute the Laplacian without factor of wf
    laplacian = 0.0d0 ! Initialise
    kloop2: do k = 1, N
       outer4: do i = 1, N-1
          inner4: do j = i+1, N
             Tijk = delta(i, k) + delta(j, k) ! Actually Tijk^2
             laplacian = laplacian + 0.5*Tijk*xi_c(i, j)
             laplacian = laplacian - delta(i, k)*xi_a(i, j)
             laplacian = laplacian - delta(j, k)*xi_b(j, i)
          end do inner4
       end do outer4
    end do kloop2

    laplacian = kappa * laplacian - 6.0d0*N*alpha ! Get factors of kappa and alpha

    ! Prefactor squared terms
    do k = 1, N
       laplacian = laplacian + dot_product(grads(k, :), grads(k, :))
    end do
    
    ! Add factors of wf
    grads = grads * wf
    laplacian = laplacian * wf
    
  end subroutine calc_wf_derivatives

end module n_drude_wf
  
    
    
    
  
    
