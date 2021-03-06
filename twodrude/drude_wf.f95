module drude_wf
  !
  ! PURPOSE: To calculate the wavefunction, gradient, and Laplacian
  !          for a system of two drude oscillators.
  !

  implicit none
  
  ! Declare parameters
  integer, parameter :: dbl = selected_real_kind(15, 307)
  real(dbl) :: alpha, delta, gamma, epsilon, mu, q ! Parameters
  
contains

  subroutine init_drude(p1, p2, p3, p4, p5, p6)
    ! Intialise the parameters alpha, delta, gamma, epsilon, mu, q
    ! to the values p1, p2, p3, p4, p5, p6

    real(dbl), intent(in) :: p1, p2, p3, p4, p5, p6

    alpha = p1
    delta = p2
    gamma = p3
    epsilon = p4
    mu = p5
    q = p6

  end subroutine init_drude
  
  subroutine calc_abc_vectors(r1, r2, R, avec, bvec, cvec)
    ! Calculates the vectors a, b, c

    real(dbl), dimension(3), intent(in) :: r1, r2, R
    real(dbl), dimension(3), intent(out) :: avec, bvec, cvec

    avec = R - r1
    bvec = -R - r2
    cvec = R - r1 + r2

  end subroutine calc_abc_vectors

  subroutine calc_abc(r1, r2, R, a, b, c, avec, bvec, cvec)
    ! Calculates the vectors a, b, c and their magnitudes

    real(dbl), dimension(3), intent(in) :: r1, r2, R
    real(dbl), intent(out) :: a, b, c
    real(dbl), dimension(3), intent(out) :: avec, bvec, cvec

    ! Get the vectors for a, b, c
    call calc_abc_vectors(r1, r2, R, avec, bvec, cvec)

    ! Compute magnitudes
    a = sqrt(dot_product(avec, avec))
    b = sqrt(dot_product(bvec, bvec))
    c = sqrt(dot_product(cvec, cvec))

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

    zeta = (theta + 3d0*z)*tau/(theta + z)
  end function zeta
    
  subroutine calc_wf_derivatives(r1, r2, R, wf, grad1_wf, grad2_wf, laplacian)
    ! Computes the value of the wavefunction, wf,
    ! the gradients of the wavefunction, grad1_wf and grad2_wf,
    ! and the value of the laplacian, given the coordinates
    ! r1, r2, and R

    real(dbl), dimension(3), intent(in) :: r1, r2, R
    real(dbl) :: wf, laplacian
    real(dbl), dimension(3), intent(out) :: grad1_wf, grad2_wf

    ! a, b, c vectors and magnitudes
    real(dbl) :: a, b, c
    real(dbl), dimension(3) :: avec, bvec, cvec
    ! tau and zeta values
    real(dbl) :: tau_a, tau_b, tau_c, zeta_a, zeta_b, zeta_c
    real(dbl), dimension(3) :: tau_c2 ! Temporary storage of 0.5tau_c*cvec
    real(dbl) :: muq2 ! and mu*q^2

    muq2 = mu*q*q
    
    ! Get the values of a, b, c
    call calc_abc(r1, r2, R, a, b, c, avec, bvec, cvec)

    ! Calculate tau and zeta values
    tau_a = tau(a, delta)
    tau_b = tau(b, epsilon)
    tau_c = tau(c, gamma)
    zeta_a = zeta(a, delta, tau_a)
    zeta_b = zeta(b, epsilon, tau_b)
    zeta_c = zeta(c, gamma, tau_c)

    ! Compute log wavefunction
    wf = 0.5d0*gamma*c/(gamma + c)
    wf = wf  -delta*a/(delta + a)
    wf = wf  -epsilon*b/(epsilon + b)
    wf = muq2*wf - alpha*(dot_product(r1, r1) + dot_product(r2, r2))
    ! Exponentiate
    wf = exp(wf)
    
    ! Compute the gradient prefactors
    tau_c2 = 0.5d0*tau_c*cvec
    grad1_wf = -2d0*alpha*r1 +muq2*(-tau_c2 + tau_a*avec)
    grad2_wf = -2d0*alpha*r2 +muq2*( tau_c2 + tau_b*bvec)

    ! Compute the Laplacian without factor of wf
    laplacian = muq2*(3d0*(tau_c - tau_a - tau_b)  - zeta_c + zeta_a + zeta_b)
    laplacian = laplacian - 12.0d0*alpha + dot_product(grad1_wf, grad1_wf)
    laplacian = laplacian + dot_product(grad2_wf, grad2_wf)
    
    ! Add factors of wf
    grad1_wf = grad1_wf * wf
    grad2_wf = grad2_wf * wf
    laplacian = laplacian * wf
    
  end subroutine calc_wf_derivatives

end module drude_wf
  
    
    
    
  
    
