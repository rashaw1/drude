subroutine ewald_real(ndim, N, rk, R, q, r_cut, alpha, box, energy)
  ! Calculates the real part of the Ewald summation
  implicit none
  
  integer,                             intent(in)  :: ndim, N
  double precision,  dimension(4, N),  intent(in)  :: rk, R
  double precision,  dimension(N),     intent(in)  :: q
  double precision,                    intent(in)  :: r_cut, alpha, box
  double precision,  dimension(N),     intent(out) :: energy

  double precision,  dimension(3)                  :: rk_i, r_i, rk_j, r_j, r_j_n
  double precision                                 :: r_ij
  integer                                          :: l_max, i, j, l_i, l_j, l_k

  ! Determine integer range, l_max
  l_max = ceiling( 2 * r_cut / box )

  ! Initialise energy value
  energy = 0d0

  ! Split by number of dimensions
  if ( ndim == 1 ) then

  else if ( ndim == 2) then

  else
     ! 3d

     ! Loop over charges
     do i = 1, N
        rk_i = rk(1:3, i)
        r_i = R(1:3, i)
        
        do j = 1, N
           rk_j = rk(1:3, j)
           r_j = R(1:3, j)
           
           ! Loop over real-space lattice vectors
           do l_i = -l_max, l_max
              do l_j = -l_max, l_max
                 do l_k = -l_max, l_max

                    ! Negative-negative
                    r_j_n(1) = rk_j(1) + l_i*box
                    r_j_n(2) = rk_j(2) + l_j*box
                    r_j_n(3) = rk_j(3) + l_k*box

                    r_j_n = rk_i - r_j_n
                    r_ij = sqrt( dot_product(r_j_n, r_j_n) )

                    if ( r_ij > 1e-10 .and. r_ij < r_cut) then 
                       energy(i) = energy(i) + q(i) * q(j) * erfc( alpha * r_ij) / r_ij
                    endif

                    ! Negative-positive
                    r_j_n(1) = r_j(1) + l_i*box
                    r_j_n(2) = r_j(2) + l_j*box
                    r_j_n(3) = r_j(3) + l_k*box

                    r_ij = sqrt( dot_product(rk_i - r_j_n, rk_i - r_j_n) )

                    if ( r_ij > 1e-10 .and. r_ij < r_cut) then
                       energy(i) = energy(i) - 2d0 * q(i) * q(j) * erfc( alpha * r_ij) / r_ij
                    endif

                    ! Positive-positive
                    r_j_n = r_i - r_j_n
                    r_ij = sqrt( dot_product(r_j_n, r_j_n) )

                    if ( r_ij > 1e-10 .and. r_ij < r_cut) then
                       energy(i) = energy(i) + q(i) * q(j) * erfc( alpha * r_ij) / r_ij
                    endif
                    
                    enddo
              enddo
           enddo

        enddo
     enddo
  endif
  
end subroutine ewald_real

subroutine ewald_reciprocal(ndim, N, rk, R, q, k_max, alpha, box, energy)
  ! Calculate the reciprocal space part of the Ewald summation
  implicit none
  
  integer,                              intent(in)  :: ndim, N, k_max
  double precision,   dimension(4, N),  intent(in)  :: rk, R
  double precision,   dimension(N),     intent(in)  :: q
  double precision,                     intent(in)  :: alpha, box
  double precision,                     intent(out) :: energy
  
  double precision                                  :: prefac, V, k_sq, ghat, h
  double precision                                  :: e_real, e_im
  double precision,   dimension(3)                  :: k
  integer                                           :: i, j, k_i, k_j, k_k
  double precision,   parameter                     :: PI = 3.14159265358979323846

  V = box**3
  prefac = 4d0 * PI / V
  energy = 0d0
  
  if (ndim == 1) then

  else if (ndim == 2) then
     
  else
     ! 3d

     do k_i = -k_max, k_max
        do k_j = -k_max, k_max
           do k_k = -k_max, k_max

              if (k_i == 0 .and. k_j == 0 .and. k_k == 0) cycle
              
              k(1) = 2d0*PI*k_i/box
              k(2) = 2d0*PI*k_j/box
              k(3) = 2d0*PI*k_k/box
              k_sq = dot_product(k, k)
              
              ghat = exp(-k_sq/(4d0*alpha**2))/k_sq

              e_real = 0d0
              e_im = 0d0
              do i = 1, N
                 ! Negative sites
                 h = dot_product(k, rk(1:3, i))
                 e_real = e_real - q(i) * cos(h)
                 e_im = e_im - q(i) * sin(h)

                 ! Positive sites
                 h = dot_product(k, R(1:3, i))
                 e_real = e_real + q(i) * cos(h)
                 e_im = e_im + q(i) * sin(h)
                 
!                 do j = 1, N
!                    ! Negative-negative
!                    h = dot_product(k, rk(1:3, i) - rk(1:3, j))
!                    e_real = e_real + q(i)*q(j) * cos(h)
!
!                    ! Negative-positive
!                    h = dot_product(k, rk(1:3, i) - R(1:3, j))
!                    e_real = e_real - 2d0*q(i)*q(j) * cos(h)
!
!                    ! Positive-positive
!                    h = dot_product(k, R(1:3, i) - R(1:3, j))
!                    e_real = e_real + q(i)*q(j) * cos(h)
!
!                 enddo
              enddo
              energy = energy + ghat * (e_real**2 + e_im**2)

           enddo
        enddo
     enddo         
     
     energy = energy * prefac
     
  endif

end subroutine ewald_reciprocal

subroutine ewald_sum(ndim, N, rk, R, q, r_cut, k_max, alpha, box, total_en)
  ! Calculates the Ewald summation, with dipole correction

  implicit none

  integer,                             intent(in)  :: ndim, N, k_max
  double precision,   dimension(4, N), intent(in)  :: rk, R
  double precision,   dimension(N),    intent(in)  :: q
  double precision,                    intent(in)  :: alpha, box, r_cut
  double precision,                    intent(out) :: total_en

  double precision,   dimension(N)                 :: real_en
  double precision,   dimension(3)                 :: dipole
  double precision                                 :: self_correct, dipole_en, recip_en
  integer                                          :: i
  double precision,   parameter                    :: PI = 3.14159265358979323846
  
  total_en = 0d0
  dipole = 0d0
  self_correct = 0d0

  call ewald_real(ndim, N, rk, R, q, r_cut, alpha, box, real_en)
  call ewald_reciprocal(ndim, N, rk, R, q, k_max, alpha, box, recip_en)

  write(*, *) real_en
  write(*, *) recip_en
  
  do i = 1, N
     total_en = total_en + 0.5d0 * real_en(i)
     self_correct = self_correct + 2d0*q(i)**2
     dipole = dipole + q(i) * ( R(1:3, i) - rk(1:3, i) ) 
  enddo

  self_correct = - ( alpha / sqrt( PI ) ) * self_correct
  dipole_en = ( 2d0 * PI / (3d0 * box**3) ) * dot_product(dipole, dipole)
  write(*, *) self_correct, dipole_en
  
  total_en = total_en + self_correct + 0.5d0 * recip_en + dipole_en

end subroutine ewald_sum

program ewald_test

  implicit none

  double precision, dimension(4, 4) :: rk, R
  double precision, dimension(4)    :: q
  double precision                  :: energy, r_cut, alpha, box
  integer                           :: k_max


  write(*, *) 'alpha, r_cut, k_max, box'
  read(*, *) alpha, r_cut, k_max, box
  
  rk = 0d0
  R = 0d0

  
! CsCl?
!  rk(1, 1) = 0.5d0*box
!  rk(2, 1) = 0.5d0*box
!  rk(3, 1) = 0.5d0*box
!  rk(4, 1) = sqrt(dot_product(rk(1:3, 1), rk(1:3, 1)))

! NaCl
!  rk(1, 1) = 0.5d0*box
!  rk(4, 1) = 0.5d0*box
!  rk(2, 2) = 0.5d0*box
!  rk(4, 2) = 0.5d0*box
!  rk(3, 3) = 0.5d0*box
!  rk(4, 3) = 0.5d0*box
!  rk(1, 4) = 0.5d0*box
!  rk(2, 4) = 0.5d0*box
!  rk(3, 4) = 0.5d0*box
!  rk(4, 4) = sqrt( dot_product(rk(1:3, 4), rk(1:3, 4)) )

! ZnS
  rk(1, 1) = 0.25d0*box
  rk(2, 1) = 0.25d0*box
  rk(3, 1) = 0.25d0*box
  rk(1, 2) = 0.25d0*box
  rk(2, 2) = 0.75d0*box
  rk(3, 2) = 0.75d0*box
  rk(1, 3) = 0.75d0*box
  rk(2, 3) = 0.25d0*box
  rk(3, 3) = 0.75d0*box
  rk(1, 4) = 0.75d0*box
  rk(2, 4) = 0.75d0*box
  rk(3, 4) = 0.25d0*box
  rk(4, 1) = sqrt( dot_product(rk(1:3, 1), rk(1:3, 1)) )
  rk(4, 2) = sqrt( dot_product(rk(1:3, 1), rk(1:3, 1)) )
  rk(4, 3) = sqrt( dot_product(rk(1:3, 1), rk(1:3, 1)) )
  rk(4, 4) = sqrt( dot_product(rk(1:3, 1), rk(1:3, 1)) )
  
  ! NaCl and ZnS
  R(1, 2) = 0.5d0*box
  R(2, 2) = 0.5d0*box
  R(1, 3) = 0.5d0*box
  R(3, 3) = 0.5d0*box
  R(2, 4) = 0.5d0*box
  R(3, 4) = 0.5d0*box
  R(4, 2) = sqrt( dot_product(R(1:3, 2), R(1:3, 2)) )
  R(4, 3) = R(4, 2)
  R(4, 4) = R(4, 2)

  q = 2d0
  
  call ewald_sum(3, 4, rk, R, q, r_cut, k_max, alpha, box, energy)

  write(*, *) energy, box*energy/8d0

end program ewald_test
        
                    
     
  
  
