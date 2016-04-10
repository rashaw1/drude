module random
  use constants
  implicit none
contains
  subroutine random_init()
    integer :: t, n, i
    integer, allocatable, dimension(:) :: seed

    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(t)

    seed = t + 37 * [ (i, i = 1, n) ]
    call random_seed(put = seed)
    deallocate(seed)
  end subroutine random_init

  real(dbl) function random_uniform(lo_in, hi_in)
    real(dbl), intent(in), optional :: lo_in, hi_in
    real(dbl) :: lo, hi

    if (present(lo_in)) then
       lo = lo_in
    else
       lo = 0
    end if

    if (present(hi_in)) then
       hi = hi_in
    else
       hi = 1
    end if

    call random_number(random_uniform)

    random_uniform = random_uniform * (hi - lo) + lo
  end function random_uniform

  real(dbl) function random_normal(mean, stddev)
    ! Gives a random number sampled from a normal
    ! distribution N(mean, stddev)
    real(dbl), intent(in) :: mean, stddev
    real(dbl) :: r, theta, sigma
    real(dbl), dimension(2) :: temp

    if (stddev == 0d0) then
       sigma = 1d0
    else if (stddev < 0d0) then
       sigma = -1*stddev
    else
       sigma = stddev
    end if

    call random_number(temp)
    r = sqrt(-2.0d0*log(temp(1)))
    theta = 2.0d0*PI*temp(2)
    random_normal = mean+sigma*r*sin(theta)
  end function random_normal

end module random
