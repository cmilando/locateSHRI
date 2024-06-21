!--------------------------------------------------------------------------
! Subroutine for swapping discrete neighbors
!--------------------------------------------------------------------------
subroutine get_neighbor(x0, x1, l)
    
    implicit none

    integer :: l
    real    :: p0, p1
    integer :: p0_int, p1_int
    integer :: x0(l), x1(l)

    ! Initialize x1 to be x0
    x1 = x0

    ! now swap two indices
    ! random_number is [0, 1)
    ! FORTRAN arrays start at 1
    ! because this will give you [0,1] you need to adjust l by 1
    call random_number(p0)
    p0_int = floor(p0 * real(l)) + 1

 10 call random_number(p1)
    p1_int = floor(p1 * real(l)) + 1
    ! these cannot be the same so, try again if so
    ! Note: This is how you ensure you have all different indices
    if (x0(p1_int) .eq. x0(p0_int) ) go to 10

    x1(p0_int) = x0(p1_int)
    x1(p1_int) = x0(p0_int)

end