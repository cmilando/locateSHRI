subroutine get_f_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)

    implicit none
    
    ! Inputs / Outputs
    integer, intent(in)        :: np                          ! total number of site-point pairs
    integer, intent(in)        :: magic_n                          ! total number of site-point pairs
    integer, intent(in)        :: nsites                      ! total number of sites
    integer, intent(in)        :: S(nsites)                   ! S, a specific iteration of the magic indicator
    real(kind=8), intent(in)   :: site_pairs_sub(np, 6)       ! X, outcome data (columns: site_id, pop_id, dist, population, metric1, score)
    real(kind=8), intent(in)        :: row_lookup(nsites, 3)       ! row lookup (columns: site_id, start, end)
    real(kind=8), intent(in)   :: penalty
    
    real(kind=8)  :: SCORE                       ! Score to be calculated
    
    ! Local variables
    integer                   :: i, j, start_idx, end_idx, s_count
    integer, allocatable      :: pop_id_list(:), S_ones(:)
    real(kind=8), allocatable :: score_list(:)
    integer                   :: freq, idx
    real(kind=8)              :: score2

    SCORE = 0.0

    ! Allocate arrays to store intermediate data
    allocate(pop_id_list(np))
    allocate(score_list(np))
    allocate(S_ones(magic_n))

    s_count = 0

    ! Loop through S to find indices where S(i) == 1
    do i = 1, nsites
        if (S(i) == 1) then
            s_count = s_count + 1
            S_ones(s_count) = i
        end if
    end do

    ! Loop through selected sites (where S = 1) and gather corresponding site pairs
    idx = 0

    do i = 1, magic_n
        ! ok this is where its starting to get weird
        start_idx = int(row_lookup(S_ones(i), 2))
        end_idx   = int(row_lookup(S_ones(i), 3))
        !write(*, *) start_idx, end_idx

        do j = start_idx, end_idx
            idx = idx + 1
            pop_id_list(idx) = int(site_pairs_sub(j, 2))
            score_list(idx) = site_pairs_sub(j, 6)
        end do
    end do

    do i = 1, idx
        write(*, *) score_list(i)
    end do

    ! Calculate frequency of each population id
    do i = 1, idx
        freq = count(pop_id_list(1:idx) == pop_id_list(i))

        ! Adjust score based on frequency and penalty
        score2 = score_list(i) * 1.0 / (penalty * freq)
 
        ! fails here
        SCORE = SCORE + score2
    end do

    write(*, *) 'SCORE: ', SCORE

    ! Deallocate arrays
    deallocate(pop_id_list)
    deallocate(score_list)
    deallocate(S_ones)

end
