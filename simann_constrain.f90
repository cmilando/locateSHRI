!--------------------------------------------------------------------------
! Author: CWM
! Date started: 6.21.2024
! Purpose: Subroutines for simulated annealing
! NOTE: This only works for S as an integer, and with no other confounders
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Subroutine for swapping discrete neighbors
!--------------------------------------------------------------------------
subroutine get_neighbor2(x0, x1, l, nS, cc, nzones, stopcond)
    
    implicit none

    integer :: l, nS, failure_count, failure_max, stopcond
    real    :: p0, p1
    integer :: p0_int, p1_int, idx0, i, count
    integer :: x0(l), x1(l), xtmp(l)
    integer :: cc(l)                   ! matrix for contraints
    integer :: nzones, zone_i                  ! 
    integer :: mult(l)

    ! integer
    failure_count = 0
    failure_max = 100

    ! Initialize x1 to be x0
    x1 = x0


    !write(*,*) "*****"

9   xtmp = x0

    ! now swap two indices
    ! random_number is [0, 1)
    ! FORTRAN arrays start at 1
    ! because this will give you [0,1] you need to adjust l by 1

    ! first get a random S=1
    call random_number(p0)
    p0_int = floor(p0 * real(nS)) + 1
    idx0 = 0
    count = 0
    !write(*,*) 'p0: ', p0_int
    
    ! Loop through x0 to find the p0_int-th 1
    do i = 1, l
        if (x0(i) == 1) then
            count = count + 1
            if (count == p0_int) then
                idx0 = i
                exit
            end if
        end if
    end do
    ! so now now idx0 = the xth 1
    !write(*,*) 'idx0: ', idx0

 10 call random_number(p1)
    p1_int = floor(p1 * real(l)) + 1
    ! these cannot be the same so, try again if so
    ! Note: This is how you ensure you have all different indices
    ! because this prevents swapping a 0 with a 0, or a 1 with a 1
    if (x0(p1_int) .eq. x0(idx0) ) go to 10
    !write(*,*) 'p1_int: ', p1_int

    ! here is where you decide what constraints to have
    ! so x0 is all 0s and 1s
    ! and cc is [1, 1, 1, ..., 2, 2, 2, ... ]
    ! and nzones is 3
    ! so x0 * cc elementwise will give [1, 0, 1, .., 2, 0, 2, ...]
    ! and I guess nzones could mean, there must be 1, 2, 3 in x0 * cc elementwise
    !do i = 1, l
    !    if(x0(i) .gt. 0) write(*,*) 'x0:', i, x0(i)
    !end do

    xtmp(idx0) = x0(p1_int)
    xtmp(p1_int) = x0(idx0)

    !do i = 1, l
    !    if(xtmp(i) .gt. 0) write(*,*) 'xtmp:', i, xtmp(i)
    !end do

    do i = 1,l
        mult(i) = 0
        if(xtmp(i) .gt. 0) mult(i) = cc(i)
    end do

    !do i = 1, l
    !    if(mult(i) .gt. 0) then
    !        write(*,*) 'mult:', i, mult(i)
    !    end if
    !end do

    failure_count = failure_count + 1
    if(failure_count >= failure_max) then
        stopcond = 1
        go to 88
    end if

    ! Check for the presence of all zone values in the result
    ! huh, if you get 2 out of sync this isn't going to work
    ! so this probably means you need to start with one point in each zone
    ! so maybe this works now that you are guaranteed to start with a new s=1 each time
    do zone_i = 1, nzones
        if (.not. any(mult == zone_i)) go to 9
    end do

    ! Return x1
    ! write(*,*) 'Found a good neighbor'
    x1 = xtmp
    
88  if(failure_count >= failure_max) write(*, *) 'fail exiting'
end

!--------------------------------------------------------------------------
! Cost function for SCORE
!--------------------------------------------------------------------------
subroutine get_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)

    implicit none
    
    ! Inputs / Outputs
    integer, intent(in)        :: np                          ! total number of site-point pairs
    integer, intent(in)        :: magic_n                     ! total number of site-point pairs
    integer, intent(in)        :: nsites                      ! total number of sites
    integer, intent(in)        :: S(nsites)                   ! S, a specific iteration of the magic indicator
    real(kind=8), intent(in)   :: site_pairs_sub(np, 6)       ! X, outcome data (columns: site_id, pop_id, dist, population, metric1, score)
    real(kind=8), intent(in)   :: row_lookup(nsites, 3)       ! row lookup (columns: site_id, start, end)
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

    ! Calculate frequency of each population id
    do i = 1, idx
        freq = count(pop_id_list(1:idx) == pop_id_list(i))

        ! Adjust score based on frequency and penalty
        score2 = score_list(i) * 1.0 / (penalty * freq)
        
        ! create a running sum
        SCORE = SCORE + score2
    end do

    SCORE = SCORE * (-1.0)

    !write(*, *) 'SCORE: ', SCORE

    ! Deallocate arrays
    deallocate(pop_id_list)
    deallocate(score_list)
    deallocate(S_ones)

end


! ------------------------------------------------------------------------------------
! Simulated Annealing algorithm
!
! informed by sci-kits sko SA.py
! ------------------------------------------------------------------------------------
subroutine simann_constrain(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE, cooling_rate, cc, nzones, verbose)
    
    !  Temp_max, Temp_min, cooling_rate

    implicit none
    
    ! Inputs / Outputs
    integer, intent(in)        :: np                     ! total number of site-point pairs
    integer, intent(in)        :: magic_n                ! total number of site-point pairs
    integer, intent(in)        :: nsites                 ! total number of sites

    real(kind=8), intent(in)   :: site_pairs_sub(np, 6)  ! X, outcome data (columns: site_id, pop_id, dist, population, metric1, score)
    real(kind=8), intent(in)   :: row_lookup(nsites, 3)  ! row lookup (columns: site_id, start, end)
    real(kind=8), intent(in)   :: penalty

    integer              :: S(nsites)              ! S, a specific iteration of the magic indicator
    integer              :: Sprev(nsites)                !  , a prevoious iteration of the magic indicator
    integer              :: Sbest(nsites)                !  , the best
    integer              :: cc(nsites)                   ! matrix for contraints
    integer              :: nzones
    integer              :: verbose
    integer :: stopcond

    real(kind=8)         :: SCORE                        !  
    real(kind=8)         :: SCOREprev                    !  
    real(kind=8)         :: SCOREbest                    !  
    real(kind=8)         :: SCOREbest1, SCOREbest2       ! generations of SCOREbest

    real(kind=8)         :: df             ! comparison of SCOREs
    real(kind=8)         :: rv             ! random variable to check

    real(kind=8) :: rel_tol, abs_tol

    ! Simulated annealing parmeters
    real(kind = 8) :: Temp_curr  ! 
    real(kind = 8) :: Temp_max   ! 
    real(kind = 8) :: Temp_min   ! 
    real(kind = 8) :: cooling_rate, curr_rel_tol, abs_score_diff
    integer :: LoC
    integer :: stay_counter
    integer :: max_stay_counter
    integer :: cycle_i
    integer :: iter_i  ! long of chain, may_stay_counter

    ! --------------------------------------------------------
    ! Initialize
    call random_seed()
    rel_tol          = 1.0e-2
    abs_tol          = 1.0e-3
    LoC              = 500
    max_stay_counter = 500
    stay_counter     = 0
    cycle_i          = 0
    SCOREbest        = 0.
    Temp_max         = 500000.
    Temp_min         = 0.0000000001
    !cooling_rate     = 0.95
    Temp_curr        = Temp_max

    ! **************
    ! INTIAL VALUES
    call get_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)
    Sprev   = S
    Sbest   = S
    SCOREprev  = SCORE
    SCOREbest  = SCORE
    SCOREbest1 = SCORE
    SCOREbest2 = SCORE
    ! **************

    stopcond = 0

    if(verbose .eq. 1) write(*, *) "Started"

    ! Simulated Annealing Loop
    ! converted to FORTRAN from python in scikit-opt on 10.4.2023 by CWM
    do 
        
        ! Update SCOREbest1 and SCOREbest2
        SCOREbest2 = SCOREbest1
        SCOREbest1 = SCORE

        ! -----------------
        ! ** PHASE 1 **
        ! -----------------
        ! Loop for all L in a specific chain
        do iter_i = 1, LoC

            !if(verbose .eq. 1) write(*, *) iter_i
            
            ! get a new x. so S becomes a new version of Sprev
            ! so x_current = Sprev
            ! and x_new = S
            call get_neighbor2(Sprev, S, nsites, magic_n, cc, nzones, stopcond)
            if(stopcond .eq. 1) go to 99
            !do i = 1, nsites
            !    if(S(i) .gt. 0) write(*,*) 'x1:', i, S(i)
            !end do

            ! calculate a new y. so this updates the value of SCORE
            ! similarly, y_current = SCOREprev
            ! and y_new = SCORE
            call get_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)

            ! ** METROPOLIS **
            df = SCORE - SCOREprev
            call random_number(rv)
            if (df < 0.0 .or. exp(-df / Temp_curr) > rv) then
                Sprev  = S             ! x_current = x_new
                SCOREprev = SCORE      ! y_current = y_new
                if (SCORE < SCOREbest) then
                    if(verbose .eq. 1) write(*, *) "best updated", cycle_i, SCORE
                    Sbest  = S
                    SCOREbest = SCORE
                end if
            end if
            
        end do

        ! -----------------
        ! ** PHASE 2 **
        ! -----------------
        ! Update cycle
        cycle_i = cycle_i + 1
        
        ! Cool down
        Temp_curr = Temp_curr * cooling_rate

        ! Get counter conditions
        ! so you have [SCOREbest2, SCOREbest1, and SCOREbest]
        ! SCOREbest was just updated, and so you should have the best past 3
        ! so now, check stay counter
        ! Need to check that this works correctly ....
        abs_score_diff = abs(SCOREbest1 - SCOREbest2)
        curr_rel_tol = rel_tol * max(abs(SCOREbest1), abs(SCOREbest2))
        if(abs_score_diff <= max(curr_rel_tol, abs_tol)) then
            stay_counter = stay_counter + 1
        else 
            stay_counter = 0
        end if

        ! check two stop conditions 
        ! Condition !: temperature   
         if (Temp_curr .lt. Temp_min) then
            if(verbose .eq. 1) write(*,*) "Cooled to final temperature"
            if(verbose .eq. 1) write(*,*) "Stay counter: ", stay_counter
            if(verbose .eq. 1) write(*,*) "abs score diff", abs_score_diff
            if(verbose .eq. 1) write(*,*) "curr_rel_tol", curr_rel_tol
            if(verbose .eq. 1) write(*,*) "abs_tol", abs_tol
            go to 99
        end if

        ! Condition 2: stability counter
        if (stay_counter .gt. max_stay_counter) then
            if(verbose .eq. 1) write(*,*) "Stayed unchanged in ", stay_counter, " iterations"
            go to 99
        end if

    end do

    ! Update the outputs
99  S = Sbest
    SCORE = SCOREbest
    if(verbose .eq. 1) write(*,*) "N. cycles = ", cycle_i

END