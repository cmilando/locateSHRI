subroutine select_non_overlapping_buffers3 (pts, n, selected_indices, buffer_dist)
    implicit none
  
    ! Define variables
    integer :: n ! Adjust the size as per the data frame size
    real(kind=8) :: pts(n, 2)
    real(kind=8) :: calc_dist
    integer :: selected_indices(n)
    integer :: i, j, selected_count
    logical :: within_radius
    real(kind=8) :: buffer_dist
  
    ! Initialize selected_indices to zero
    selected_indices = 0
    selected_count = 0

    !write(*, *) "    --", pts(1, 1)
    !write(*, *) "    --", pts(1, 2)

    do i = 1, n

        !write(*, *) "*****", i
        if (selected_count == 0) then
            within_radius = .false.
        else
            within_radius = .false.
            do j = 1, selected_count
                !write(*, *) "    --", pts(i, 1)
                !write(*, *) "    --", pts(i, 2)
                !write(*, *) "    --", pts(selected_indices(j), 1)
                !write(*, *) "    --", pts(selected_indices(j), 2)

                calc_dist = sqrt((pts(i, 1) - pts(selected_indices(j), 1))**2 + &
                             (pts(i, 2) - pts(selected_indices(j), 2))**2)
                !write(*, *) "   ", j, calc_dist
                if (calc_dist <= buffer_dist) then
                    within_radius = .true.
                    exit
                end if
            end do
        end if
        !write(*, *) "assess:", selected_count, within_radius
        if (selected_count == 0 .or. .not. within_radius) then
            !write(*, *) ">> UPDATED"
            selected_count = selected_count + 1
            selected_indices(selected_count) = i
        end if
    end do
  
    ! Output the selected indices
    !print *, "Selected indices: "
    !do i = 1, selected_count
    !  print *, selected_indices(i)
    !end do
  
  end 
  