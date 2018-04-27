module reduce_mod

  use data_mod
  use mesh_mod
  use log_mod
  use params_mod
  use parallel_mod
  use types_mod

  implicit none

  private

  public reduce_run
  public reduce_lb
  public reduce_ub
  public raw_array_to_reduced_array
  public reduced_tend_to_raw_tend

contains

  subroutine reduce_run(state, static, dt)

    type(state_type), intent(inout) :: state
    type(static_type), intent(inout) :: static
    real, intent(in) :: dt

    real cfl, mean_dlon
    integer i, j, k

    ! Calculate maximum CFL along each zonal circle.
    state%reduce_factor(:) = 1
    mean_dlon = sum(coef%full_dlon) / mesh%num_full_lat
    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        state%max_cfl(j) = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          cfl = dt * state%iap%gd(i,j) / coef%full_dlon(j)
          if (state%max_cfl(j) < cfl) state%max_cfl(j) = cfl
        end do
        ! Calculate reduced state.
        if (state%max_cfl(j) > 0.2) then
          ! Find reduce_factor based on excess of CFL.
          do k = 1, size(zonal_reduce_factors)
            !if (state%max_cfl(j) < zonal_reduce_factors(k)) then
            if (mean_dlon / coef%full_dlon(j) < zonal_reduce_factors(k)) then
              state%reduce_factor(j) = zonal_reduce_factors(k)
              exit
            end if
            if (zonal_reduce_factors(k) == 0) then
              ! call log_warning('Insufficient zonal_reduce_factors or time step size is too large!')
              state%reduce_factor(j) = zonal_reduce_factors(k - 1)
              exit
            end if
          end do
          call raw_array_to_reduced_array(state%gd(:,j), state%reduce_gd(:,j), state%reduce_factor(j))
          call raw_array_to_reduced_array(state%iap%u(:,j), state%iap%reduce_u(:,j), state%reduce_factor(j))
          call raw_array_to_reduced_array(state%iap%v(:,j-1), state%iap%reduce_v(:,j,1), state%reduce_factor(j))
          call raw_array_to_reduced_array(state%iap%v(:,j), state%iap%reduce_v(:,j,2), state%reduce_factor(j))
          call raw_array_to_reduced_array(static%ghs(:,j), static%reduce_ghs(:,j), state%reduce_factor(j))
          state%iap%reduce_gd(:,j) = sqrt(state%reduce_gd(:,j))
          state%reduce_u(:,j) = state%iap%reduce_u(:,j) / state%iap%reduce_gd(:,j)
          ! print *, j, state%max_cfl(j)
        end if
      end do
    end if
    call log_add_diag('num_reduce_zonal', count(state%reduce_factor /= 1))

  end subroutine reduce_run


  integer function reduce_lb(reduce_factor) result(lb)

    integer, intent(in) :: reduce_factor

    lb = 1 - parallel%lon_halo_width

  end function reduce_lb

  integer function reduce_ub(reduce_factor) result(ub)

    integer, intent(in) :: reduce_factor

    ub = mesh%num_full_lon / reduce_factor + parallel%lon_halo_width

  end function reduce_ub

  subroutine raw_array_to_reduced_array(raw_array, reduce_array, reduce_factor)

    ! NOTE: Here we allocate more-than-need space for reduce_array.
    real, intent(inout) :: raw_array(:)
    real, intent(out) :: reduce_array(:)
    integer, intent(in) :: reduce_factor

    integer i, j, count, m, n

    n = parallel%lon_halo_width
    reduce_array(:) = 0.0
    j = n + 1
    count = 0
    do i = 1 + n, size(raw_array) - n
      count = count + 1
      ! write(*, '(F20.5)', advance='no') raw_array(i)
      reduce_array(j) = reduce_array(j) + raw_array(i)
      if (count == reduce_factor) then
        reduce_array(j) = reduce_array(j) / reduce_factor
        ! write(*, '(" | ", F20.5)') reduce_array(j)
        j = j + 1
        count = 0
      end if
    end do

    ! Fill halo for reduce_array.
    m = (size(raw_array) - 2 * n) / reduce_factor + 2 * n
    reduce_array(1 : n) = reduce_array(m - 2 * n + 1 : m - n)
    reduce_array(m - n + 1 : m) = reduce_array(1 + n : 2 * n)

  end subroutine raw_array_to_reduced_array

  subroutine reduced_tend_to_raw_tend(reduce_tend, raw_tend, reduce_factor)

    real, intent(in) :: reduce_tend(:)
    real, intent(out) :: raw_tend(:)
    integer, intent(in) :: reduce_factor

    integer i, j, count, n

    n = parallel%lon_halo_width
    j = n + 1
    count = 0
    do i = 1 + n, size(raw_tend) - n
      count = count + 1
      raw_tend(i) = reduce_tend(j)
      if (count == reduce_factor) then
        j = j + 1
        count = 0
      end if
    end do

  end subroutine reduced_tend_to_raw_tend

end module reduce_mod