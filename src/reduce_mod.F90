module reduce_mod

  use data_mod
  use mesh_mod
  use log_mod
  use params_mod
  use parallel_mod
  use types_mod

  implicit none

  private

  public reduce_init
  public reduce_final
  public reduce_run
  public reduced_start_idx
  public reduced_end_idx
  public average_raw_array_to_reduced_array
  public append_reduced_tend_to_raw_tend
  public reduce_factor

  integer, allocatable :: reduce_factor(:)

contains

  subroutine reduce_init()

    real mean_dlon
    integer i, j, k

    if (.not. allocated(reduce_factor)) allocate(reduce_factor(mesh%num_full_lat))

    if (zonal_reduce_start_lat == -999) zonal_reduce_start_lat = 85.0

    reduce_factor(:) = 1
    mean_dlon = sum(coef%full_dlon) / mesh%num_full_lat / 2.0
    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (abs(mesh%full_lat_deg(j)) >= zonal_reduce_start_lat) then
          do k = 1, size(zonal_reduce_factors)
            if (mean_dlon / coef%full_dlon(j) < zonal_reduce_factors(k)) then
              reduce_factor(j) = zonal_reduce_factors(k)
              exit
            end if
            if (zonal_reduce_factors(k) == 0) then
              reduce_factor(j) = zonal_reduce_factors(k-1)
              exit
            end if
          end do
        end if
      end do
    end if

  end subroutine reduce_init

  subroutine reduce_final()

    if (allocated(reduce_factor)) deallocate(reduce_factor)

  end subroutine reduce_final

  subroutine reduce_run(state, static)

    type(state_type), intent(inout) :: state
    type(static_type), intent(inout) :: static

    integer j, k

    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        ! Calculate reduced state.
        if (reduce_factor(j) > 1) then
          do k = 1, reduce_factor(j)
            call average_raw_array_to_reduced_array(reduce_factor(j), k, state%gd(:,j), state%reduced_gd(:,k,j))
            call average_raw_array_to_reduced_array(reduce_factor(j), k, state%iap%u(:,j), state%iap%reduced_u(:,k,j))
            call average_raw_array_to_reduced_array(reduce_factor(j), k, state%iap%v(:,j-1), state%iap%reduced_v(:,k,j,1))
            call average_raw_array_to_reduced_array(reduce_factor(j), k, state%iap%v(:,j), state%iap%reduced_v(:,k,j,2))
            call average_raw_array_to_reduced_array(reduce_factor(j), k, static%ghs(:,j), static%reduced_ghs(:,k,j))
            state%iap%reduced_gd(:,k,j) = sqrt(state%reduced_gd(:,k,j))
            state%reduced_u(:,k,j) = state%iap%reduced_u(:,k,j) / state%iap%reduced_gd(:,k,j)
          end do
        end if
      end do
    end if

  end subroutine reduce_run

  integer function reduced_start_idx(j) result(start_idx)

    integer, intent(in) :: j

    start_idx = 1

  end function reduced_start_idx

  integer function reduced_end_idx(j) result(end_idx)

    integer, intent(in) :: j

    end_idx = mesh%num_full_lon / reduce_factor(j)

  end function reduced_end_idx

  subroutine average_raw_array_to_reduced_array(reduce_factor, k, raw_array, reduced_array)

    integer, intent(in) :: reduce_factor
    integer, intent(in) :: k
    real, intent(in) :: raw_array(:)
    real, intent(out) :: reduced_array(:)

    integer i, l, count, m, n

    n = parallel%lon_halo_width
    reduced_array(:) = 0.0
    l = n + 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_array) - parallel%lon_halo_width_for_reduce
      count = count + 1
      reduced_array(l) = reduced_array(l) + raw_array(i)
      if (count == reduce_factor) then
        reduced_array(l) = reduced_array(l) / reduce_factor
        l = l + 1
        count = 0
      end if
    end do

    ! Fill halo for reduced_array.
    m = (size(raw_array) - 2 * parallel%lon_halo_width_for_reduce) / reduce_factor + 2 * n
    reduced_array(1:n) = reduced_array(m-2*n+1:m-n)
    reduced_array(m-n+1:m) = reduced_array(1+n:2*n)

  end subroutine average_raw_array_to_reduced_array

  subroutine append_reduced_array_to_raw_array(reduce_factor, k, reduced_array, raw_array)

    integer, intent(in) :: reduce_factor
    integer, intent(in) :: k
    real, intent(in) :: reduced_array(:)
    real, intent(out) :: raw_array(:)

    integer i, l, n, count

    n = parallel%lon_halo_width
    l = n + 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_array) - parallel%lon_halo_width_for_reduce
      count = count + 1
      raw_array(i) = raw_array(i) + reduced_array(l)
      if (count == reduce_factor) then
        l = l + 1
        count = 0
      end if
    end do

  end subroutine append_reduced_array_to_raw_array

  subroutine append_reduced_tend_to_raw_tend(reduce_factor, k, reduced_tend, raw_tend)

    integer, intent(in) :: reduce_factor
    integer, intent(in) :: k
    real, intent(in) :: reduced_tend(:)
    real, intent(inout) :: raw_tend(:)

    integer i, l, count

    l = 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_tend) - parallel%lon_halo_width_for_reduce
      count = count + 1
      raw_tend(i) = raw_tend(i) + reduced_tend(l)
      if (count == reduce_factor) then
        l = l + 1
        count = 0
      end if
    end do

  end subroutine append_reduced_tend_to_raw_tend

end module reduce_mod