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
  public reduce_test
  public reduced_full_start_idx
  public reduced_full_end_idx
  public reduced_half_start_idx
  public reduced_half_end_idx
  public average_raw_array_to_reduced_array
  public append_reduced_tend_to_raw_tend
  public full_reduce_factor
  public half_reduce_factor

  integer, allocatable :: full_reduce_factor(:)
  integer, allocatable :: half_reduce_factor(:)

contains

  subroutine reduce_init()

    real mean_dlon
    integer i, j, k

    if (.not. allocated(full_reduce_factor)) allocate(full_reduce_factor(mesh%num_full_lat))
    if (.not. allocated(half_reduce_factor)) allocate(half_reduce_factor(mesh%num_half_lat))

    if (zonal_reduce_start_lat == -999) zonal_reduce_start_lat = 85.0

    full_reduce_factor(:) = 1
    mean_dlon = sum(coef%full_dlon) / mesh%num_full_lat / 2.0
    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (abs(mesh%full_lat_deg(j)) >= zonal_reduce_start_lat) then
          do k = 1, size(zonal_reduce_factors)
            if (mean_dlon / coef%full_dlon(j) < zonal_reduce_factors(k)) then
              full_reduce_factor(j) = zonal_reduce_factors(k)
              exit
            end if
            if (zonal_reduce_factors(k) == 0) then
              full_reduce_factor(j) = zonal_reduce_factors(k-1)
              exit
            end if
          end do
        end if
      end do
    end if

    half_reduce_factor(:) = 1
    mean_dlon = sum(coef%half_dlon) / mesh%num_half_lat / 2.0
    if (use_zonal_reduce) then
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (abs(mesh%half_lat_deg(j)) >= zonal_reduce_start_lat) then
          do k = 1, size(zonal_reduce_factors)
            if (mean_dlon / coef%half_dlon(j) < zonal_reduce_factors(k)) then
              half_reduce_factor(j) = zonal_reduce_factors(k)
              exit
            end if
            if (zonal_reduce_factors(k) == 0) then
              half_reduce_factor(j) = zonal_reduce_factors(k-1)
              exit
            end if
          end do
        end if
      end do
    end if

  end subroutine reduce_init

  subroutine reduce_final()

    if (allocated(full_reduce_factor)) deallocate(full_reduce_factor)
    if (allocated(half_reduce_factor)) deallocate(half_reduce_factor)

  end subroutine reduce_final

  subroutine reduce_run(state, static)

    type(state_type), intent(inout) :: state
    type(static_type), intent(inout) :: static

    integer j, k

    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        ! Calculate reduced state.
        if (full_reduce_factor(j) > 1) then
          do k = 1, full_reduce_factor(j)
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, state%gd(:,j), state%reduced_gd(:,k,j))
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, state%iap%u(:,j), state%iap%reduced_u(:,k,j))
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, state%iap%v(:,j-1), state%iap%full_reduced_v(:,k,j,1))
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, state%iap%v(:,j), state%iap%full_reduced_v(:,k,j,2))
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, static%ghs(:,j), static%reduced_ghs(:,k,j))
            state%iap%reduced_gd(:,k,j) = sqrt(state%reduced_gd(:,k,j))
            state%reduced_u(:,k,j) = state%iap%reduced_u(:,k,j) / state%iap%reduced_gd(:,k,j)
          end do
        end if
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) > 1) then
          do k = 1, half_reduce_factor(j)
            call average_raw_array_to_reduced_array(half_reduce_factor(j), k, state%iap%v(:,j), state%iap%reduced_v(:,k,j))
            call average_raw_array_to_reduced_array(half_reduce_factor(j), k, state%u(:,j), state%half_reduced_u(:,k,j,1))
            call average_raw_array_to_reduced_array(half_reduce_factor(j), k, state%u(:,j+1), state%half_reduced_u(:,k,k,2))
          end do
        end if
      end do
    end if

  end subroutine reduce_run

  subroutine reduce_test(state)

    type(state_type), intent(inout) :: state

    integer i, j, k

    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) > 1) then
          do k = 1, full_reduce_factor(j)
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, state%gd(:,j), state%reduced_gd(:,k,j))
            call average_raw_array_to_reduced_array(full_reduce_factor(j), k, state%iap%u(:,j), state%iap%reduced_u(:,k,j))
          end do
          state%gd(:,j) = 0.0
          state%iap%u(:,j) = 0.0
          do k = 1, full_reduce_factor(j)
            call append_reduced_array_to_raw_array(full_reduce_factor(j), k, state%reduced_gd(:,k,j), state%gd(:,j))
            call append_reduced_array_to_raw_array(full_reduce_factor(j), k, state%iap%reduced_u(:,k,j), state%iap%u(:,j))
          end do
          call parallel_overlay_inner_halo(state%gd(:,j), left_halo=.true.)
          call parallel_overlay_inner_halo(state%iap%u(:,j), left_halo=.true.)
          state%gd(:,j) = state%gd(:,j) / full_reduce_factor(j)
          state%iap%u(:,j) = state%iap%u(:,j) / full_reduce_factor(j)
        end if
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) > 1) then
          do k = 1, half_reduce_factor(j)
            call average_raw_array_to_reduced_array(half_reduce_factor(j), k, state%iap%v(:,j), state%iap%reduced_v(:,k,j))
          end do
          state%iap%v(:,j) = 0.0
          do k = 1, half_reduce_factor(j)
            call append_reduced_array_to_raw_array(half_reduce_factor(j), k, state%iap%reduced_v(:,k,j), state%iap%v(:,j))
          end do
          call parallel_overlay_inner_halo(state%iap%v(:,j), left_halo=.true.)
          state%iap%v(:,j) = state%iap%v(:,j) / half_reduce_factor(j)
        end if
      end do
      call parallel_fill_halo(state%gd, all_halo=.true.)
      call parallel_fill_halo(state%iap%u, all_halo=.true.)
      call parallel_fill_halo(state%iap%v, all_halo=.true.)
      state%iap%gd(:,:) = sqrt(state%gd(:,:))
      state%u(:,:) = state%iap%u(:,:) / state%iap%gd(:,:)
      state%v(:,:) = state%iap%v(:,:) / state%iap%gd(:,:)
    end if

  end subroutine reduce_test

  integer function reduced_full_start_idx(j) result(start_idx)

    integer, intent(in) :: j

    start_idx = 1

  end function reduced_full_start_idx

  integer function reduced_full_end_idx(j) result(end_idx)

    integer, intent(in) :: j

    end_idx = mesh%num_full_lon / full_reduce_factor(j)

  end function reduced_full_end_idx

  integer function reduced_half_start_idx(j) result(start_idx)

    integer, intent(in) :: j

    start_idx = 1

  end function reduced_half_start_idx

  integer function reduced_half_end_idx(j) result(end_idx)

    integer, intent(in) :: j

    end_idx = mesh%num_half_lon / half_reduce_factor(j)

  end function reduced_half_end_idx

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