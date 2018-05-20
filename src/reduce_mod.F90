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
  public reduced_start_idx_at_full_lat
  public reduced_end_idx_at_full_lat
  public reduced_start_idx_at_half_lat
  public reduced_end_idx_at_half_lat
  public append_reduced_tend_to_raw_tend_at_full_lat
  public append_reduced_tend_to_raw_tend_at_half_lat
  public full_reduce_factor
  public half_reduce_factor

  integer, allocatable :: full_reduce_factor(:)
  integer, allocatable :: half_reduce_factor(:)
  real, allocatable :: full_reduce_weight(:,:)
  real, allocatable :: half_reduce_weight(:,:)

contains

  subroutine reduce_init()

    real mean_dlon
    integer i, j, k

    if (.not. allocated(full_reduce_factor)) allocate(full_reduce_factor(mesh%num_full_lat))
    if (.not. allocated(full_reduce_weight)) allocate(full_reduce_weight(maxval(zonal_reduce_factors),mesh%num_full_lat))
    if (.not. allocated(half_reduce_factor)) allocate(half_reduce_factor(mesh%num_half_lat))
    if (.not. allocated(half_reduce_weight)) allocate(half_reduce_weight(maxval(zonal_reduce_factors),mesh%num_half_lat))

    if (zonal_reduce_start_lat == -999) zonal_reduce_start_lat = 85.0

    full_reduce_factor(:) = 1
    half_reduce_factor(:) = 1
    mean_dlon = sum(coef%full_dlon) / mesh%num_full_lat / 2.0
    if (use_zonal_reduce) then
      ! do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      !   if (abs(mesh%full_lat_deg(j)) >= zonal_reduce_start_lat) then
      !     do k = 1, size(zonal_reduce_factors)
      !       if (mean_dlon / coef%full_dlon(j) < zonal_reduce_factors(k)) then
      !         full_reduce_factor(j) = zonal_reduce_factors(k)
      !         exit
      !       end if
      !       if (zonal_reduce_factors(k) == 0) then
      !         full_reduce_factor(j) = zonal_reduce_factors(k-1)
      !         exit
      !       end if
      !     end do
      !   end if
      ! end do
      full_reduce_factor(2) = 6
      full_reduce_factor(3) = 5
      full_reduce_factor(4) = 4
      full_reduce_factor(5) = 3
      full_reduce_factor(6) = 2
      full_reduce_factor(177) = 2
      full_reduce_factor(177) = 3
      full_reduce_factor(178) = 4
      full_reduce_factor(179) = 5
      full_reduce_factor(180) = 6
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) /= 1) then
          ! do i = 1, full_reduce_factor(j)
          !   full_reduce_weight(i,j) = exp(-(i - 1 - (full_reduce_factor(j) - 1) * 0.5)**2 / (full_reduce_factor(j) * 0.5))
          ! end do
          ! full_reduce_weight(:full_reduce_factor(j),j) = full_reduce_weight(:full_reduce_factor(j),j) / sum(full_reduce_weight(:full_reduce_factor(j),j))
          full_reduce_weight(:,j) = 1.0 / full_reduce_factor(j)
        end if
      end do
      ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      !   if (abs(mesh%half_lat_deg(j)) >= zonal_reduce_start_lat) then
      !     do k = 1, size(zonal_reduce_factors)
      !       if (mean_dlon / coef%half_dlon(j) < zonal_reduce_factors(k)) then
      !         half_reduce_factor(j) = zonal_reduce_factors(k)
      !         exit
      !       end if
      !       if (zonal_reduce_factors(k) == 0) then
      !         half_reduce_factor(j) = zonal_reduce_factors(k-1)
      !         exit
      !       end if
      !     end do
      !   end if
      ! end do
      half_reduce_factor(1) = 6
      half_reduce_factor(2) = 5
      half_reduce_factor(3) = 4
      half_reduce_factor(4) = 3
      half_reduce_factor(5) = 2
      half_reduce_factor(176) = 2
      half_reduce_factor(177) = 3
      half_reduce_factor(178) = 4
      half_reduce_factor(179) = 5
      half_reduce_factor(180) = 6
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) /= 1) then
          ! do i = 1, half_reduce_factor(j)
          !   half_reduce_weight(i,j) = exp(-(i - 1 - (half_reduce_factor(j) - 1) * 0.5)**2 / (half_reduce_factor(j) * 0.5))
          ! end do
          ! half_reduce_weight(:half_reduce_factor(j),j) = half_reduce_weight(:half_reduce_factor(j),j) / sum(half_reduce_weight(:half_reduce_factor(j),j))
          half_reduce_weight(:,j) = 1.0 / half_reduce_factor(j)
        end if
      end do
    end if

  end subroutine reduce_init

  subroutine reduce_final()

    if (allocated(full_reduce_factor)) deallocate(full_reduce_factor)
    if (allocated(full_reduce_weight)) deallocate(full_reduce_weight)
    if (allocated(half_reduce_factor)) deallocate(half_reduce_factor)
    if (allocated(half_reduce_weight)) deallocate(half_reduce_weight)

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
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j), state%reduced_gd(:,k,j))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%u(:,j), state%iap%reduced_u(:,k,j))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%v(:,j-1), state%iap%reduced_v4u(:,k,j,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%v(:,j), state%iap%reduced_v4u(:,k,j,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, static%ghs(:,j), static%reduced_ghs(:,k,j))
            state%iap%reduced_gd(:,k,j) = sqrt(state%reduced_gd(:,k,j))
            state%reduced_u(:,k,j) = state%iap%reduced_u(:,k,j) / state%iap%reduced_gd(:,k,j)
          end do
        end if
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) > 1) then
          do k = 1, half_reduce_factor(j)
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%v(:,j), state%iap%reduced_v(:,k,j))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%u(:,j), state%reduced_u4v(:,k,j,1))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%u(:,j+1), state%reduced_u4v(:,k,j,2))
          end do
        end if
      end do
    end if

  end subroutine reduce_run

  integer function reduced_start_idx_at_full_lat(j) result(start_idx)

    integer, intent(in) :: j

    start_idx = 1

  end function reduced_start_idx_at_full_lat

  integer function reduced_end_idx_at_full_lat(j) result(end_idx)

    integer, intent(in) :: j

    end_idx = mesh%num_full_lon / full_reduce_factor(j)

  end function reduced_end_idx_at_full_lat

  integer function reduced_start_idx_at_half_lat(j) result(start_idx)

    integer, intent(in) :: j

    start_idx = 1

  end function reduced_start_idx_at_half_lat

  integer function reduced_end_idx_at_half_lat(j) result(end_idx)

    integer, intent(in) :: j

    end_idx = mesh%num_full_lon / half_reduce_factor(j)

  end function reduced_end_idx_at_half_lat

  subroutine average_raw_array_to_reduced_array_at_full_lat(j, k, raw_array, reduced_array)

    integer, intent(in) :: j
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
      reduced_array(l) = reduced_array(l) + raw_array(i) * full_reduce_weight(count,j)
      if (count == full_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if
    end do

    ! Fill halo for reduced_array.
    m = (size(raw_array) - 2 * parallel%lon_halo_width_for_reduce) / full_reduce_factor(j) + 2 * n
    reduced_array(1:n) = reduced_array(m-2*n+1:m-n)
    reduced_array(m-n+1:m) = reduced_array(1+n:2*n)

  end subroutine average_raw_array_to_reduced_array_at_full_lat

  subroutine average_raw_array_to_reduced_array_at_half_lat(j, k, raw_array, reduced_array)

    integer, intent(in) :: j
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
      reduced_array(l) = reduced_array(l) + raw_array(i) * half_reduce_weight(count,j)
      if (count == half_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if
    end do

    ! Fill halo for reduced_array.
    m = (size(raw_array) - 2 * parallel%lon_halo_width_for_reduce) / half_reduce_factor(j) + 2 * n
    reduced_array(1:n) = reduced_array(m-2*n+1:m-n)
    reduced_array(m-n+1:m) = reduced_array(1+n:2*n)

  end subroutine average_raw_array_to_reduced_array_at_half_lat

  subroutine append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, raw_tend)

    integer, intent(in) :: j
    integer, intent(in) :: k
    real, intent(in) :: reduced_tend(:)
    real, intent(inout) :: raw_tend(:)

    integer i, l, count

    l = 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_tend) - parallel%lon_halo_width_for_reduce
      count = count + 1
      raw_tend(i) = raw_tend(i) + reduced_tend(l) * full_reduce_weight(count,j)
      if (count == full_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if
    end do

  end subroutine append_reduced_tend_to_raw_tend_at_full_lat

  subroutine append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, raw_tend)

    integer, intent(in) :: j
    integer, intent(in) :: k
    real, intent(in) :: reduced_tend(:)
    real, intent(inout) :: raw_tend(:)

    integer i, l, count

    l = 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_tend) - parallel%lon_halo_width_for_reduce
      count = count + 1
      raw_tend(i) = raw_tend(i) + reduced_tend(l) * half_reduce_weight(count,j)
      if (count == half_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if
    end do

  end subroutine append_reduced_tend_to_raw_tend_at_half_lat

end module reduce_mod