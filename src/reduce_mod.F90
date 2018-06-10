module reduce_mod

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
  public reduced_static_type
  public reduced_iap_type
  public reduced_state_type
  public full_reduced_state
  public half_reduced_state
  public full_reduced_static
  public half_reduced_static

  integer, allocatable :: full_reduce_factor(:)
  integer, allocatable :: half_reduce_factor(:)
  real, allocatable :: full_reduce_weight(:,:)
  real, allocatable :: half_reduce_weight(:,:)

  type reduced_static_type
    real, allocatable :: ghs(:,:,:)
  end type reduced_static_type

  type reduced_iap_type
    real, allocatable :: u(:,:,:)
    real, allocatable :: v(:,:,:)
    real, allocatable :: gd(:,:,:)
  end type reduced_iap_type

  type reduced_state_type
    real, allocatable :: u(:,:,:)
    real, allocatable :: v(:,:,:)
    real, allocatable :: gd(:,:,:)
    type(reduced_iap_type) iap
  end type reduced_state_type

  type(reduced_state_type), allocatable :: full_reduced_state(:)
  type(reduced_state_type), allocatable :: half_reduced_state(:)
  type(reduced_static_type), allocatable :: full_reduced_static(:)
  type(reduced_static_type), allocatable :: half_reduced_static(:)

contains

  subroutine reduce_init()

    integer i, j

    if (.not. allocated(full_reduce_factor)) allocate(full_reduce_factor(mesh%num_full_lat))
    if (.not. allocated(full_reduce_weight)) allocate(full_reduce_weight(maxval(zonal_reduce_factors),mesh%num_full_lat))
    if (.not. allocated(half_reduce_factor)) allocate(half_reduce_factor(mesh%num_half_lat))
    if (.not. allocated(half_reduce_weight)) allocate(half_reduce_weight(maxval(zonal_reduce_factors),mesh%num_half_lat))
    if (.not. allocated(full_reduced_state)) allocate(full_reduced_state(parallel%full_lat_start_idx_no_pole:parallel%full_lat_end_idx_no_pole))
    if (.not. allocated(half_reduced_state)) allocate(half_reduced_state(parallel%half_lat_start_idx:parallel%half_lat_end_idx))
    if (.not. allocated(full_reduced_static)) allocate(full_reduced_static(parallel%full_lat_start_idx_no_pole:parallel%full_lat_end_idx_no_pole))
    if (.not. allocated(half_reduced_static)) allocate(half_reduced_static(parallel%half_lat_start_idx:parallel%half_lat_end_idx))

    if (zonal_reduce_start_lat == -999) zonal_reduce_start_lat = 85.0

    full_reduce_factor(:) = 1
    half_reduce_factor(:) = 1
    if (use_zonal_reduce) then
      if (parallel%has_south_pole) then
        do j = 1, size(zonal_reduce_factors)
          if (zonal_reduce_factors(j) == 0) exit
          full_reduce_factor(parallel%full_lat_start_idx+j) = zonal_reduce_factors(j)
          half_reduce_factor(parallel%half_lat_start_idx+j-1) = zonal_reduce_factors(j)
        end do
      end if
      if (parallel%has_north_pole) then
        do j = 1, size(zonal_reduce_factors)
          if (zonal_reduce_factors(j) == 0) exit
          full_reduce_factor(parallel%full_lat_end_idx-j) = zonal_reduce_factors(j)
          half_reduce_factor(parallel%half_lat_end_idx-j+1) = zonal_reduce_factors(j)
        end do
      end if
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) /= 1) then
          full_reduce_weight(:,j) = 1.0 / full_reduce_factor(j)
        end if
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) /= 1) then
          half_reduce_weight(:,j) = 1.0 / half_reduce_factor(j)
        end if
      end do

      ! Allocate reduced data arrays.
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) /= 1) then
          call parallel_allocate(full_reduced_state(j)%u,      dim=[2,3], size=[full_reduce_factor(j),3], half_lon=.true.)
          call parallel_allocate(full_reduced_state(j)%v,      dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(full_reduced_state(j)%gd,     dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(full_reduced_state(j)%iap%u,  dim=[2,3], size=[full_reduce_factor(j),3], half_lon=.true.)
          call parallel_allocate(full_reduced_state(j)%iap%v,  dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(full_reduced_state(j)%iap%gd, dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(full_reduced_static(j)%ghs,   dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true.)
        end if
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) /= 1) then
          call parallel_allocate(half_reduced_state(j)%u,      dim=[2,3], size=[half_reduce_factor(j),3], half_lon=.true.)
          call parallel_allocate(half_reduced_state(j)%v,      dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(half_reduced_state(j)%gd,     dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(half_reduced_state(j)%iap%u,  dim=[2,3], size=[half_reduce_factor(j),3], half_lon=.true.)
          call parallel_allocate(half_reduced_state(j)%iap%v,  dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(half_reduced_state(j)%iap%gd, dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
          call parallel_allocate(half_reduced_static(j)%ghs,   dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
        end if
      end do
    end if

  end subroutine reduce_init

  subroutine reduce_final()

    if (allocated(full_reduce_factor)) deallocate(full_reduce_factor)
    if (allocated(full_reduce_weight)) deallocate(full_reduce_weight)
    if (allocated(half_reduce_factor)) deallocate(half_reduce_factor)
    if (allocated(half_reduce_weight)) deallocate(half_reduce_weight)
    if (allocated(full_reduced_state)) deallocate(full_reduced_state)
    if (allocated(half_reduced_state)) deallocate(half_reduced_state)
    if (allocated(full_reduced_static)) deallocate(full_reduced_static)
    if (allocated(half_reduced_static)) deallocate(half_reduced_static)

  end subroutine reduce_final

  subroutine reduce_run(state, static)

    type(state_type), intent(inout) :: state
    type(static_type), intent(inout) :: static

    integer j, k

    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) > 1) then
          do k = 1, full_reduce_factor(j)
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j-1),    full_reduced_state(j)%gd(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j  ),    full_reduced_state(j)%gd(:,k,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j+1),    full_reduced_state(j)%gd(:,k,3))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%u(:,j-1), full_reduced_state(j)%iap%u(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%u(:,j  ), full_reduced_state(j)%iap%u(:,k,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%u(:,j+1), full_reduced_state(j)%iap%u(:,k,3))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%v(:,j-1), full_reduced_state(j)%iap%v(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%v(:,j  ), full_reduced_state(j)%iap%v(:,k,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%iap%v(:,j+1), full_reduced_state(j)%iap%v(:,k,3))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, static%ghs(:,j-1),  full_reduced_static(j)%ghs(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, static%ghs(:,j  ),  full_reduced_static(j)%ghs(:,k,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, static%ghs(:,j+1),  full_reduced_static(j)%ghs(:,k,3))
            full_reduced_state(j)%iap%gd(:,k,:) = sqrt(full_reduced_state(j)%gd(:,k,:))
            full_reduced_state(j)%u(:,k,:) = full_reduced_state(j)%iap%u(:,k,:) / full_reduced_state(j)%iap%gd(:,k,:)
            full_reduced_state(j)%v(:,k,:) = full_reduced_state(j)%iap%v(:,k,:) / full_reduced_state(j)%iap%gd(:,k,:)
          end do
        end if
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if (half_reduce_factor(j) > 1) then
          do k = 1, half_reduce_factor(j)
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%gd(:,j-1),    half_reduced_state(j)%gd(:,k,1))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%gd(:,j  ),    half_reduced_state(j)%gd(:,k,2))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%gd(:,j+1),    half_reduced_state(j)%gd(:,k,3))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%u(:,j-1), half_reduced_state(j)%iap%u(:,k,1))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%u(:,j  ), half_reduced_state(j)%iap%u(:,k,2))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%u(:,j+1), half_reduced_state(j)%iap%u(:,k,3))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%v(:,j-1), half_reduced_state(j)%iap%v(:,k,1))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%v(:,j  ), half_reduced_state(j)%iap%v(:,k,2))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, state%iap%v(:,j+1), half_reduced_state(j)%iap%v(:,k,3))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, static%ghs(:,j-1),  half_reduced_static(j)%ghs(:,k,1))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, static%ghs(:,j  ),  half_reduced_static(j)%ghs(:,k,2))
            call average_raw_array_to_reduced_array_at_half_lat(j, k, static%ghs(:,j+1),  half_reduced_static(j)%ghs(:,k,3))
            half_reduced_state(j)%iap%gd(:,k,:) = sqrt(half_reduced_state(j)%gd(:,k,:))
            half_reduced_state(j)%u(:,k,:) = half_reduced_state(j)%iap%u(:,k,:) / half_reduced_state(j)%iap%gd(:,k,:)
            half_reduced_state(j)%v(:,k,:) = half_reduced_state(j)%iap%v(:,k,:) / half_reduced_state(j)%iap%gd(:,k,:)
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