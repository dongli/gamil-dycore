module filter_mod

  use mesh_mod
  use parallel_mod
  use log_mod
  use params_mod

  implicit none

  private

  public filter_init
  public filter_run
  public filter_final

  public filter_full_zonal_tend
  public filter_half_zonal_tend

  real, allocatable :: wave_array(:)
  real, allocatable :: work_array(:)
  real, allocatable :: filter_factor(:)

  logical, allocatable :: filter_full_zonal_tend(:)
  logical, allocatable :: filter_half_zonal_tend(:)

contains

  subroutine filter_init()

    integer i, j, n, ierr

    if (.not. allocated(filter_full_zonal_tend)) allocate(filter_full_zonal_tend(mesh%num_full_lat))
    if (.not. allocated(filter_half_zonal_tend)) allocate(filter_half_zonal_tend(mesh%num_half_lat))

    filter_full_zonal_tend(:) = .false.
    filter_half_zonal_tend(:) = .false.
    if (parallel%has_south_pole) then
      do j = 1, size(zonal_tend_filter_tags)
        if (zonal_tend_filter_tags(j) == 1) then
          filter_full_zonal_tend(parallel%full_lat_start_idx+j) = .true.
          filter_half_zonal_tend(parallel%half_lat_start_idx+j-1) = .true.
        end if
      end do
    end if
    if (parallel%has_north_pole) then
      do j = 1, size(zonal_tend_filter_tags)
        if (zonal_tend_filter_tags(j) == 1) then
          filter_full_zonal_tend(parallel%full_lat_end_idx-j) = .true.
          filter_half_zonal_tend(parallel%full_lat_end_idx-j+1) = .true.
        end if
      end do
    end if
    ! do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
    !   print *, j, filter_full_zonal_tend(j)
    ! end do
    ! stop

    ! N + INT(LOG(REAL(N))) + 4
    if (.not. allocated(wave_array)) allocate(wave_array(mesh%num_full_lon + int(log(real(mesh%num_full_lon)) / log(2.0)) + 4))
    if (.not. allocated(work_array)) allocate(work_array(mesh%num_full_lon))
    if (.not. allocated(filter_factor)) allocate(filter_factor(mesh%num_full_lon))

    call rfft1i(mesh%num_full_lon, wave_array, size(wave_array), ierr)
    if (ierr /= 0) then
      call log_error('Failed to initialize FFTPACK!')
    end if

    filter_factor = 1.0
    n = mesh%num_full_lon / 2 - 1
    do i = zonal_tend_filter_cutoff_wavenumber + 2, n
      filter_factor(2 + 2 * i - 1) = 0.0
      filter_factor(2 + 2 * i) = 0.0
    end do

    call log_notice('Filter module is initialized.')

  end subroutine filter_init

  subroutine filter_run(x)

    real, intent(inout) :: x(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce)

    real local_x(mesh%num_full_lon)
    integer i, j, ierr

    j = 1
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
      local_x(j) = x(i)
      j = j + 1
    end do

    call rfft1f(mesh%num_full_lon, 1, local_x, mesh%num_full_lon, wave_array, size(wave_array), work_array, size(work_array), ierr)
    if (ierr /= 0) then
      call log_error('Failed to do forward FFT!')
    end if
    local_x(:) = local_x(:) * filter_factor(:)
    call rfft1b(mesh%num_full_lon, 1, local_x, mesh%num_full_lon, wave_array, size(wave_array), work_array, size(work_array), ierr)
        if (ierr /= 0) then
      call log_error('Failed to do backward FFT!')
    end if

    j = 1
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
      x(i) = local_x(j)
      j = j + 1
    end do

  end subroutine filter_run

  subroutine filter_final()

    if (allocated(filter_full_zonal_tend)) deallocate(filter_full_zonal_tend)
    if (allocated(filter_half_zonal_tend)) deallocate(filter_half_zonal_tend)
    if (allocated(wave_array)) deallocate(wave_array)
    if (allocated(work_array)) deallocate(work_array)

  end subroutine filter_final

end module filter_mod