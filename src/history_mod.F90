module history_mod

  use io_mod
  use log_mod
  use mesh_mod
  use params_mod
  use parallel_mod
  use types_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write

  ! A-grid velocity
  real, allocatable :: u(:,:)
  real, allocatable :: v(:,:)

contains

  subroutine history_init()

    call io_create_dataset(desc=case_desc, file_prefix=case_name // '.h0')
    call io_add_meta('use_zonal_reduce', use_zonal_reduce)
    call io_add_meta('zonal_reduce_factors', pack(zonal_reduce_factors, zonal_reduce_factors /= 0))
    call io_add_meta('time_step_size', time_step_size)
    call io_add_meta('time_scheme', time_scheme)
    call io_add_meta('split_scheme', split_scheme)
    call io_add_meta('subcycles', subcycles)
    call io_add_dim('lon', size=mesh%num_full_lon)
    call io_add_dim('lat', size=mesh%num_full_lat)
    call io_add_dim('time')
    call io_add_var('u', long_name='u wind component', units='m s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('v', long_name='v wind component', units='m s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('gd', long_name='geopotential depth', units='m2 s-2', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('ghs', long_name='surface geopotential', units='m2 s-2', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('rf', long_name='reduce factor', units='1', dim_names=['lat ', 'time'])

    if (.not. allocated(u)) call parallel_allocate(u)
    if (.not. allocated(v)) call parallel_allocate(v)

    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write(state, static)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static

    integer i, j

    ! Convert wind from C grid to A grid.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        u(i,j) = 0.5 * (state%u(i,j) + state%u(i-1,j))
      end do
    end do
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        v(i,j) = 0.5 * (state%v(i,j) + state%v(i,j-1))
      end do
    end do
    call io_start_output()
    call io_output('lon', mesh%full_lon_deg)
    call io_output('lat', mesh%full_lat_deg)
    call io_output('u', u(:,:))
    call io_output('v', v(:,:))
    call io_output('gd', state%gd(:,:))
    call io_output('ghs', static%ghs(:,:))
    call io_output('rf', state%reduce_factor(:))
    call io_end_output()

  end subroutine history_write

end module history_mod
