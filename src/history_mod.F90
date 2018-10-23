module history_mod

  use io_mod
  use log_mod
  use mesh_mod
  use params_mod
  use parallel_mod
  use types_mod
  use diag_mod
  use string_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write

  ! A-grid velocity
  real, allocatable :: u(:,:)
  real, allocatable :: v(:,:)

  interface history_write
    module procedure history_write_state
    module procedure history_write_tendency
  end interface history_write

contains

  subroutine history_init()

    call io_create_dataset(desc=case_desc, file_prefix=case_name // '.h0')
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
    call io_add_var('vor', long_name='relative vorticity', units='s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('div', long_name='divergence', units='s-1', dim_names=['lon ', 'lat ', 'time'])

    call io_create_dataset(name='debug', desc=case_desc, file_prefix=case_name // '.debug')
    call io_add_dim('lon', 'debug', size=mesh%num_full_lon)
    call io_add_dim('lat', 'debug', size=mesh%num_full_lat)
    call io_add_dim('ilon', 'debug', size=mesh%num_half_lon)
    call io_add_dim('ilat', 'debug', size=mesh%num_half_lat)
    call io_add_dim('time', 'debug')
    call io_add_var('u_adv_lon', 'debug', long_name='u_adv_lon', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('u_adv_lat', 'debug', long_name='u_adv_lat', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('v_adv_lon', 'debug', long_name='v_adv_lon', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('v_adv_lat', 'debug', long_name='v_adv_lat', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('fv', 'debug', long_name='fv', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('fu', 'debug', long_name='fu', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('u_pgf', 'debug', long_name='u_pgf', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('v_pgf', 'debug', long_name='v_pgf', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('mass_div_lon', 'debug', long_name='mass_div_lon', units='', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('mass_div_lat', 'debug', long_name='mass_div_lat', units='', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('du', 'debug', long_name='du', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('dv', 'debug', long_name='dv', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('dgd', 'debug', long_name='dgd', units='', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('iap_u', 'debug', long_name='IAP u', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('iap_v', 'debug', long_name='IAP v', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('iap_gd', 'debug', long_name='geopotential depth', units='', dim_names=['lon ', 'lat ', 'time'])

    if (.not. allocated(u)) call parallel_allocate(u)
    if (.not. allocated(v)) call parallel_allocate(v)

    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write_state(state, static, diag)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static
    type(diag_type), intent(in) :: diag

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
    call io_output('lon', mesh%full_lon_deg(:))
    call io_output('lat', mesh%full_lat_deg(:))
    call io_output('u', u(:,:))
    call io_output('v', v(:,:))
    call io_output('gd', state%gd(:,:))
    call io_output('ghs', static%ghs(:,:))
    call io_output('vor', diag%vor(:,:))
    call io_output('div', diag%div(:,:))
    call io_end_output()

  end subroutine history_write_state

  subroutine history_write_tendency(state, tend, tag)

    type(state_type), intent(in) :: state
    type(tend_type), intent(in) :: tend
    integer, intent(in) :: tag

    call io_start_output('debug', trim(to_string(tag)))
    call io_output('lon', mesh%full_lon_deg(:), 'debug')
    call io_output('lat', mesh%full_lat_deg(:), 'debug')
    call io_output('ilon', mesh%half_lon_deg(:), 'debug')
    call io_output('ilat', mesh%half_lat_deg(:), 'debug')
    call io_output('u_adv_lon', tend%u_adv_lon(:,:), 'debug')
    call io_output('u_adv_lat', tend%u_adv_lat(:,:), 'debug')
    call io_output('v_adv_lon', tend%v_adv_lon(:,:), 'debug')
    call io_output('v_adv_lat', tend%v_adv_lat(:,:), 'debug')
    call io_output('fv', tend%fv(:,:), 'debug')
    call io_output('fu', tend%fu(:,:), 'debug')
    call io_output('u_pgf', tend%u_pgf(:,:), 'debug')
    call io_output('v_pgf', tend%v_pgf(:,:), 'debug')
    call io_output('mass_div_lon', tend%mass_div_lon(:,:), 'debug')
    call io_output('mass_div_lat', tend%mass_div_lat(:,:), 'debug')
    call io_output('du', tend%du(:,:), 'debug')
    call io_output('dv', tend%dv(:,:), 'debug')
    call io_output('dgd', tend%dgd(:,:), 'debug')
    call io_output('iap_u', state%iap%u(:,:), 'debug')
    call io_output('iap_v', state%iap%v(:,:), 'debug')
    call io_output('iap_gd', state%iap%gd(:,:), 'debug')
    call io_end_output('debug')

  end subroutine history_write_tendency

end module history_mod
