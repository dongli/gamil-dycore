module restart_mod

  use datetime_mod
  use time_mod
  use io_mod
  use log_mod
  use mesh_mod
  use params_mod
  use parallel_mod
  use types_mod

  implicit none

  private

  public restart_init
  public restart_write
  public restart_read

contains

  subroutine restart_init()

    call io_create_dataset('restart', desc='Restart file', file_prefix=trim(case_name) // '.r', period=restart_period)
    call io_add_dim('lon', size=mesh%num_full_lon, dataset_name='restart')
    call io_add_dim('ilon', size=mesh%num_half_lon, dataset_name='restart')
    call io_add_dim('lat', size=mesh%num_full_lat, dataset_name='restart')
    call io_add_dim('ilat', size=mesh%num_half_lat, dataset_name='restart')
    call io_add_dim('time', dataset_name='restart')
    call io_add_var('u', long_name='u wind component', units='m s-1', dim_names=['ilon', 'lat ', 'time'], dataset_name='restart')
    call io_add_var('v', long_name='v wind component', units='m s-1', dim_names=['lon ', 'ilat', 'time'], dataset_name='restart')
    call io_add_var('gd', long_name='geopotential depth', units='m2 s-2', dim_names=['lon ', 'lat ', 'time'], dataset_name='restart')
    call io_add_var('ghs', long_name='surface geopotential', units='m2 s-2', dim_names=['lon ', 'lat ', 'time'], dataset_name='restart')

  end subroutine restart_init

  subroutine restart_read(state, static)

    type(state_type), intent(inout) :: state
    type(static_type), intent(inout) :: static

    ! call io_create_dataset('restart', file_path=restart_file, mode='input')
    ! call io_start_input('restart')
    ! call time_reset_start_time(datetime(io_get_meta('restart_time', 'restart')))
    ! call log_notice('Reset time to ' // trim(curr_time_format) // '.')
    ! call io_input('u', state%u, 'restart')
    ! call parallel_fill_halo(state%u(:,:), all_halo=.true.)
    ! call io_input('v', state%v, 'restart')
    ! call parallel_fill_halo(state%v(:,:), all_halo=.true.)
    ! call io_input('gd', state%gd, 'restart')
    ! call parallel_fill_halo(state%gd(:,:), all_halo=.true.)
    ! call io_input('ghs', static%ghs, 'restart')
    ! call parallel_fill_halo(static%ghs, all_halo=.true.)

  end subroutine restart_read

  subroutine restart_write(state, static)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static

    ! call io_add_meta('restart_time', curr_time_format, 'restart')
    ! call io_add_meta('elapsed_seconds', time_elapsed_seconds(), 'restart')
    ! call io_start_output('restart')
    ! call io_output('lon', mesh%full_lon_deg(:), 'restart')
    ! call io_output('ilon', mesh%half_lon_deg(:), 'restart')
    ! call io_output('lat', mesh%full_lat_deg(:), 'restart')
    ! call io_output('ilat', mesh%half_lat_deg(:), 'restart')
    ! call io_output('u', state%u(:,:), 'restart')
    ! call io_output('v', state%v(:,:), 'restart')
    ! call io_output('gd', state%gd(:,:), 'restart')
    ! call io_output('ghs', static%ghs(:,:), 'restart')
    ! call io_end_output('restart')

  end subroutine restart_write

end module restart_mod
