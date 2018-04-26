module params_mod

  implicit none

  ! Contant parameters
  real, parameter :: pi = atan(1.0) * 4.0
  real, parameter :: rad_to_deg = 180.0 / pi
  real, parameter :: deg_to_rad = pi / 180.0
  real, parameter :: omega = 2.0 * pi / 86400.0
  real, parameter :: radius = 6.37122e6
  real, parameter :: g = 9.80616

  integer num_lon
  integer num_lat
  integer :: subcycles = 4
  real(8) time_step_size
 
  integer :: days = 0
  integer :: hours = 0
  integer :: minutes = 0
  integer :: start_time(5) = [0, 0, 0, 0, 0]
  integer :: end_time(5) = [0, 0, 0, 0, 0]
  character(30) :: time_units = 'days'

  character(256) case_name
  character(256) case_desc
  character(256) author
  character(30) :: history_periods(1) = ['6 hours']
  character(30) :: restart_period = ''

  character(256) :: restart_file = ''

  ! Options:
  ! - predict-correct
  ! - runge-kutta
  character(30) time_scheme ! Time integration scheme
  integer time_order ! Time integration order (different schemes will have different meanings)
  logical qcon_modified ! Switch whether quadratic conservation modification is added

  ! Options:
  ! - csp1
  ! - csp2
  ! - isp
  ! - none
  character(30) split_scheme

  logical use_zonal_reduce
  integer zonal_reduce_factors(10)

  logical is_restart_run

  namelist /qconswm_params/ &
    num_lon, &
    num_lat, &
    subcycles, &
    days, &
    hours, &
    minutes, &
    start_time, &
    end_time, &
    time_units, &
    time_step_size, &
    case_name, &
    case_desc, &
    author, &
    history_periods, &
    restart_period, &
    restart_file, &
    time_scheme, &
    time_order, &
    qcon_modified, &
    split_scheme, &
    use_zonal_reduce, &
    zonal_reduce_factors

contains

  subroutine params_read(file_path)

    character(*), intent(in) :: file_path

    use_zonal_reduce = .true.
    zonal_reduce_factors(:) = 0

    open(10, file=file_path)
    read(10, nml=qconswm_params)
    close(10)

    is_restart_run = restart_file /= ''
    if (restart_period == '') restart_period = history_periods(1)

  end subroutine params_read

end module params_mod
