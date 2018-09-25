module dycore_mod

  use params_mod, time_scheme_in => time_scheme, split_scheme_in => split_scheme
  use data_mod
  use reduce_mod
  use log_mod
  use types_mod
  use mesh_mod
  use time_mod
  use parallel_mod
  use io_mod
  use diag_mod
  use history_mod
  use restart_mod
  use diffusion_mod
  use filter_mod
  use pole_a_grid_mod
  use space_operators_mod

  implicit none

  private

  public dycore_init
  public dycore_restart
  public dycore_run
  public dycore_final

  ! 1: predict_correct
  ! 2: runge_kutta
  ! 3: leap_frog
  ! 4: middle_point
  integer time_scheme
  ! 1: csp1: first order conservative split
  ! 2: csp2: second order conservative split
  ! 3: isp: inproved second order conservative split
  integer split_scheme

  integer, parameter :: half_time_idx = 0

  interface
    subroutine integrator_interface(time_step_size, old_time_idx, new_time_idx, pass)
      real, intent(in) :: time_step_size
      integer, intent(in), optional :: old_time_idx
      integer, intent(in), optional :: new_time_idx
      integer, intent(in), optional :: pass
    end subroutine
  end interface

  procedure(integrator_interface), pointer :: integrator

  ! For debug
  integer :: tag = 0

contains

  subroutine dycore_init()

    if (case_name == '') then
      call log_error('case_name is not set!')
    end if

    call log_init()
    call mesh_init()
    call time_init()
    call parallel_init()
    call io_init()
    call diag_init()
    call history_init()
    call restart_init()
    call data_init()
    call reduce_init()
    call diffusion_init()
    call filter_init()
    call pole_a_grid_init()

    select case (time_scheme_in)
    case ('predict_correct')
      time_scheme = 1
    case ('runge_kutta')
      time_scheme = 2
    case ('leap_frog')
      time_scheme = 3
    case ('middle_point')
      time_scheme = 4
    case default
      call log_error('Unknown time_scheme ' // trim(time_scheme_in) // '!')
    end select

    select case (split_scheme_in)
    case ('csp1')
      split_scheme = 1
    case ('csp2')
      split_scheme = 2
    case ('isp')
      split_scheme = 3
    case default
      split_scheme = 0
      call log_notice('No fast-slow split.')
    end select

    select case (time_scheme)
    case (1)
      integrator => predict_correct
    case (2)
      integrator => runge_kutta
    case (3)
      integrator => leap_frog
    case (4)
      integrator => middle_point
    case default
      call log_error('Unknown time scheme!')
    end select

    call log_notice('Dycore module is initialized.')

  end subroutine dycore_init

  subroutine dycore_restart()

    call restart_read(state(old_time_idx), static)

  end subroutine dycore_restart

  subroutine dycore_run()

    call reset_cos_lat_at_poles()

    call pole_a_grid_init_state(state(old_time_idx))
    call iap_transform(state(old_time_idx))

    call diag_run(state(old_time_idx))
    call output(state(old_time_idx))
    call log_add_diag('total_mass', diag%total_mass)
    call log_add_diag('total_energy', diag%total_energy)
    call log_step()

    do while (.not. time_is_finished())
      tag = 0
      call time_integrate()
      call time_advance()
      call diag_run(state(old_time_idx))
      call output(state(old_time_idx))
      call log_add_diag('total_mass', diag%total_mass)
      call log_add_diag('total_energy', diag%total_energy)
      call log_step()
    end do

  end subroutine dycore_run

  subroutine dycore_final()

    call mesh_final()
    call parallel_final()
    call diag_final()
    call history_final()
    call data_final()
    call reduce_final()
    call diffusion_final()
    call filter_final()

    call log_notice('Dycore module is finalized.')

  end subroutine dycore_final

  subroutine reset_cos_lat_at_poles()

    integer j

    j = parallel%full_lat_south_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_south_pole_idx) * 0.25
    coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)

    j = parallel%full_lat_north_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_north_pole_idx) * 0.25
    coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)

  end subroutine reset_cos_lat_at_poles

  subroutine output(state)

    type(state_type), intent(inout) :: state

    if (time_is_alerted('hist0.output')) call history_write(state, static, diag)
    if (time_is_alerted('restart.output')) call restart_write(state, static)

  end subroutine output

  subroutine update_state(dt, tend, old_state, new_state)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
      end do
    end do

    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)

    do j = parallel%full_lat_lb, parallel%full_lat_ub
      do i = parallel%full_lon_lb_for_reduce, parallel%full_lon_ub_for_reduce
        new_state%iap%gd(i,j) = sqrt(new_state%gd(i,j))
      end do
    end do

    ! Update IAP wind state on A grids at South Pole.
    do j = full_lat_start_idx_a_grid(south_pole), full_lat_end_idx_a_grid(south_pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%u_a(i,j) = old_state%iap%u_a(i,j) + dt * tend%du_a(i,j)
        new_state%iap%v_a(i,j) = old_state%iap%v_a(i,j) + dt * tend%dv_a(i,j)
      end do
    end do

    ! Update IAP wind state on A grids at North Pole.
    do j = full_lat_start_idx_a_grid(north_pole), full_lat_end_idx_a_grid(north_pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%u_a(i,j) = old_state%iap%u_a(i,j) + dt * tend%du_a(i,j)
        new_state%iap%v_a(i,j) = old_state%iap%v_a(i,j) + dt * tend%dv_a(i,j)
      end do
    end do

    ! Update IAP wind state on C grids.
    do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%iap%u_c(i,j) = old_state%iap%u_c(i,j) + dt * tend%du_c(i,j)
      end do
    end do
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%v_c(i,j) = old_state%iap%v_c(i,j) + dt * tend%dv_c(i,j)
      end do
    end do

    ! Transform from IAP to normal state on A grids at South Pole.
    do j = full_lat_start_idx_a_grid(south_pole), full_lat_end_idx_a_grid(south_pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%u_a(i,j) = new_state%iap%u_a(i,j) / new_state%iap%gd(i,j)
        new_state%v_a(i,j) = new_state%iap%v_a(i,j) / new_state%iap%gd(i,j)
      end do
    end do

    ! Transform from IAP to normal state on A grids at North Pole.
    do j = full_lat_start_idx_a_grid(north_pole), full_lat_end_idx_a_grid(north_pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%u_a(i,j) = new_state%iap%u_a(i,j) / new_state%iap%gd(i,j)
        new_state%v_a(i,j) = new_state%iap%v_a(i,j) / new_state%iap%gd(i,j)
      end do
    end do

    ! Transform from IAP to normal state on C grids.
    do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u_c(i,j) = new_state%iap%u_c(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i+1,j))
      end do
    end do
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v_c(i,j) = new_state%iap%v_c(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i,j+1))
      end do
    end do

    call parallel_fill_halo(new_state%iap%u_a(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v_a(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u_a(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v_a(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%u_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v_c(:,:), all_halo=.true.)

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0
    do j = full_lat_start_idx_a_grid(south_pole), full_lat_end_idx_a_grid(south_pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%du_a(i,j) * tend2%du_a(i,j) * mesh%full_cos_lat(j)
        res = res + tend1%dv_a(i,j) * tend2%dv_a(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = full_lat_start_idx_a_grid(north_pole), full_lat_end_idx_a_grid(north_pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%du_a(i,j) * tend2%du_a(i,j) * mesh%full_cos_lat(j)
        res = res + tend1%dv_a(i,j) * tend2%dv_a(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend1%du_c(i,j) * tend2%du_c(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dv_c(i,j) * tend2%dv_c(i,j) * mesh%half_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dgd(i,j) * tend2%dgd(i,j) * mesh%full_cos_lat(j)
      end do
    end do

  end function inner_product

  subroutine time_integrate()

    real subcycle_time_step_size
    integer subcycle, time_idx1, time_idx2

    subcycle_time_step_size = time_step_size / subcycles
    time_idx1 = 0
    time_idx2 = old_time_idx

    select case (split_scheme)
    case (2) ! csp2
      call integrator(0.5 * time_step_size, old_time_idx, time_idx1, slow_pass)
      do subcycle = 1, subcycles
        call integrator(subcycle_time_step_size, time_idx1, time_idx2, fast_pass)
        call time_swap_indices(time_idx1, time_idx2)
      end do
      call integrator(0.5 * time_step_size, time_idx1, new_time_idx, slow_pass)
    case default
      call integrator(time_step_size)
    end select

    if (use_diffusion) then
      call diffusion_run(time_step_size, state(new_time_idx))
    end if

  end subroutine time_integrate

  subroutine middle_point(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

    integer old, new, half, pass, iteration
    real dt, e1, e2

    if (present(old_time_idx_)) then
      old = old_time_idx_
    else
      old = old_time_idx
    end if
    if (present(new_time_idx_)) then
      new = new_time_idx_
    else
      new = new_time_idx
    end if
    if (present(pass_)) then
      pass = pass_
    else
      pass = all_pass
    end if
    half = half_time_idx
    dt = time_step_size

    call copy_state(state(old), state(new))

    e1 = diag_total_energy(state(old))
    do iteration = 1, 8
      call average_state(state(old), state(new), state(half))
      call space_operators_run(state(half), tend(old), dt, pass)
      call update_state(dt, tend(old), state(old), state(new))
      e2 = diag_total_energy(state(new))
      if (abs(e2 - e1) * 2 / (e2 + e1) < 5.0e-15) then
        exit
      end if
    end do

  end subroutine middle_point

  subroutine predict_correct(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

    integer old, new, pass
    real dt, ip1, ip2, beta

    if (present(old_time_idx_)) then
      old = old_time_idx_
    else
      old = old_time_idx
    end if
    if (present(new_time_idx_)) then
      new = new_time_idx_
    else
      new = new_time_idx
    end if
    if (present(pass_)) then
      pass = pass_
    else
      pass = all_pass
    end if
    dt = time_step_size * 0.5

    ! Do first predict step.
    call space_operators_run(state(old), tend(old), dt, pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do second predict step.
    call space_operators_run(state(new), tend(old), dt, pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do correct step.
    call space_operators_run(state(new), tend(new), dt, pass)

    ip1 = inner_product(tend(old), tend(new))
    ip2 = inner_product(tend(new), tend(new))
    beta = merge(ip1 / ip2, 1.0, qcon_modified .and. ip1 /= 0.0 .and. ip2 /= 0.0)
    call log_add_diag('beta', beta)

    dt = time_step_size * beta
    call update_state(dt, tend(new), state(old), state(new))

  end subroutine predict_correct

  subroutine runge_kutta(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

    integer old, new, tmp, pass, i, j
    real dt, ip00, ip12, ip23, ip34, beta

    tmp = -1
    if (present(old_time_idx_)) then
      old = old_time_idx_
    else
      old = old_time_idx
    end if
    if (present(new_time_idx_)) then
      new = new_time_idx_
    else
      new = new_time_idx
    end if
    if (present(pass_)) then
      pass = pass_
    else
      pass = all_pass
    end if
    dt = time_step_size * 0.5d0

    ! Compute RK1
    call space_operators_run(state(old), tend_rk(1), dt, pass)

    ! Compute RK2
    call update_state(dt, tend_rk(1), state(old), state(tmp))
    call space_operators_run(state(tmp), tend_rk(2), dt, pass)

    ! Compute RK3
    call update_state(dt, tend_rk(2), state(old), state(tmp))
    call space_operators_run(state(tmp), tend_rk(3), dt, pass)

    ! Compute RK4
    call update_state(time_step_size, tend_rk(3), state(old), state(tmp))
    call space_operators_run(state(tmp), tend_rk(4), dt, pass)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend_rk(0)%dgd(i,j) = tend_rk(1)%dgd(i,j) / 6.0d0 + tend_rk(2)%dgd(i,j) / 3.0d0 + tend_rk(3)%dgd(i,j) / 3.0d0 + tend_rk(4)%dgd(i,j) / 6.0d0
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend_rk(0)%du_c(i,j) = tend_rk(1)%du_c(i,j) / 6.0d0 + tend_rk(2)%du_c(i,j) / 3.0d0 + tend_rk(3)%du_c(i,j) / 3.0d0 + tend_rk(4)%du_c(i,j) / 6.0d0
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend_rk(0)%dv_c(i,j) = tend_rk(1)%dv_c(i,j) / 6.0d0 + tend_rk(2)%dv_c(i,j) / 3.0d0 + tend_rk(3)%dv_c(i,j) / 3.0d0 + tend_rk(4)%dv_c(i,j) / 6.0d0
      end do
    end do

    ip00 = inner_product(tend_rk(0), tend_rk(0))
    ip12 = inner_product(tend_rk(1), tend_rk(2))
    ip23 = inner_product(tend_rk(2), tend_rk(3))
    ip34 = inner_product(tend_rk(3), tend_rk(4))
    beta = 1.0d0 / (3.0d0 * ip00) * (ip12 + ip23 + ip34)
    call log_add_diag('beta', beta)

    dt = time_step_size * beta

    call update_state(dt, tend_rk(0), state(old), state(new))

  end subroutine runge_kutta

  subroutine leap_frog(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

  end subroutine leap_frog
  
end module dycore_mod
