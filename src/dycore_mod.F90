module dycore_mod

  use params_mod, split_scheme_in => split_scheme
  use data_mod
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

  implicit none

  private

  public dycore_init
  public dycore_restart
  public dycore_run
  public dycore_final

  ! 1: csp1: first order conservative split
  ! 2: csp2: second order conservative split
  ! 3: isp: inproved second order conservative split
  integer split_scheme
  integer, parameter :: all_pass = 0
  integer, parameter :: fast_pass = 1
  integer, parameter :: slow_pass = 2

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
    call diffusion_init()
    call filter_init()

    select case (time_scheme)
    case ('predict_correct')
      integrator => predict_correct
    case ('runge_kutta')
      integrator => runge_kutta
    case default
      call log_error('Unknown time_scheme ' // trim(time_scheme) // '!')
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

    call log_notice('Dycore module is initialized.')

  end subroutine dycore_init

  subroutine dycore_restart()

    call restart_read(state(old_time_idx), static)

  end subroutine dycore_restart

  subroutine dycore_run()

    call reset_cos_lat_at_poles()

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
    call diffusion_final()
    call filter_final()

    call log_notice('Dycore module is finalized.')

  end subroutine dycore_final

  subroutine reset_cos_lat_at_poles()

    integer j

    j = parallel%full_lat_south_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_south_pole_idx) * 0.25
    coef%full_dlon(j) = radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = radius * mesh%dlat * mesh%full_cos_lat(j)

    j = parallel%full_lat_north_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_north_pole_idx) * 0.25
    coef%full_dlon(j) = radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = radius * mesh%dlat * mesh%full_cos_lat(j)

  end subroutine reset_cos_lat_at_poles

  subroutine output(state)

    type(state_type), intent(in) :: state

    if (time_is_alerted('hist0.output')) call history_write(state, static, diag)
    if (time_is_alerted('restart.output')) call restart_write(state, static)

  end subroutine output

  subroutine space_operators(state, tend, dt, pass)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real, intent(in) :: dt
    integer, intent(in) :: pass

    integer i, j, i1, i2
    real s1, s2

    ! Allow me to use i1 and i2 as shorthands.
    i1 = parallel%full_lon_start_idx
    i2 = parallel%full_lon_end_idx

    select case (pass)
    case (all_pass)
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)
      call coriolis_operator(state, tend)
      call zonal_pressure_gradient_force_operator(state, tend)
      call meridional_pressure_gradient_force_operator(state, tend)
      call zonal_mass_divergence_operator(state, tend)
      call meridional_mass_divergence_operator(state, tend)

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j) + tend%fv(i,j) - tend%u_pgf(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j) - tend%fu(i,j) - tend%v_pgf(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
    case (slow_pass)
#ifndef NDEBUG
      tend%fv = 0.0
      tend%u_pgf = 0.0
      tend%fu = 0.0
      tend%v_pgf = 0.0
      tend%mass_div_lon = 0.0
      tend%mass_div_lat = 0.0
#endif
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du(:,j))
            s2 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
            tend%du(i1:i2,j) = tend%du(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dv(:,j))
            s2 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
            tend%dv(i1:i2,j) = tend%dv(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      tend%dgd = 0.0
    case (fast_pass)
#ifndef NDEBUG
      tend%u_adv_lon = 0.0
      tend%u_adv_lat = 0.0
      tend%v_adv_lon = 0.0
      tend%v_adv_lat = 0.0
#endif
      call coriolis_operator(state, tend)
      call zonal_pressure_gradient_force_operator(state, tend)
      call meridional_pressure_gradient_force_operator(state, tend)
      call zonal_mass_divergence_operator(state, tend)
      call meridional_mass_divergence_operator(state, tend)

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%dgd(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dgd(:,j))
            s2 = sum(tend%dgd(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
            tend%dgd(i1:i2,j) = tend%dgd(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du(i,j) = tend%fv(i,j) - tend%u_pgf(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du(:,j))
            s2 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
            tend%du(i1:i2,j) = tend%du(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%fu(i,j) - tend%v_pgf(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dv(:,j))
            s2 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
            tend%dv(i1:i2,j) = tend%dv(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
    end select

    tag = tag + 1
    if (time_is_alerted('hist0.output') .and. (tag == 3 .or. tag == 6)) then
     call history_write(tend, tag)
    end if

    ! call check_antisymmetry(tend, state)

  end subroutine space_operators

  subroutine zonal_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_adv_lon(i,j) = 0.25 / coef%full_dlon(j) * &
          ((state%u(i,j) + state%u(i+1,j)) * state%iap%u(i+1,j) - &
           (state%u(i,j) + state%u(i-1,j)) * state%iap%u(i-1,j))
      end do
    end do

    ! V
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lon(i,j) = 0.25 / coef%full_dlon(j) * &
          ((state%u(i,j) + state%u(i+1,j)) * state%iap%v(i+1,j) - &
           (state%u(i,j) + state%u(i-1,j)) * state%iap%v(i-1,j))
      end do
    end do

  end subroutine zonal_momentum_advection_operator

  subroutine meridional_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_adv_lat(i,j) = 0.25 / coef%full_dlat(j) * &
          ((state%v(i,j) * mesh%full_cos_lat(j) + state%v(i,j+1) * mesh%full_cos_lat(j+1)) * state%iap%u(i,j+1) - &
           (state%v(i,j) * mesh%full_cos_lat(j) + state%v(i,j-1) * mesh%full_cos_lat(j-1)) * state%iap%u(i,j-1))
                              
      end do
    end do

    ! V
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat(i,j) = 0.25 / coef%full_dlat(j) * &
          ((state%v(i,j) * mesh%full_cos_lat(j) + state%v(i,j+1) * mesh%full_cos_lat(j+1)) * state%iap%v(i,j+1) - &
           (state%v(i,j) * mesh%full_cos_lat(j) + state%v(i,j-1) * mesh%full_cos_lat(j-1)) * state%iap%v(i,j-1))
      end do
    end do

  end subroutine meridional_momentum_advection_operator

  subroutine coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

#ifdef AVERAGE_CORIOLIS
    real c1, c2
#endif
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
#ifdef AVERAGE_CORIOLIS
      c1 = mesh%full_cos_lat(j-1) / mesh%full_cos_lat(j)
      c2 = mesh%full_cos_lat(j+1) / mesh%full_cos_lat(j)
#endif
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
#ifdef AVERAGE_CORIOLIS
        tend%fv(i,j) = 0.25 * (coef%full_f(j) + coef%full_c(j) * state%u(i,j)) * &
          (c1 * (state%iap%v(i-1,j-1) + state%iap%v(i+1,j-1)) + &
           c2 * (state%iap%v(i-1,j+1) + state%iap%v(i+1,j+1)))
#else
        tend%fv(i,j) = (coef%full_f(j) + coef%full_c(j) * state%u(i,j)) * state%iap%v(i,j)
#endif
      end do
    end do

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
#ifdef AVERAGE_CORIOLIS
        tend%fu(i,j) = 0.25 * ((coef%full_f(j-1) + coef%full_c(j-1) * state%u(i-1,j-1)) * state%iap%u(i-1,j-1) + &
                               (coef%full_f(j-1) + coef%full_c(j-1) * state%u(i+1,j-1)) * state%iap%u(i+1,j-1) + &
                               (coef%full_f(j+1) + coef%full_c(j+1) * state%u(i-1,j+1)) * state%iap%u(i-1,j+1) + &
                               (coef%full_f(j+1) + coef%full_c(j+1) * state%u(i+1,j+1)) * state%iap%u(i+1,j+1))
#else
        tend%fu(i,j) = (coef%full_f(j) + coef%full_c(j) * state%u(i,j)) * state%iap%u(i,j)
#endif
      end do
    end do

  end subroutine coriolis_operator

  subroutine zonal_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_pgf(i,j) = 0.5 * state%iap%gd(i,j) / coef%full_dlon(j) * &
          (state%gd(i+1,j) + static%ghs(i+1,j) - state%gd(i-1,j) - static%ghs(i-1,j))
      end do
    end do

  end subroutine zonal_pressure_gradient_force_operator

  subroutine meridional_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = 0.5 * state%iap%gd(i,j) / coef%full_dlat(j) * mesh%full_cos_lat(j) * &
          (state%gd(i,j+1) + static%ghs(i,j+1) - state%gd(i,j-1) - static%ghs(i,j-1))
      end do
    end do

  end subroutine meridional_pressure_gradient_force_operator

  subroutine zonal_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lon(i,j) = (state%iap%gd(i+1,j) * state%iap%u(i+1,j) - &
                                  state%iap%gd(i-1,j) * state%iap%u(i-1,j)) &
                                 * 0.5 / coef%full_dlon(j)
      end do
    end do

  end subroutine zonal_mass_divergence_operator

  subroutine meridional_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real sp, np
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = (state%iap%gd(i,j+1) * state%iap%v(i,j+1) * mesh%full_cos_lat(j+1) - &
                                  state%iap%gd(i,j-1) * state%iap%v(i,j-1) * mesh%full_cos_lat(j-1)) &
                                 * 0.5 / coef%full_dlat(j)
      end do
    end do

    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + state%iap%gd(i,j+1) * state%iap%v(i,j+1) * mesh%full_cos_lat(j+1)
      end do
      call parallel_zonal_sum(sp)
      sp = sp / mesh%num_full_lon * 0.5 / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sp
      end do
    end if

    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - state%iap%gd(i,j-1) * state%iap%v(i,j-1) * mesh%full_cos_lat(j-1)
      end do
      call parallel_zonal_sum(np)
      np = np / mesh%num_full_lon * 0.5 / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = np
      end do
    end if

  end subroutine meridional_mass_divergence_operator

  subroutine update_state(dt, tend, old_state, new_state)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j

    ! Update IAP state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
        new_state%iap%u(i,j) = old_state%iap%u(i,j) + dt * tend%du(i,j)
        new_state%iap%v(i,j) = old_state%iap%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v(:,:), all_halo=.true.)

    ! Transform from IAP to normal state.
    do j = parallel%full_lat_lb, parallel%full_lat_ub
      do i = parallel%full_lon_lb, parallel%full_lon_ub
        new_state%iap%gd(i,j) = sqrt(new_state%gd(i,j))
        new_state%u(i,j) = new_state%iap%u(i,j) / new_state%iap%gd(i,j)
        new_state%v(i,j) = new_state%iap%v(i,j) / new_state%iap%gd(i,j)
      end do
    end do

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%du(i,j) * tend2%du(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dv(i,j) * tend2%dv(i,j) * mesh%full_cos_lat(j)
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
      call ordinary_diffusion(time_step_size, state(new_time_idx))
    end if

  end subroutine time_integrate

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
    call space_operators(state(old), tend(old), dt, pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do second predict step.
    call space_operators(state(new), tend(old), dt, pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do correct step.
    call space_operators(state(new), tend(new), dt, pass)

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
    call space_operators(state(old), tend_rk(1), dt, pass)

    ! Compute RK2
    call update_state(dt, tend_rk(1), state(old), state(tmp))
    call space_operators(state(tmp), tend_rk(2), dt, pass)

    ! Compute RK3
    call update_state(dt, tend_rk(2), state(old), state(tmp))
    call space_operators(state(tmp), tend_rk(3), dt, pass)

    ! Compute RK4
    call update_state(time_step_size, tend_rk(3), state(old), state(tmp))
    call space_operators(state(tmp), tend_rk(4), dt, pass)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend_rk(0)%dgd(i,j) = tend_rk(1)%dgd(i,j) / 6.0d0 + tend_rk(2)%dgd(i,j) / 3.0d0 + tend_rk(3)%dgd(i,j) / 3.0d0 + tend_rk(4)%dgd(i,j) / 6.0d0
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend_rk(0)%du(i,j) = tend_rk(1)%du(i,j) / 6.0d0 + tend_rk(2)%du(i,j) / 3.0d0 + tend_rk(3)%du(i,j) / 3.0d0 + tend_rk(4)%du(i,j) / 6.0d0
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend_rk(0)%dv(i,j) = tend_rk(1)%dv(i,j) / 6.0d0 + tend_rk(2)%dv(i,j) / 3.0d0 + tend_rk(3)%dv(i,j) / 3.0d0 + tend_rk(4)%dv(i,j) / 6.0d0
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

  subroutine check_antisymmetry(tend, state)

    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: state

    integer i, j
    real ip_u_adv_lon
    real ip_u_adv_lat
    real ip_fv
    real ip_u_pgf
    real ip_v_adv_lon
    real ip_v_adv_lat
    real ip_fu
    real ip_v_pgf
    real ip_mass_div_lon
    real ip_mass_div_lat

    ip_u_adv_lon = 0.0
    ip_u_adv_lat = 0.0
    ip_fv = 0.0
    ip_u_pgf = 0.0
    ip_v_adv_lon = 0.0
    ip_v_adv_lat = 0.0
    ip_fu = 0.0
    ip_v_pgf = 0.0
    ip_mass_div_lon = 0.0
    ip_mass_div_lat = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_u_adv_lon = ip_u_adv_lon + tend%u_adv_lon(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_u_adv_lat = ip_u_adv_lat + tend%u_adv_lat(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_fv = ip_fv + tend%fv(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_u_pgf = ip_u_pgf + tend%u_pgf(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_v_adv_lon = ip_v_adv_lon + tend%v_adv_lon(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_v_adv_lat = ip_v_adv_lat + tend%v_adv_lat(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_fu = ip_fu + tend%fu(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_v_pgf = ip_v_pgf + tend%v_pgf(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_mass_div_lon = ip_mass_div_lon + tend%mass_div_lon(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%full_cos_lat(j)
        ip_mass_div_lat = ip_mass_div_lat + tend%mass_div_lat(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%full_cos_lat(j)
      end do
    end do

    print *, &
      ip_u_adv_lon + ip_v_adv_lon + ip_u_adv_lat + ip_v_adv_lat, &
      ip_u_pgf + ip_mass_div_lon, &
      ip_fv - ip_fu, &
      ip_v_pgf + ip_mass_div_lat

  end subroutine check_antisymmetry
  
end module dycore_mod
