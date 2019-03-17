module dycore_mod

  use params_mod, split_scheme_in => split_scheme, uv_adv_scheme_in => uv_adv_scheme
  use data_mod
  use log_mod
  use types_mod
  use mesh_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use parallel_mod
  use io_mod
  use diag_mod
  use history_mod
  use restart_mod
  use diffusion_mod
  use filter_mod
  use weno_mod

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
  integer, parameter :: csp1 = 1
  integer, parameter :: csp2 = 2
  integer, parameter :: isp  = 3
  integer, parameter :: all_pass = 0
  integer, parameter :: fast_pass = 1
  integer, parameter :: slow_pass = 2

  integer uv_adv_scheme
  integer, parameter :: center_diff = 0
  integer, parameter :: upwind = 1
  integer, parameter :: weno = 2

  interface
    subroutine integrator_interface(time_step_size, old, new, pass, qcon_modified_)
      real, intent(in) :: time_step_size
      integer, intent(in) :: old
      integer, intent(in) :: new
      integer, intent(in) :: pass
      logical, intent(in), optional :: qcon_modified_
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
    case default
      call log_error('Unknown time_scheme ' // trim(time_scheme) // '!')
    end select

    select case (split_scheme_in)
    case ('csp1')
      split_scheme = csp1
    case ('csp2')
      split_scheme = csp2
    case ('isp')
      split_scheme = isp
    case default
      split_scheme = 0
      call log_notice('No fast-slow split.')
    end select

    select case (uv_adv_scheme_in)
    case ('center_diff')
      uv_adv_scheme = center_diff
    case ('upwind')
      uv_adv_scheme = upwind
    case ('weno')
      uv_adv_scheme = weno
      call weno_init()
    case default
      call log_error('Unknown uv_adv_scheme ' // trim(uv_adv_scheme_in) // '!')
    end select

    call log_notice('Dycore module is initialized.')

  end subroutine dycore_init

  subroutine dycore_restart()

    call restart_read(state(old), static)

  end subroutine dycore_restart

  subroutine dycore_run()

    call reset_cos_lat_at_poles()

    call iap_transform(state(old))

    call diag_run(state(old))
    call output(state(old))
    call log_add_diag('total_mass', diag%total_mass)
    call log_add_diag('total_energy', diag%total_energy)
    call log_step()

    do while (.not. time_is_finished())
      tag = 0
      call time_integrate()
      call time_advance()
      call diag_run(state(old))
      call output(state(old))
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
    call weno_final()

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

  subroutine space_operators(state, tend, pass)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
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

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j) + tend%fv(i,j) - tend%u_pgf(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_full_lat(j, tend%du(:,j))
            s2 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
            tend%du(i1:i2,j) = tend%du(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j) - tend%fu(i,j) - tend%v_pgf(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
          s1 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_half_lat(j, tend%dv(:,j))
            s2 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
            tend%dv(i1:i2,j) = tend%dv(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%dgd(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_full_lat(j, tend%dgd(:,j))
            s2 = sum(tend%dgd(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
            tend%dgd(i1:i2,j) = tend%dgd(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_full_lat(j, tend%du(:,j))
            s2 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
            tend%du(i1:i2,j) = tend%du(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
          s1 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_half_lat(j, tend%dv(:,j))
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
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_full_lat(j, tend%dgd(:,j))
            s2 = sum(tend%dgd(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
            tend%dgd(i1:i2,j) = tend%dgd(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%fv(i,j) - tend%u_pgf(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_full_lat(j, tend%du(:,j))
            s2 = sum(tend%du(i1:i2,j) * state%iap%u(i1:i2,j))
            tend%du(i1:i2,j) = tend%du(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%fu(i,j) - tend%v_pgf(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
          s1 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
          if (abs(s1) > filter_inner_product_threshold) then
            call filter_array_at_half_lat(j, tend%dv(:,j))
            s2 = sum(tend%dv(i1:i2,j) * state%iap%v(i1:i2,j))
            tend%dv(i1:i2,j) = tend%dv(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
    end select

    ! tag = tag + 1
    ! if (time_is_alerted('hist0.output') .and. (tag == 3 .or. tag == 6)) then
    !  call history_write(state, tend, tag)
    ! end if

    ! call check_antisymmetry(tend, state)

  end subroutine space_operators

  subroutine zonal_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real u1, u2
    integer i, j

    select case (uv_adv_scheme)
    case (center_diff)
      ! U
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          u1 = state%u(i,j) + state%u(i-1,j)
          u2 = state%u(i,j) + state%u(i+1,j)
          tend%u_adv_lon(i,j) = 0.25 / coef%full_dlon(j) * (u2 * state%iap%u(i+1,j) - u1 * state%iap%u(i-1,j))
        end do
      end do
      ! V
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          u1 = state%u(i-1,j) + state%u(i-1,j+1)
          u2 = state%u(i,  j) + state%u(i,  j+1)
          tend%v_adv_lon(i,j) = 0.25 / coef%half_dlon(j) * (u2 * state%iap%v(i+1,j) - u1 * state%iap%v(i-1,j))
        end do
      end do
    case (upwind)
      ! U
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          u1 = state%u(i,j) + state%u(i-1,j)
          u2 = state%u(i,j) + state%u(i+1,j)
          tend%u_adv_lon(i,j) = 0.25 / coef%full_dlon(j) * &
            (u2 * (state%iap%u(i,j) + state%iap%u(i+1,j)) - uv_adv_upwind_lon_beta * abs(u2) * (state%iap%u(i+1,j) - state%iap%u(i,  j)) - &
             u1 * (state%iap%u(i,j) + state%iap%u(i-1,j)) + uv_adv_upwind_lon_beta * abs(u1) * (state%iap%u(i,  j) - state%iap%u(i-1,j)) - &
             (u2 - u1) * state%iap%u(i,j))
        end do
      end do
      ! V
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          u1 = state%u(i-1,j) + state%u(i-1,j+1)
          u2 = state%u(i,  j) + state%u(i,  j+1)
          tend%v_adv_lon(i,j) = 0.25 / coef%half_dlon(j) * &
            (u2 * (state%iap%v(i,j) + state%iap%v(i+1,j)) - uv_adv_upwind_lon_beta * abs(u2) * (state%iap%v(i+1,j) - state%iap%v(i,  j)) - &
             u1 * (state%iap%v(i,j) + state%iap%v(i-1,j)) + uv_adv_upwind_lon_beta * abs(u1) * (state%iap%v(i,  j) - state%iap%v(i-1,j)) - &
             (u2 - u1) * state%iap%v(i,j))
        end do
      end do
    case (weno)
      call weno_zonal_momentum_advection_operator(state, tend)
    end select

  end subroutine zonal_momentum_advection_operator

  subroutine meridional_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real v1, v2
    integer i, j

    select case (uv_adv_scheme)
    case (center_diff)
      ! U
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          v1 = (state%v(i,j-1) + state%v(i+1,j-1)) * mesh%half_cos_lat(j-1)
          v2 = (state%v(i,j  ) + state%v(i+1,j  )) * mesh%half_cos_lat(j  )
          tend%u_adv_lat(i,j) = 0.25 / coef%full_dlat(j) * (v2 * state%iap%u(i,j+1) - v1 * state%iap%u(i,j-1))
        end do
      end do
      ! V
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          v1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
          v2 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
          tend%v_adv_lat(i,j) = 0.25 / coef%half_dlat(j) * (v2 * state%iap%v(i,j+1) - v1 * state%iap%v(i,j-1))
        end do
      end do
    case (upwind)
      ! U
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          v1 = (state%v(i,j-1) + state%v(i+1,j-1)) * mesh%half_cos_lat(j-1)
          v2 = (state%v(i,j  ) + state%v(i+1,j  )) * mesh%half_cos_lat(j  )
          tend%u_adv_lat(i,j) = 0.25 / coef%full_dlat(j) * &
            (v2 * (state%iap%u(i,j) + state%iap%u(i,j+1)) - uv_adv_upwind_lat_beta * abs(v2) * (state%iap%u(i,j+1) - state%iap%u(i,j  )) - &
             v1 * (state%iap%u(i,j) + state%iap%u(i,j-1)) + uv_adv_upwind_lat_beta * abs(v1) * (state%iap%u(i,j  ) - state%iap%u(i,j-1)) - &
             (v2 - v1) * state%iap%u(i,j))
        end do
      end do
      ! V
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          v1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
          v2 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
          tend%v_adv_lat(i,j) = 0.25 / coef%half_dlat(j) * &
            (v2 * (state%iap%v(i,j) + state%iap%v(i,j+1)) - uv_adv_upwind_lat_beta * abs(v2) * (state%iap%v(i,j+1) - state%iap%v(i,j  )) - &
             v1 * (state%iap%v(i,j) + state%iap%v(i,j-1)) + uv_adv_upwind_lat_beta * abs(v1) * (state%iap%v(i,j  ) - state%iap%v(i,j-1)) - &
             (v2 - v1) * state%iap%v(i,j))
        end do
      end do
    case (weno)
      call weno_meridional_momentum_advection_operator(state, tend)
    end select

  end subroutine meridional_momentum_advection_operator

  subroutine coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real c1, c2
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      c1 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%fv(i,j) = 0.25 * (coef%full_f(j) + coef%full_c(j) * state%u(i,j)) * &
          (c1 * (state%iap%v(i,j-1) + state%iap%v(i+1,j-1)) + &
           c2 * (state%iap%v(i,j  ) + state%iap%v(i+1,j  )))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%fu(i,j) = 0.25 * &
          ((coef%full_f(j  ) + coef%full_c(j  ) * state%u(i,  j  )) * state%iap%u(i,  j  ) + &
           (coef%full_f(j  ) + coef%full_c(j  ) * state%u(i-1,j  )) * state%iap%u(i-1,j  ) + &
           (coef%full_f(j+1) + coef%full_c(j+1) * state%u(i,  j+1)) * state%iap%u(i,  j+1) + &
           (coef%full_f(j+1) + coef%full_c(j+1) * state%u(i-1,j+1)) * state%iap%u(i-1,j+1))
      end do
    end do

  end subroutine coriolis_operator

  subroutine zonal_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_pgf(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) / coef%full_dlon(j) * &
          (state%gd(i+1,j) + static%ghs(i+1,j) - state%gd(i,j) - static%ghs(i,j))
      end do
    end do

  end subroutine zonal_pressure_gradient_force_operator

  subroutine meridional_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) / coef%half_dlat(j) * mesh%half_cos_lat(j) * &
          (state%gd(i,j+1) + static%ghs(i,j+1) - state%gd(i,j) - static%ghs(i,j))
      end do
    end do

  end subroutine meridional_pressure_gradient_force_operator

  subroutine zonal_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lon(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%iap%u(i,  j) - &
                                  (state%iap%gd(i,j) + state%iap%gd(i-1,j)) * state%iap%u(i-1,j)) &
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
        tend%mass_div_lat(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v(i,j  ) * mesh%half_cos_lat(j  ) - &
                                  (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v(i,j-1) * mesh%half_cos_lat(j-1)) &
                                 * 0.5 / coef%full_dlat(j)
      end do
    end do

    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v(i,j)
      end do
      call parallel_zonal_sum(sp)
      sp = sp * 2.0 / mesh%num_full_lon / radius / mesh%dlat
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sp
      end do
    end if

    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v(i,j-1)
      end do
      call parallel_zonal_sum(np)
      np = np * 2.0 / mesh%num_full_lon / radius / mesh%dlat
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

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
      end do
    end do

    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)

    do j = parallel%full_lat_lb, parallel%full_lat_ub
      do i = parallel%full_lon_lb, parallel%full_lon_ub
        new_state%iap%gd(i,j) = sqrt(new_state%gd(i,j))
      end do
    end do

    ! Update IAP wind state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%iap%u(i,j) = old_state%iap%u(i,j) + dt * tend%du(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%v(i,j) = old_state%iap%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    ! Transform from IAP to normal state.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = new_state%iap%u(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v(i,j) = new_state%iap%v(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i,j+1))
      end do
    end do

    call parallel_fill_halo(new_state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v(:,:), all_halo=.true.)

  end subroutine update_state

  subroutine time_integrate()

    select case (split_scheme)
    case (csp2)
      call csp2_splitting()
    case (isp)
      call isp_splitting()
    case default
      call integrator(time_step_size, old, new, all_pass)
    end select

    if (use_diffusion) then
      call ordinary_diffusion(time_step_size, state(new))
    end if

  end subroutine time_integrate

  subroutine csp2_splitting()

    real fast_dt
    integer subcycle, t1, t2

    fast_dt = time_step_size / subcycles
    t1 = 0
    t2 = old

    call integrator(0.5 * time_step_size, old, t1, slow_pass)
    do subcycle = 1, subcycles
      call integrator(fast_dt, t1, t2, fast_pass)
      call time_swap_indices(t1, t2)
    end do
    call integrator(0.5 * time_step_size, t1, new, slow_pass)

  end subroutine csp2_splitting

  subroutine isp_splitting()

    integer subcycle, t1, t2
    real fast_dt, half_dt, ip1, ip2, beta

    type(state_type), pointer :: saved_state
    type(tend_type), pointer :: slow_tend
    type(tend_type), pointer :: accum_fast_tend

    fast_dt = time_step_size / subcycles
    half_dt = time_step_size * 0.5
    t1 = old
    t2 = 0
    saved_state => state(-1)
    slow_tend => tend(-2)
    accum_fast_tend => tend(-1)

    call copy_data(state(old), saved_state)                            ! F^n
    call zero_data(accum_fast_tend)
    call space_operators(state(old), slow_tend, slow_pass)             ! L2(F^n)

    ! Fast pass
    do subcycle = 1, subcycles
      call space_operators(state(t1), tend(t1), fast_pass)             ! L1(P^{n,k})
      call add_data(slow_tend, tend(t1))                               ! L1(P^{n,k}) + L2(F^n)
      call update_state(fast_dt * 0.5, tend(t1), state(t1), state(t2)) ! P1^{n,k+1/2}

      call space_operators(state(t2), tend(t1), fast_pass)             ! L1(P1^{n,k+1/2})
      call add_data(slow_tend, tend(t1))                               ! L1(P1^{n,k+1/2}) + L2(F^n)
      call update_state(fast_dt * 0.5, tend(t1), state(t1), state(t2)) ! P2^{n,k+1/2}

      call space_operators(state(t2), tend(t2), fast_pass)             ! L1(P2^{n,k+1/2})
      call add_data(tend(t2), accum_fast_tend)                         ! sum(L1(P2^{n,k+1/2}))
      call add_data(slow_tend, tend(t2))                               ! L1(P2^{n,k+1/2}) + L2(F^n)
      call update_state(fast_dt, tend(t2), state(t1), state(t2))       ! P^{n,k+1}

      call time_swap_indices(t1, t2)
    end do
    call scale_data(2.0 / subcycles, accum_fast_tend)                  ! 2 / m sum(L1(P2^{n,k+1/2}))

    ! Slow pass
    call space_operators(state(t1), tend(t1), slow_pass)               ! L2(Q^n)
    call sub_data(slow_tend, tend(t1))                                 ! L2(Q^n) - L2(F^n)
    call update_state(half_dt, tend(t1), state(t1), state(new))        ! Q1^{n+1}

    call space_operators(state(new), tend(t1), slow_pass)              ! L2(Q1^{n+1})
    call sub_data(slow_tend, tend(t1))                                 ! L2(Q1^{n+1}) - L2(F^n)
    call update_state(half_dt, tend(t1), state(t1), state(new))        ! Q2^{n+1}

    call space_operators(state(new), tend(new), slow_pass)             ! L2(Q2^{n+1})

    call add_data(slow_tend, tend(new))                                ! L2(Q2^{n+1}) + L2(F^n)
    call add_data(accum_fast_tend, tend(new))                          ! L2(Q2^{n+1}) + L2(F^n) + 2 / m sum(L1(P2^{n,k+1/2})) -> R2^n 

    ip1 = inner_product(tend(new), saved_state)                        ! (R2^n, F^n)
    ip2 = inner_product(tend(new), tend(new))                          ! (R2^n, R2^n)
    beta = merge(ip1 / ip2, 1.0, qcon_modified .and. ip1 /= 0.0 .and. ip2 /= 0.0)
    beta = beta * 4.0 / time_step_size
    call log_add_diag('beta', beta)

    half_dt = half_dt * beta
    call update_state(half_dt, tend(new), saved_state, state(new))

  end subroutine isp_splitting

  subroutine predict_correct(time_step_size, old, new, pass, qcon_modified_)

    real, intent(in) :: time_step_size
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass
    logical, intent(in), optional :: qcon_modified_

    real dt, ip1, ip2, beta
    logical flag

    if (present(qcon_modified_)) then
      flag = qcon_modified_
    else
      flag = qcon_modified
    end if

    dt = time_step_size * 0.5

    ! Do first predict step.
    call space_operators(state(old), tend(old), pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do second predict step.
    call space_operators(state(new), tend(old), pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do correct step.
    call space_operators(state(new), tend(new), pass)

    ip1 = inner_product(tend(old), tend(new))
    ip2 = inner_product(tend(new), tend(new))
    beta = merge(ip1 / ip2, 1.0, flag .and. ip1 /= 0.0 .and. ip2 /= 0.0)
    call log_add_diag('beta', beta)

    dt = time_step_size * beta
    call update_state(dt, tend(new), state(old), state(new))

  end subroutine predict_correct

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
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ip_u_adv_lon = ip_u_adv_lon + tend%u_adv_lon(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_u_adv_lat = ip_u_adv_lat + tend%u_adv_lat(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_fv = ip_fv + tend%fv(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_u_pgf = ip_u_pgf + tend%u_pgf(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_v_adv_lon = ip_v_adv_lon + tend%v_adv_lon(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_v_adv_lat = ip_v_adv_lat + tend%v_adv_lat(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_fu = ip_fu + tend%fu(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_v_pgf = ip_v_pgf + tend%v_pgf(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
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
