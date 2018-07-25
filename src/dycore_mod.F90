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
    call reduce_init()
    call diffusion_init()

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

    call log_notice('Dycore module is finalized.')

  end subroutine dycore_final

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

    integer i, j

    call reduce_run(state, static)

    select case (pass)
    case (all_pass)
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)
      call coriolis_operator(state, tend)
      call curvature_operator(state, tend)
      call zonal_pressure_gradient_force_operator(state, tend)
      call meridional_pressure_gradient_force_operator(state, tend)
      call zonal_mass_divergence_operator(state, tend)
      call meridional_mass_divergence_operator(state, tend)

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j) + tend%fv(i,j) + tend%cv(i,j) - tend%u_pgf(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j) - tend%fu(i,j) - tend%cu(i,j) - tend%v_pgf(i,j)
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
    case (slow_pass)
#ifndef NDEBUG
      tend%fv = 0.0
      tend%cv = 0.0
      tend%u_pgf = 0.0
      tend%fu = 0.0
      tend%cu = 0.0
      tend%v_pgf = 0.0
      tend%mass_div_lon = 0.0
      tend%mass_div_lat = 0.0
#endif
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j)
        end do
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
      call curvature_operator(state, tend)
      call zonal_pressure_gradient_force_operator(state, tend)
      call meridional_pressure_gradient_force_operator(state, tend)
      call zonal_mass_divergence_operator(state, tend)
      call meridional_mass_divergence_operator(state, tend)

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%fv(i,j) + tend%cv(i,j) - tend%u_pgf(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%fu(i,j) - tend%cu(i,j) - tend%v_pgf(i,j)
        end do
      end do
    end select

    tag = tag + 1
    if (time_is_alerted('hist0.output') .and. (tag == 3 .or. tag == 6 .or. tag == 18 .or. tag == 21)) then
      call history_write(tend, tag)
    end if

    ! call check_antisymmetry(tend, state)

  end subroutine space_operators

  subroutine zonal_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    integer i, j, k

    ! U
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (half_reduce_factor(j) == 1 .or. .not. reduce_adv_lon) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%u_adv_lon(i,j) = 0.5 / coef%half_dlon(j) * &
            ((state%u(i,j) + state%u(i+1,j)) * state%iap%u(i+1,j) - &
             (state%u(i,j) + state%u(i-1,j)) * state%iap%u(i-1,j))
        end do
      else
        tend%u_adv_lon(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = 0.5 / coef%half_dlon(j) / half_reduce_factor(j) * &
              ((half_reduced_state(j)%u(i,k,2) + half_reduced_state(j)%u(i+1,k,2)) * half_reduced_state(j)%iap%u(i+1,k,2) - &
               (half_reduced_state(j)%u(i,k,2) + half_reduced_state(j)%u(i-1,k,2)) * half_reduced_state(j)%iap%u(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%u_adv_lon(:,j))
        end do
        call parallel_overlay_inner_halo(tend%u_adv_lon(:,j), left_halo=.true.)
      end if
    end do

    ! V
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) == 1 .or. .not. reduce_adv_lon) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%v_adv_lon(i,j) = 0.5 / coef%full_dlon(j) * &
            ((state%u(i,  j) + state%u(i,  j-1)) * state%iap%v(i+1,j) - &
             (state%u(i-1,j) + state%u(i-1,j-1)) * state%iap%v(i-1,j))
        end do
      else
        tend%v_adv_lon(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.5 / coef%full_dlon(j) / full_reduce_factor(j) * &
              ((full_reduced_state(j)%u(i,  k,2) + full_reduced_state(j)%u(i,  k,1)) * full_reduced_state(j)%iap%v(i+1,k,2) - &
               (full_reduced_state(j)%u(i-1,k,2) + full_reduced_state(j)%u(i-1,k,1)) * full_reduced_state(j)%iap%v(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%v_adv_lon(:,j))
        end do
        call parallel_overlay_inner_halo(tend%v_adv_lon(:,j), left_halo=.true.)
      end if
    end do

  end subroutine zonal_momentum_advection_operator

  subroutine meridional_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_adv_lat(i,j) = (state%iap%u(i,j+1) * (state%v(i,j+1) + state%v(i+1,j+1)) * mesh%full_cos_lat(j+1) - &
                               state%iap%u(i,j-1) * (state%v(i,j  ) + state%v(i+1,j  )) * mesh%full_cos_lat(j  )) * &
                              0.5 / coef%half_dlat(j)
      end do
    end do

    ! V
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat(i,j) = (state%iap%v(i,j+1) * (state%v(i,j) * mesh%full_cos_lat(j) + state%v(i,j+1) * mesh%full_cos_lat(j+1)) - &
                               state%iap%v(i,j-1) * (state%v(i,j) * mesh%full_cos_lat(j) + state%v(i,j-1) * mesh%full_cos_lat(j-1))) * &
                              0.5 / coef%full_dlat(j)
      end do
    end do

  end subroutine meridional_momentum_advection_operator

  subroutine coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    real c1, c2
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c1, c2, i, k, reduced_tend)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      c1 = mesh%full_cos_lat(j  ) / mesh%half_cos_lat(j)
      c2 = mesh%full_cos_lat(j+1) / mesh%half_cos_lat(j)
      if (half_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%fv(i,j) = 0.25 * coef%cori(j) * &
            (c1 * (state%iap%v(i,j  ) + state%iap%v(i+1,j  )) + &
             c2 * (state%iap%v(i,j+1) + state%iap%v(i+1,j+1)))
        end do
      else
        tend%fv(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = 0.25 * coef%cori(j) * &
              (c1 * (half_reduced_state(j)%iap%v(i,k,2) + half_reduced_state(j)%iap%v(i+1,k,2)) + &
               c2 * (half_reduced_state(j)%iap%v(i,k,3) + half_reduced_state(j)%iap%v(i+1,k,3)))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%fv(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fv(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      tend%fu(:,j) = 0.0
      if (half_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%fu(i,j) = 0.25 * coef%cori(j) * (state%iap%u(i,j) + state%iap%u(i-1,j))
        end do
      else
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = 0.25 * coef%cori(j) * (half_reduced_state(j)%iap%u(i,k,2) + half_reduced_state(j)%iap%u(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%fu(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fu(:,j), left_halo=.true.)
      end if
      if (half_reduce_factor(j-1) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%fu(i,j) = tend%fu(i,j) + 0.25 * coef%cori(j-1) * (state%iap%u(i,j-1) + state%iap%u(i-1,j-1))
        end do
      else
        ! Clear out right halo of tendency, because we will overlay them with left inner halo below.
        tend%fu(parallel%full_lon_end_idx+1:,j) = 0.0
        do k = 1, half_reduce_factor(j-1)
          do i = reduced_start_idx_at_half_lat(j-1), reduced_end_idx_at_half_lat(j-1)
            reduced_tend(i) = 0.25 * coef%cori(j-1) * (half_reduced_state(j-1)%iap%u(i,k,2) + half_reduced_state(j-1)%iap%u(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j-1, k, reduced_tend, tend%fu(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fu(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine coriolis_operator

  subroutine curvature_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    real c1, c2
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c1, c2, i, k, reduced_tend)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      c1 = mesh%full_cos_lat(j  ) / mesh%half_cos_lat(j)
      c2 = mesh%full_cos_lat(j+1) / mesh%half_cos_lat(j)
      if (half_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%cv(i,j) = 0.25 * coef%curv(j) * state%u(i,j) * &
            (c1 * (state%iap%v(i,j  ) + state%iap%v(i+1,j  )) + &
             c2 * (state%iap%v(i,j+1) + state%iap%v(i+1,j+1)))
        end do
      else
        tend%cv(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = 0.25 * coef%curv(j) * half_reduced_state(j)%u(i,k,2) * &
              (c1 * (half_reduced_state(j)%iap%v(i,k,2) + half_reduced_state(j)%iap%v(i+1,k,2)) + &
               c2 * (half_reduced_state(j)%iap%v(i,k,3) + half_reduced_state(j)%iap%v(i+1,k,3)))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%cv(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cv(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      tend%cu(:,j) = 0.0
      if (half_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%cu(i,j) = 0.25 * coef%curv(j) * (state%u(i,  j) * state%iap%u(i,  j) + &
                                                state%u(i-1,j) * state%iap%u(i-1,j))
        end do
      else
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = 0.25 * coef%curv(j) * &
              (half_reduced_state(j)%u(i,  k,2) * half_reduced_state(j)%iap%u(i,  k,2) + &
               half_reduced_state(j)%u(i-1,k,2) * half_reduced_state(j)%iap%u(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%cu(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cu(:,j), left_halo=.true.)
      end if
      if (half_reduce_factor(j-1) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%cu(i,j) = tend%cu(i,j) + 0.25 * coef%curv(j-1) * (state%u(i,  j-1) * state%iap%u(i,  j-1) + &
                                                                 state%u(i-1,j-1) * state%iap%u(i-1,j-1))
        end do
      else
        ! Clear out right halo of tendency, because we will overlay them with left inner halo below.
        tend%cu(parallel%full_lon_end_idx+1:,j) = 0.0
        do k = 1, half_reduce_factor(j-1)
          do i = reduced_start_idx_at_half_lat(j-1), reduced_end_idx_at_half_lat(j-1)
            reduced_tend(i) = 0.25 * coef%curv(j-1) * &
              (half_reduced_state(j-1)%u(i,  k,2) * half_reduced_state(j-1)%iap%u(i,  k,2) + &
               half_reduced_state(j-1)%u(i-1,k,2) * half_reduced_state(j-1)%iap%u(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j-1, k, reduced_tend, tend%cu(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cu(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine curvature_operator

  subroutine zonal_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (half_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%u_pgf(i,j) = (state%iap%gd(i,j) + state%iap%gd(i+1,j)) / coef%half_dlon(j) * &
            (state%gd(i+1,j) - state%gd(i,j) + static%ghs(i+1,j) - static%ghs(i,j))
        end do
      else
        tend%u_pgf(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = (half_reduced_state(j)%iap%gd(i,k,2) + half_reduced_state(j)%iap%gd(i+1,k,2)) * &
              (half_reduced_state(j)%gd(i+1,k,2) - half_reduced_state(j)%gd(i,k,2) + &
               half_reduced_static(j)%ghs(i+1,k,2) - half_reduced_static(j)%ghs(i,k,2)) / &
              coef%half_dlon(j) / half_reduce_factor(j)
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%u_pgf(:,j))
        end do
        call parallel_overlay_inner_halo(tend%u_pgf(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine zonal_pressure_gradient_force_operator

  subroutine meridional_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * &
          (state%gd(i,j) + static%ghs(i,j) - state%gd(i,j-1) - static%ghs(i,j-1)) * &
          mesh%full_cos_lat(j) / coef%full_dlat(j)
      end do
    end do

  end subroutine meridional_pressure_gradient_force_operator

  subroutine zonal_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (half_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%mass_div_lon(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%iap%u(i,  j) - &
                                    (state%iap%gd(i,j) + state%iap%gd(i-1,j)) * state%iap%u(i-1,j)) &
                                   / coef%half_dlon(j)
        end do
      else
        tend%mass_div_lon(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = ((half_reduced_state(j)%iap%gd(i,k,2) + half_reduced_state(j)%iap%gd(i+1,k,2)) * half_reduced_state(j)%iap%u(i,  k,2) - &
                               (half_reduced_state(j)%iap%gd(i,k,2) + half_reduced_state(j)%iap%gd(i-1,k,2)) * half_reduced_state(j)%iap%u(i-1,k,2)) / &
                              coef%half_dlon(j) / half_reduce_factor(j)
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%mass_div_lon(:,j))
        end do
        call parallel_overlay_inner_halo(tend%mass_div_lon(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine zonal_mass_divergence_operator

  subroutine meridional_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real sp, np
    integer i, j

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v(i,j+1) * mesh%full_cos_lat(j+1) - &
                                  (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v(i,j  ) * mesh%full_cos_lat(j  )) / &
                                 coef%half_dlat(j)
      end do
    end do

  end subroutine meridional_mass_divergence_operator

  subroutine update_state(dt, tend, old_state, new_state)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
        new_state%iap%gd(i,j) = sqrt(new_state%gd(i,j))
      end do
    end do

    call parallel_fill_halo(new_state%iap%gd(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)

    ! Update IAP wind state.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%iap%u(i,j) = old_state%iap%u(i,j) + dt * tend%du(i,j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%v(i,j) = old_state%iap%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    ! Transform from IAP to normal state.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = new_state%iap%u(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v(i,j) = new_state%iap%v(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i,j-1))
      end do
    end do

    call parallel_fill_halo(new_state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v(:,:), all_halo=.true.)

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend1%du(i,j) * tend2%du(i,j) * mesh%half_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dv(i,j) * tend2%dv(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dgd(i,j) * tend2%dgd(i,j) * mesh%half_cos_lat(j)
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
      call space_operators(state(half), tend(old), dt, pass)
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

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend_rk(0)%dgd(i,j) = tend_rk(1)%dgd(i,j) / 6.0d0 + tend_rk(2)%dgd(i,j) / 3.0d0 + tend_rk(3)%dgd(i,j) / 3.0d0 + tend_rk(4)%dgd(i,j) / 6.0d0
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
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

  subroutine leap_frog(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

  end subroutine leap_frog

  subroutine check_antisymmetry(tend, state)

    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: state

    integer i, j
    real ip_u_adv_lon
    real ip_u_adv_lat
    real ip_fv
    real ip_cv
    real ip_u_pgf
    real ip_v_adv_lon
    real ip_v_adv_lat
    real ip_fu
    real ip_cu
    real ip_v_pgf
    real ip_mass_div_lon
    real ip_mass_div_lat

    ip_u_adv_lon = 0.0
    ip_u_adv_lat = 0.0
    ip_fv = 0.0
    ip_cv = 0.0
    ip_u_pgf = 0.0
    ip_v_adv_lon = 0.0
    ip_v_adv_lat = 0.0
    ip_fu = 0.0
    ip_cu = 0.0
    ip_v_pgf = 0.0
    ip_mass_div_lon = 0.0
    ip_mass_div_lat = 0.0

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ip_u_adv_lon = ip_u_adv_lon + tend%u_adv_lon(i,j) * state%iap%u(i,j) * mesh%half_cos_lat(j)
        ip_u_adv_lat = ip_u_adv_lat + tend%u_adv_lat(i,j) * state%iap%u(i,j) * mesh%half_cos_lat(j)
        ip_fv = ip_fv + tend%fv(i,j) * state%iap%u(i,j) * mesh%half_cos_lat(j)
        ip_cv = ip_cv + tend%cv(i,j) * state%iap%u(i,j) * mesh%half_cos_lat(j)
        ip_u_pgf = ip_u_pgf + tend%u_pgf(i,j) * state%iap%u(i,j) * mesh%half_cos_lat(j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_v_adv_lon = ip_v_adv_lon + tend%v_adv_lon(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_v_adv_lat = ip_v_adv_lat + tend%v_adv_lat(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_fu = ip_fu + tend%fu(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_cu = ip_cu + tend%cu(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
        ip_v_pgf = ip_v_pgf + tend%v_pgf(i,j) * state%iap%v(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_mass_div_lon = ip_mass_div_lon + tend%mass_div_lon(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%half_cos_lat(j)
        ip_mass_div_lat = ip_mass_div_lat + tend%mass_div_lat(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%half_cos_lat(j)
      end do
    end do

    print *, &
      ip_u_adv_lon, ip_v_adv_lon, ip_u_adv_lat, ip_v_adv_lat, &
      ip_u_pgf + ip_mass_div_lon, &
      ip_fv - ip_fu, &
      ip_cv - ip_cu, &
      ip_v_pgf + ip_mass_div_lat

  end subroutine check_antisymmetry
  
end module dycore_mod
