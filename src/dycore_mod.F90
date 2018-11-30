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
      integer, intent(in) :: old_time_idx
      integer, intent(in) :: new_time_idx
      integer, intent(in) :: pass
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

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du_c(i,j) = - tend%u_adv_lon_c(i,j) - tend%u_adv_lat_c(i,j) + tend%fv_c(i,j) - tend%u_pgf_c(i,j)
          tend%dv_d(i,j) = - tend%v_adv_lon_d(i,j) - tend%v_adv_lat_d(i,j) - tend%fu_d(i,j) - tend%v_pgf_d(i,j)
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du_d(i,j) = - tend%u_adv_lon_d(i,j) - tend%u_adv_lat_d(i,j) + tend%fv_d(i,j) - tend%u_pgf_d(i,j)
          tend%dv_c(i,j) = - tend%v_adv_lon_c(i,j) - tend%v_adv_lat_c(i,j) - tend%fu_c(i,j) - tend%v_pgf_c(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
    case (slow_pass)
#ifndef NDEBUG
      tend%fv_c = 0.0
      tend%fu_c = 0.0
      tend%u_pgf_c = 0.0
      tend%v_pgf_c = 0.0
      tend%fv_d = 0.0
      tend%fu_d = 0.0
      tend%u_pgf_d = 0.0
      tend%v_pgf_d = 0.0
      tend%mass_div_lon = 0.0
      tend%mass_div_lat = 0.0
#endif
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du_c(i,j) = - tend%u_adv_lon_c(i,j) - tend%u_adv_lat_c(i,j)
          tend%dv_d(i,j) = - tend%v_adv_lon_d(i,j) - tend%v_adv_lat_d(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du_c(:,j))
            s2 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
            tend%du_c(i1:i2,j) = tend%du_c(i1:i2,j) * s1 / s2
          end if
          s1 = sum(tend%dv_d(i1:i2,j) * state%iap%v_d(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dv_d(:,j))
            s2 = sum(tend%dv_d(i1:i2,j) * state%iap%v_d(i1:i2,j))
            tend%dv_d(i1:i2,j) = tend%dv_d(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du_d(i,j) = - tend%u_adv_lon_d(i,j) - tend%u_adv_lat_d(i,j)
          tend%dv_c(i,j) = - tend%v_adv_lon_c(i,j) - tend%v_adv_lat_c(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
          s1 = sum(tend%du_d(i1:i2,j) * state%iap%u_d(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_half_lat(j, tend%du_d(:,j))
            s2 = sum(tend%du_d(i1:i2,j) * state%iap%u_d(i1:i2,j))
            tend%du_d(i1:i2,j) = tend%du_d(i1:i2,j) * s1 / s2
          end if
          s1 = sum(tend%dv_c(i1:i2,j) * state%iap%v_c(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_half_lat(j, tend%dv_c(:,j))
            s2 = sum(tend%dv_c(i1:i2,j) * state%iap%v_c(i1:i2,j))
            tend%dv_c(i1:i2,j) = tend%dv_c(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      tend%dgd = 0.0
    case (fast_pass)
#ifndef NDEBUG
      tend%u_adv_lon_c = 0.0
      tend%u_adv_lat_c = 0.0
      tend%v_adv_lon_c = 0.0
      tend%v_adv_lat_c = 0.0
      tend%u_adv_lon_d = 0.0
      tend%u_adv_lat_d = 0.0
      tend%v_adv_lon_d = 0.0
      tend%v_adv_lat_d = 0.0
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

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du_c(i,j) =   tend%fv_c(i,j) - tend%u_pgf_c(i,j)
          tend%dv_d(i,j) = - tend%fu_d(i,j) - tend%v_pgf_d(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du_c(:,j))
            s2 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
            tend%du_c(i1:i2,j) = tend%du_c(i1:i2,j) * s1 / s2
          end if
          s1 = sum(tend%dv_d(i1:i2,j) * state%iap%v_d(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dv_d(:,j))
            s2 = sum(tend%dv_d(i1:i2,j) * state%iap%v_d(i1:i2,j))
            tend%dv_d(i1:i2,j) = tend%dv_d(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du_d(i,j) =   tend%fv_d(i,j) - tend%u_pgf_d(i,j)
          tend%dv_c(i,j) = - tend%fu_c(i,j) - tend%v_pgf_c(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
          s1 = sum(tend%du_d(i1:i2,j) * state%iap%u_d(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_half_lat(j, tend%du_d(:,j))
            s2 = sum(tend%du_d(i1:i2,j) * state%iap%u_d(i1:i2,j))
            tend%du_d(i1:i2,j) = tend%du_d(i1:i2,j) * s1 / s2
          end if
          s1 = sum(tend%dv_c(i1:i2,j) * state%iap%v_c(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_half_lat(j, tend%dv_c(:,j))
            s2 = sum(tend%dv_c(i1:i2,j) * state%iap%v_c(i1:i2,j))
            tend%dv_c(i1:i2,j) = tend%dv_c(i1:i2,j) * s1 / s2
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

    integer i, j

    ! U on C grid
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_adv_lon_c(i,j) = 0.25 / coef%full_dlon(j) * &
          ((state%u_c(i,j) + state%u_c(i+1,j)) * state%iap%u_c(i+1,j) - &
           (state%u_c(i,j) + state%u_c(i-1,j)) * state%iap%u_c(i-1,j))
      end do
    end do

    ! V on C grid
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lon_c(i,j) = 0.25 / coef%half_dlon(j) * &
          ((state%u_c(i,  j) + state%u_c(i,  j+1)) * state%iap%v_c(i+1,j) - &
           (state%u_c(i-1,j) + state%u_c(i-1,j+1)) * state%iap%v_c(i-1,j))
      end do
    end do

    ! U on D grid
    !
    !       o u_d (i-1,j)    o u_d (i,j)    o u_d (i+1,j)
    !
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_adv_lon_d(i,j) = 0.25 / coef%half_dlon(j) * &
          ((state%u_d(i,j) + state%u_d(i+1,j)) * state%iap%u_d(i+1,j) - &
           (state%u_d(i,j) + state%u_d(i-1,j)) * state%iap%u_d(i-1,j))
      end do
    end do

    ! V on D grid
    !
    !                      x u_d (i,j-1)                   x u_d (i+1,j-1)
    !
    !       o v_d (i-1,j)                   o v_d (i,j)                      o v_d (i+1,j)
    !
    !                      x u_d (i, j )                   x u_d (i+1, j )
    !
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%v_adv_lon_d(i,j) = 0.25 / coef%full_dlon(j) * &
          ((state%u_d(i+1,j-1) + state%u_d(i+1,j)) * state%iap%v_d(i+1,j) - &
           (state%u_d(i,  j-1) + state%u_d(i,  j)) * state%iap%v_d(i-1,j))
      end do
    end do

  end subroutine zonal_momentum_advection_operator

  subroutine meridional_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U on C grid
    !
    !                        o u_c (i,j-1)
    !
    !       x v_c (i,j-1)                    x v_c (i+1,j-1)
    !
    !                        o u_c (i, j )
    !
    !       x v_c (i,j  )                    x v_c (i+1,j  )
    !
    !                        o u_c (i,j+1)
    !
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_adv_lat_c(i,j) = 0.25 / coef%full_dlat(j) * &
          ((state%v_c(i,j  ) + state%v_c(i+1,j  )) * mesh%half_cos_lat(j  ) * state%iap%u_c(i,j+1) - &
           (state%v_c(i,j-1) + state%v_c(i+1,j-1)) * mesh%half_cos_lat(j-1) * state%iap%u_c(i,j-1))
      end do
    end do

    ! V on C grid
    !
    !                        o v_c (i,j-1)
    !
    !                        o v_c (i, j )
    !
    !                        o v_c (i,j+1)
    !
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat_c(i,j) = 0.25 / coef%half_dlat(j) * &
          ((state%v_c(i,j) * mesh%half_cos_lat(j) + state%v_c(i,j+1) * mesh%half_cos_lat(j+1)) * state%iap%v_c(i,j+1) - &
           (state%v_c(i,j) * mesh%half_cos_lat(j) + state%v_c(i,j-1) * mesh%half_cos_lat(j-1)) * state%iap%v_c(i,j-1))
      end do
    end do

    ! U on D grid
    !
    !                        o u_d (i,j-1)
    !
    !     x v_d (i-1,j  )                    x v_d (i,j  )
    !
    !                        o u_d (i, j )
    !
    !     x v_d (i-1,j+1)                    x v_d (i,j+1)
    !
    !                        o u_d (i,j+1)
    !
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_adv_lat_d(i,j) = 0.25 / coef%half_dlat(j) * &
          ((state%v_d(i-1,j+1) + state%v_d(i,j+1)) * mesh%full_cos_lat(j+1) * state%iap%u_d(i,j+1) - &
           (state%v_d(i-1,j  ) + state%v_d(i,j  )) * mesh%full_cos_lat(j  ) * state%iap%u_d(i,j-1))
      end do
    end do

    ! V on D grid
    !
    !                        o v_d (i,j-1)
    !
    !                        o v_d (i, j )
    !
    !                        o v_d (i,j+1)
    !
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%v_adv_lat_d(i,j) = 0.25 / coef%full_dlat(j) * &
          ((state%v_d(i,j) * mesh%full_cos_lat(j) + state%v_d(i,j+1) * mesh%full_cos_lat(j+1)) * state%iap%v_d(i,j+1) - &
           (state%v_d(i,j) * mesh%full_cos_lat(j) + state%v_d(i,j-1) * mesh%full_cos_lat(j-1)) * state%iap%v_d(i,j-1))
      end do
    end do

  end subroutine meridional_momentum_advection_operator

  subroutine coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U on C grid
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%fv_c(i,j) = (coef%full_f(j) + coef%full_c(j) * state%u_c(i,j)) * state%iap%v_d(i,j)
      end do
    end do

    ! V on C grid
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%fu_c(i,j) = (coef%half_f(j) + coef%half_c(j) * state%u_d(i,j)) * state%iap%u_d(i,j)
      end do
    end do

    ! U on D grid
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%fv_d(i,j) = (coef%half_f(j) + coef%half_c(j) * state%u_d(i,j)) * state%iap%v_c(i,j)
      end do
    end do

    ! V on D grid
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%fu_d(i,j) = (coef%full_f(j) + coef%full_c(j) * state%u_c(i,j)) * state%iap%u_c(i,j)
      end do
    end do

  end subroutine coriolis_operator

  subroutine zonal_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U on C grid
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_pgf_c(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) / coef%full_dlon(j) * &
          (state%gd(i+1,j) + static%ghs(i+1,j) - state%gd(i,j) - static%ghs(i,j))
      end do
    end do

    ! U on D grid
    !
    !       x gd  (i-1, j )    x gd  (i, j )    x gd (i+1, j )
    !
    !                          o u_d (i, j )
    !
    !       x gd  (i-1,j+1)    x gd  (i,j+1)    x gd (i+1,j+1)
    !
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_pgf_d(i,j) = 0.125 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) / coef%half_dlon(j) * &
          ((state%gd(i+1,j) + state%gd(i+1,j+1) + static%ghs(i+1,j) + static%ghs(i+1,j+1)) - &
           (state%gd(i-1,j) + state%gd(i-1,j+1) + static%ghs(i-1,j) + static%ghs(i-1,j+1)))
      end do
    end do

  end subroutine zonal_pressure_gradient_force_operator

  subroutine meridional_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! V on C grid
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf_c(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) / coef%half_dlat(j) * mesh%half_cos_lat(j) * &
          (state%gd(i,j+1) + static%ghs(i,j+1) - state%gd(i,j) - static%ghs(i,j))
      end do
    end do

    ! V on D grid
    !
    !       x gd (i,j-1)                            x gd (i+1,j-1)
    !
    !       x gd (i, j)          o v_d (i,j)        x gd (i+1, j )
    !
    !       x gd (i,j+1)                            x gd (i+1,j+1)
    !
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%v_pgf_d(i,j) = 0.125 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) / coef%full_dlat(j) * mesh%full_cos_lat(j) * &
          ((state%gd(i,j+1) + state%gd(i+1,j+1) + static%ghs(i,j+1) + static%ghs(i+1,j+1)) - &
           (state%gd(i,j-1) + state%gd(i+1,j-1) + static%ghs(i,j-1) + static%ghs(i+1,j-1)))
      end do
    end do

  end subroutine meridional_pressure_gradient_force_operator

  subroutine zonal_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lon(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%iap%u_c(i,  j) - &
                                  (state%iap%gd(i,j) + state%iap%gd(i-1,j)) * state%iap%u_c(i-1,j)) &
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
        tend%mass_div_lat(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v_c(i,j  ) * mesh%half_cos_lat(j  ) - &
                                  (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v_c(i,j-1) * mesh%half_cos_lat(j-1)) &
                                 * 0.5 / coef%full_dlat(j)
      end do
    end do

    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
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
        np = np - (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v_c(i,j-1) * mesh%half_cos_lat(j-1)
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

    integer i, j, i1, i2

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
        new_state%iap%u_c(i,j) = old_state%iap%u_c(i,j) + dt * tend%du_c(i,j)
        new_state%iap%v_d(i,j) = old_state%iap%v_d(i,j) + dt * tend%dv_d(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%u_d(i,j) = old_state%iap%u_d(i,j) + dt * tend%du_d(i,j)
        new_state%iap%v_c(i,j) = old_state%iap%v_c(i,j) + dt * tend%dv_c(i,j)
      end do
    end do

    ! Transform from IAP to normal state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u_c(i,j) = new_state%iap%u_c(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i+1,j))
        new_state%v_d(i,j) = new_state%iap%v_d(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%u_d(i,j) = new_state%iap%u_d(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i,j+1))
        new_state%v_c(i,j) = new_state%iap%v_c(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i,j+1))
      end do
    end do

    call parallel_fill_halo(new_state%iap%u_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%u_d(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v_d(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v_c(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u_d(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v_d(:,:), all_halo=.true.)

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend1%du_c(i,j) * tend2%du_c(i,j) * mesh%full_cos_lat(j)
        res = res + tend1%dv_d(i,j) * tend2%dv_d(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%du_d(i,j) * tend2%du_d(i,j) * mesh%half_cos_lat(j)
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
      call integrator(time_step_size, old_time_idx, new_time_idx, all_pass)
    end select

    if (use_diffusion) then
      call ordinary_diffusion(time_step_size, state(new_time_idx))
    end if

  end subroutine time_integrate

  subroutine predict_correct(time_step_size, old, new, pass)

    real, intent(in) :: time_step_size
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    real dt, ip1, ip2, beta

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
    beta = merge(ip1 / ip2, 1.0, qcon_modified .and. ip1 /= 0.0 .and. ip2 /= 0.0)
    call log_add_diag('beta', beta)

    dt = time_step_size * beta
    call update_state(dt, tend(new), state(old), state(new))

  end subroutine predict_correct

  subroutine check_antisymmetry(tend, state)

    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: state

    integer i, j
    real ip_u_adv_lon_c
    real ip_u_adv_lat_c
    real ip_v_adv_lon_c
    real ip_v_adv_lat_c
    real ip_fv_c
    real ip_fu_c
    real ip_u_pgf_c
    real ip_v_pgf_c
    real ip_u_adv_lon_d
    real ip_u_adv_lat_d
    real ip_v_adv_lon_d
    real ip_v_adv_lat_d
    real ip_fv_d
    real ip_fu_d
    real ip_u_pgf_d
    real ip_v_pgf_d
    real ip_mass_div_lon
    real ip_mass_div_lat

    ip_u_adv_lon_c = 0.0
    ip_u_adv_lat_c = 0.0
    ip_v_adv_lon_c = 0.0
    ip_v_adv_lat_c = 0.0
    ip_fu_c = 0.0
    ip_fv_c = 0.0
    ip_u_pgf_c = 0.0
    ip_v_pgf_c = 0.0
    ip_u_adv_lon_d = 0.0
    ip_u_adv_lat_d = 0.0
    ip_v_adv_lon_d = 0.0
    ip_v_adv_lat_d = 0.0
    ip_fu_d = 0.0
    ip_fv_d = 0.0
    ip_u_pgf_d = 0.0
    ip_v_pgf_d = 0.0
    ip_mass_div_lon = 0.0
    ip_mass_div_lat = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ip_u_adv_lon_c = ip_u_adv_lon_c + tend%u_adv_lon_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_u_adv_lat_c = ip_u_adv_lat_c + tend%u_adv_lat_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_v_adv_lon_d = ip_v_adv_lon_d + tend%v_adv_lon_d(i,j) * state%iap%v_d(i,j) * mesh%full_cos_lat(j)
        ip_v_adv_lat_d = ip_v_adv_lat_d + tend%v_adv_lat_d(i,j) * state%iap%v_d(i,j) * mesh%full_cos_lat(j)
        ip_fv_c = ip_fv_c + tend%fv_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_fu_d = ip_fu_d + tend%fu_d(i,j) * state%iap%v_d(i,j) * mesh%full_cos_lat(j)
        ip_u_pgf_c = ip_u_pgf_c + tend%u_pgf_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_v_pgf_d = ip_v_pgf_d + tend%v_pgf_d(i,j) * state%iap%v_d(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_u_adv_lon_d = ip_u_adv_lon_d + tend%u_adv_lon_d(i,j) * state%iap%u_d(i,j) * mesh%half_cos_lat(j)
        ip_u_adv_lat_d = ip_u_adv_lat_d + tend%u_adv_lat_d(i,j) * state%iap%u_d(i,j) * mesh%half_cos_lat(j)
        ip_v_adv_lon_c = ip_v_adv_lon_c + tend%v_adv_lon_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_v_adv_lat_c = ip_v_adv_lat_c + tend%v_adv_lat_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_fv_d = ip_fv_d + tend%fv_d(i,j) * state%iap%u_d(i,j) * mesh%half_cos_lat(j)
        ip_fu_c = ip_fu_c + tend%fu_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_u_pgf_d = ip_u_pgf_d + tend%u_pgf_d(i,j) * state%iap%u_d(i,j) * mesh%half_cos_lat(j)
        ip_v_pgf_c = ip_v_pgf_c + tend%v_pgf_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_mass_div_lon = ip_mass_div_lon + tend%mass_div_lon(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%full_cos_lat(j)
        ip_mass_div_lat = ip_mass_div_lat + tend%mass_div_lat(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%full_cos_lat(j)
      end do
    end do

    print *, &
      ip_u_adv_lon_c + ip_v_adv_lon_c + ip_u_adv_lat_c + ip_v_adv_lat_c, &
      ip_u_adv_lon_d + ip_v_adv_lon_d + ip_u_adv_lat_d + ip_v_adv_lat_d, &
      ip_u_pgf_c + ip_mass_div_lon, &
      ip_fv_c - ip_fu_d, &
      ip_fu_c - ip_fv_d, &
      ip_v_pgf_c + ip_mass_div_lat

  end subroutine check_antisymmetry
  
end module dycore_mod
