module space_operators_mod

  use data_mod
  use reduce_mod
  use types_mod
  use params_mod
  use mesh_mod
  use parallel_mod
  use filter_mod
  use pole_a_grid_mod

  implicit none

  private

  public space_operators_run

contains

  subroutine zonal_momentum_advection_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

    ! U
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1 .or. .not. reduce_adv_lon) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%u_adv_lon_a(i,j) = 0.5 / coef%full_dlon(j) * &
            ((state%u_a(i,j) + state%u_a(i+1,j)) * state%iap%u_a(i+1,j) - &
             (state%u_a(i,j) + state%u_a(i-1,j)) * state%iap%u_a(i-1,j))
        end do
      else
        tend%u_adv_lon_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.5 / coef%full_dlon(j) / full_reduce_factor(j) * &
              ((full_reduced_state(j)%u_a(i,k,2) + full_reduced_state(j)%u_a(i+1,k,2)) * full_reduced_state(j)%iap%u_a(i+1,k,2) - &
               (full_reduced_state(j)%u_a(i,k,2) + full_reduced_state(j)%u_a(i-1,k,2)) * full_reduced_state(j)%iap%u_a(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%u_adv_lon_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%u_adv_lon_a(:,j), left_halo=.true.)
      end if
    end do

    ! V
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1 .or. .not. reduce_adv_lon) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%v_adv_lon_a(i,j) = 0.5 / coef%full_dlon(j) * &
            ((state%u_a(i,j) + state%u_a(i+1,j)) * state%iap%v_a(i+1,j) - &
             (state%u_a(i,j) + state%u_a(i-1,j)) * state%iap%v_a(i-1,j))
        end do
      else
        tend%v_adv_lon_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.5 / coef%full_dlon(j) / full_reduce_factor(j) * &
              ((full_reduced_state(j)%u_a(i,k,2) + full_reduced_state(j)%u_a(i+1,k,2)) * full_reduced_state(j)%iap%v_a(i+1,k,2) - &
               (full_reduced_state(j)%u_a(i,k,2) + full_reduced_state(j)%u_a(i-1,k,2)) * full_reduced_state(j)%iap%v_a(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%v_adv_lon_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%v_adv_lon_a(:,j), left_halo=.true.)
      end if
    end do

  end subroutine zonal_momentum_advection_operator_on_a_grid

  subroutine meridional_momentum_advection_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    integer i, j

    ! U
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%u_adv_lat_a(i,j) = 0.5 / coef%full_dlat(j) * &
          ((state%v_a(i,j) * mesh%full_cos_lat(j) + state%v_a(i,j+1) * mesh%full_cos_lat(j+1)) * state%iap%u_a(i,j+1) - &
           (state%v_a(i,j) * mesh%full_cos_lat(j) + state%v_a(i,j-1) * mesh%full_cos_lat(j-1)) * state%iap%u_a(i,j-1))
                              
      end do
    end do

    ! V
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat_a(i,j) = 0.5 / coef%full_dlat(j) * &
          ((state%v_a(i,j) * mesh%full_cos_lat(j) + state%v_a(i,j+1) * mesh%full_cos_lat(j+1)) * state%iap%v_a(i,j+1) - &
           (state%v_a(i,j) * mesh%full_cos_lat(j) + state%v_a(i,j-1) * mesh%full_cos_lat(j-1)) * state%iap%v_a(i,j-1))
      end do
    end do

  end subroutine meridional_momentum_advection_operator_on_a_grid

  subroutine coriolis_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%fv_a(i,j) = coef%cori(j) * state%iap%v_a(i,j)
        end do
      else
        tend%fv_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = coef%cori(j) * full_reduced_state(j)%iap%v_a(i,k,2)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%fv_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fv_a(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%fu_a(i,j) = coef%cori(j) * state%iap%u_a(i,j)
        end do
      else
        tend%fu_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = coef%cori(j) * full_reduced_state(j)%iap%u_a(i,k,2)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%fu_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fu_a(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine coriolis_operator_on_a_grid

  subroutine curvature_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%cv_a(i,j) = coef%curv(j) * state%u_a(i,j) * state%iap%v_a(i,j)
        end do
      else
        tend%cv_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = coef%curv(j) * full_reduced_state(j)%u_a(i,k,2) * full_reduced_state(j)%iap%v_a(i,k,2)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%cv_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cv_a(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%cu_a(i,j) = coef%curv(j) * state%u_a(i,j) * state%iap%u_a(i,j)
        end do
      else
        tend%cu_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = coef%curv(j) * full_reduced_state(j)%u_a(i,k,2) * full_reduced_state(j)%iap%u_a(i,k,2)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%cu_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cu_a(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine curvature_operator_on_a_grid

  subroutine zonal_pressure_gradient_force_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%u_pgf_a(i,j) = state%iap%gd(i,j) / coef%full_dlon(j) * &
            (state%gd(i+1,j) + static%ghs(i+1,j) - state%gd(i-1,j) - static%ghs(i-1,j))
        end do
      else
        tend%u_pgf_a(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = full_reduced_state(j)%iap%gd(i,k,2) * &
              ((full_reduced_state(j)%gd(i+1,k,2) + full_reduced_static(j)%ghs(i+1,k,2)) - &
               (full_reduced_state(j)%gd(i-1,k,2) + full_reduced_static(j)%ghs(i-1,k,2))) / &
              coef%full_dlon(j) / full_reduce_factor(j)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%u_pgf_a(:,j))
        end do
        call parallel_overlay_inner_halo(tend%u_pgf_a(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine zonal_pressure_gradient_force_operator_on_a_grid

  subroutine meridional_pressure_gradient_force_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    integer i, j

    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf_a(i,j) = state%iap%gd(i,j) / coef%full_dlat(j) * mesh%full_cos_lat(j) * &
          (state%gd(i,j+1) + static%ghs(i,j+1) - state%gd(i,j-1) - static%ghs(i,j-1))
      end do
    end do

  end subroutine meridional_pressure_gradient_force_operator_on_a_grid

  subroutine zonal_mass_divergence_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%mass_div_lon(i,j) = (state%iap%gd(i+1,j) * state%iap%u_a(i+1,j) - &
                                    state%iap%gd(i-1,j) * state%iap%u_a(i-1,j)) &
                                   / coef%full_dlon(j)
        end do
      else
        tend%mass_div_lon(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = (full_reduced_state(j)%iap%gd(i+1,k,2) * full_reduced_state(j)%iap%u_a(i+1,k,2) - &
                               full_reduced_state(j)%iap%gd(i-1,k,2) * full_reduced_state(j)%iap%u_a(i-1,k,2)) / &
                              coef%full_dlon(j) / full_reduce_factor(j)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%mass_div_lon(:,j))
        end do
        call parallel_overlay_inner_halo(tend%mass_div_lon(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine zonal_mass_divergence_operator_on_a_grid

  subroutine meridional_mass_divergence_operator_on_a_grid(state, tend, pole)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pole

    real sp, np
    integer i, j

    do j = full_lat_start_idx_no_pole_a_grid(pole), full_lat_end_idx_no_pole_a_grid(pole)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = (state%iap%gd(i,j+1) * state%iap%v_a(i,j+1) * mesh%full_cos_lat(j+1) - &
                                  state%iap%gd(i,j-1) * state%iap%v_a(i,j-1) * mesh%full_cos_lat(j-1)) / &
                                 coef%full_dlat(j)
      end do
    end do

    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + state%iap%gd(i,j+1) * state%iap%v_a(i,j+1) * mesh%full_cos_lat(j+1)
      end do
      call parallel_zonal_sum(sp)
      sp = sp / mesh%num_full_lon / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sp
      end do
    end if

    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - state%iap%gd(i,j-1) * state%iap%v_a(i,j-1) * mesh%full_cos_lat(j-1)
      end do
      call parallel_zonal_sum(np)
      np = np / mesh%num_full_lon / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = np
      end do
    end if

  end subroutine meridional_mass_divergence_operator_on_a_grid

  subroutine zonal_momentum_advection_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    integer i, j, k

    ! U
    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      if (full_reduce_factor(j) == 1 .or. .not. reduce_adv_lon) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%u_adv_lon_c(i,j) = 0.5 / coef%full_dlon(j) * &
            ((state%u_c(i,j) + state%u_c(i+1,j)) * state%iap%u_c(i+1,j) - &
             (state%u_c(i,j) + state%u_c(i-1,j)) * state%iap%u_c(i-1,j))
        end do
      else
        tend%u_adv_lon_c(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.5 / coef%full_dlon(j) / full_reduce_factor(j) * &
              ((full_reduced_state(j)%u_c(i,k,2) + full_reduced_state(j)%u_c(i+1,k,2)) * full_reduced_state(j)%iap%u_c(i+1,k,2) - &
               (full_reduced_state(j)%u_c(i,k,2) + full_reduced_state(j)%u_c(i-1,k,2)) * full_reduced_state(j)%iap%u_c(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%u_adv_lon_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%u_adv_lon_c(:,j), left_halo=.true.)
      end if
    end do

    ! V
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      if (half_reduce_factor(j) == 1 .or. .not. reduce_adv_lon) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%v_adv_lon_c(i,j) = 0.5 / coef%half_dlon(j) * &
            ((state%u_c(i,  j) + state%u_c(i,  j+1)) * state%iap%v_c(i+1,j) - &
             (state%u_c(i-1,j) + state%u_c(i-1,j+1)) * state%iap%v_c(i-1,j))
        end do
      else
        tend%v_adv_lon_c(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = 0.5 / coef%half_dlon(j) / half_reduce_factor(j) * &
              ((half_reduced_state(j)%u_c(i,  k,2) + half_reduced_state(j)%u_c(i,  k,3)) * half_reduced_state(j)%iap%v_c(i+1,k,2) - &
               (half_reduced_state(j)%u_c(i-1,k,2) + half_reduced_state(j)%u_c(i-1,k,3)) * half_reduced_state(j)%iap%v_c(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, tend%v_adv_lon_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%v_adv_lon_c(:,j), left_halo=.true.)
      end if
    end do

  end subroutine zonal_momentum_advection_operator_on_c_grid

  subroutine meridional_momentum_advection_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U
    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_adv_lat_c(i,j) = 0.5 / coef%full_dlat(j) * &
          ((state%v_c(i,j  ) + state%v_c(i+1,j  )) * mesh%half_cos_lat(j  ) * state%iap%u_c(i,j+1) - &
           (state%v_c(i,j-1) + state%v_c(i+1,j-1)) * mesh%half_cos_lat(j-1) * state%iap%u_c(i,j-1))
      end do
    end do

    ! V
    do j = half_lat_start_idx_no_pole_c_grid, half_lat_end_idx_no_pole_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat_c(i,j) = 0.5 / coef%half_dlat(j) * &
          ((state%v_c(i,j) * mesh%half_cos_lat(j) + state%v_c(i,j+1) * mesh%half_cos_lat(j+1)) * state%iap%v_c(i,j+1) - &
           (state%v_c(i,j) * mesh%half_cos_lat(j) + state%v_c(i,j-1) * mesh%half_cos_lat(j-1)) * state%iap%v_c(i,j-1))
      end do
    end do

    ! Handle meridional advection at South Pole.
    if (parallel%has_south_pole) then
      j = parallel%half_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat_c(i,j) = 0.5 / coef%half_dlat(j) * &
          (state%v_c(i,j) * mesh%half_cos_lat(j) + state%v_c(i,j+1) * mesh%half_cos_lat(j+1)) * state%iap%v_c(i,j+1)
      end do
    end if

    ! Handle meridional advection at North Pole.
    if (parallel%has_north_pole) then
      j = parallel%half_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat_c(i,j) = - 0.5 / coef%half_dlat(j) * &
          (state%v_c(i,j) * mesh%half_cos_lat(j) + state%v_c(i,j-1) * mesh%half_cos_lat(j-1)) * state%iap%v_c(i,j-1)
      end do
    end if

  end subroutine meridional_momentum_advection_operator_on_c_grid

  subroutine coriolis_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    real c1, c2
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c1, c2, i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      c1 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%fv_c(i,j) = 0.25 * coef%cori(j) * &
            (c1 * (state%iap%v_c(i,j-1) + state%iap%v_c(i+1,j-1)) + &
             c2 * (state%iap%v_c(i,j  ) + state%iap%v_c(i+1,j  )))
        end do
      else
        tend%fv_c(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.25 * coef%cori(j) * &
              (c1 * (full_reduced_state(j)%iap%v_c(i,k,1) + full_reduced_state(j)%iap%v_c(i+1,k,1)) + &
               c2 * (full_reduced_state(j)%iap%v_c(i,k,2) + full_reduced_state(j)%iap%v_c(i+1,k,2)))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%fv_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fv_c(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      tend%fu_c(:,j) = 0.0
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%fu_c(i,j) = 0.25 * coef%cori(j) * (state%iap%u_c(i,j) + state%iap%u_c(i-1,j))
        end do
      else
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.25 * coef%cori(j) * (full_reduced_state(j)%iap%u_c(i,k,2) + full_reduced_state(j)%iap%u_c(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%fu_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fu_c(:,j), left_halo=.true.)
      end if
      if (full_reduce_factor(j+1) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%fu_c(i,j) = tend%fu_c(i,j) + 0.25 * coef%cori(j+1) * (state%iap%u_c(i,j+1) + state%iap%u_c(i-1,j+1))
        end do
      else
        ! Clear out right halo of tendency, because we will overlay them with left inner halo below.
        tend%fu_c(parallel%full_lon_end_idx+1:,j) = 0.0
        do k = 1, full_reduce_factor(j+1)
          do i = reduced_start_idx_at_full_lat(j+1), reduced_end_idx_at_full_lat(j+1)
            reduced_tend(i) = 0.25 * coef%cori(j+1) * (full_reduced_state(j+1)%iap%u_c(i,k,2) + full_reduced_state(j+1)%iap%u_c(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j+1, k, reduced_tend, tend%fu_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%fu_c(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine coriolis_operator_on_c_grid

  subroutine curvature_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    real c1, c2
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c1, c2, i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      c1 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      if (full_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%cv_c(i,j) = 0.25 * coef%curv(j) * state%u_c(i,j) * &
            (c1 * (state%iap%v_c(i,j-1) + state%iap%v_c(i+1,j-1)) + &
             c2 * (state%iap%v_c(i,j  ) + state%iap%v_c(i+1,j  )))
        end do
      else
        tend%cv_c(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.25 * coef%curv(j) * full_reduced_state(j)%u_c(i,k,2) * &
              (c1 * (full_reduced_state(j)%iap%v_c(i,k,1) + full_reduced_state(j)%iap%v_c(i+1,k,1)) + &
               c2 * (full_reduced_state(j)%iap%v_c(i,k,2) + full_reduced_state(j)%iap%v_c(i+1,k,2)))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%cv_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cv_c(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      tend%cu_c(:,j) = 0.0
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%cu_c(i,j) = 0.25 * coef%curv(j) * (state%u_c(i,  j) * state%iap%u_c(i,  j) + &
                                                state%u_c(i-1,j) * state%iap%u_c(i-1,j))
        end do
      else
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = 0.25 * coef%curv(j) * &
              (full_reduced_state(j)%u_c(i,  k,2) * full_reduced_state(j)%iap%u_c(i,  k,2) + &
               full_reduced_state(j)%u_c(i-1,k,2) * full_reduced_state(j)%iap%u_c(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%cu_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cu_c(:,j), left_halo=.true.)
      end if
      if (full_reduce_factor(j+1) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%cu_c(i,j) = tend%cu_c(i,j) + 0.25 * coef%curv(j+1) * (state%u_c(i,  j+1) * state%iap%u_c(i,  j+1) + &
                                                                 state%u_c(i-1,j+1) * state%iap%u_c(i-1,j+1))
        end do
      else
        ! Clear out right halo of tendency, because we will overlay them with left inner halo below.
        tend%cu_c(parallel%full_lon_end_idx+1:,j) = 0.0
        do k = 1, full_reduce_factor(j+1)
          do i = reduced_start_idx_at_full_lat(j+1), reduced_end_idx_at_full_lat(j+1)
            reduced_tend(i) = 0.25 * coef%curv(j+1) * &
              (full_reduced_state(j+1)%u_c(i,  k,2) * full_reduced_state(j+1)%iap%u_c(i,  k,2) + &
               full_reduced_state(j+1)%u_c(i-1,k,2) * full_reduced_state(j+1)%iap%u_c(i-1,k,2))
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j+1, k, reduced_tend, tend%cu_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%cu_c(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine curvature_operator_on_c_grid

  subroutine zonal_pressure_gradient_force_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%half_lon_start_idx:parallel%half_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      if (full_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%u_pgf_c(i,j) = (state%iap%gd(i,j) + state%iap%gd(i+1,j)) / coef%full_dlon(j) * &
            (state%gd(i+1,j) + static%ghs(i+1,j) - state%gd(i,j) - static%ghs(i,j))
        end do
      else
        tend%u_pgf_c(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = (full_reduced_state(j)%iap%gd(i,k,2) + full_reduced_state(j)%iap%gd(i+1,k,2)) * &
              ((full_reduced_state(j)%gd(i+1,k,2) + full_reduced_static(j)%ghs(i+1,k,2)) - &
               (full_reduced_state(j)%gd(i,  k,2) + full_reduced_static(j)%ghs(i,  k,2))) / &
              coef%full_dlon(j) / full_reduce_factor(j)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%u_pgf_c(:,j))
        end do
        call parallel_overlay_inner_halo(tend%u_pgf_c(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine zonal_pressure_gradient_force_operator_on_c_grid

  subroutine meridional_pressure_gradient_force_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf_c(i,j) = (state%iap%gd(i,j) + state%iap%gd(i,j+1)) / coef%half_dlat(j) * mesh%half_cos_lat(j) * &
          (state%gd(i,j+1) + static%ghs(i,j+1) - state%gd(i,j) - static%ghs(i,j))
      end do
    end do

  end subroutine meridional_pressure_gradient_force_operator_on_c_grid

  subroutine zonal_mass_divergence_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    integer i, j, k

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, reduced_tend)
    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%mass_div_lon(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%iap%u_c(i,  j) - &
                                    (state%iap%gd(i,j) + state%iap%gd(i-1,j)) * state%iap%u_c(i-1,j)) &
                                   / coef%full_dlon(j)
        end do
      else
        tend%mass_div_lon(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = ((full_reduced_state(j)%iap%gd(i,k,2) + full_reduced_state(j)%iap%gd(i+1,k,2)) * full_reduced_state(j)%iap%u_c(i,  k,2) - &
                               (full_reduced_state(j)%iap%gd(i,k,2) + full_reduced_state(j)%iap%gd(i-1,k,2)) * full_reduced_state(j)%iap%u_c(i-1,k,2)) / &
                              coef%full_dlon(j) / full_reduce_factor(j)
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, tend%mass_div_lon(:,j))
        end do
        call parallel_overlay_inner_halo(tend%mass_div_lon(:,j), left_halo=.true.)
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine zonal_mass_divergence_operator_on_c_grid

  subroutine meridional_mass_divergence_operator_on_c_grid(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real sp, np
    integer i, j

    do j = full_lat_start_idx_no_pole_c_grid, full_lat_end_idx_no_pole_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = ((state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v_c(i,j  ) * mesh%half_cos_lat(j  ) - &
                                  (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v_c(i,j-1) * mesh%half_cos_lat(j-1)) / &
                                 coef%full_dlat(j)
      end do
    end do

    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
      end do
      call parallel_zonal_sum(sp)
      sp = sp / mesh%num_full_lon / coef%full_dlat(j)
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
      np = np / mesh%num_full_lon / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = np
      end do
    end if

  end subroutine meridional_mass_divergence_operator_on_c_grid

  subroutine space_operators_run_on_a_grids(state, tend, dt, pass, pole)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real, intent(in) :: dt
    integer, intent(in) :: pass
    integer, intent(in) :: pole

    integer i, j, i1, i2
    real s1, s2

    ! call reduce_run(state, static)

    ! Allow me to use i1 and i2 as shorthands.
    i1 = parallel%full_lon_start_idx
    i2 = parallel%full_lon_end_idx

    select case (pass)
    case (all_pass)
      call zonal_momentum_advection_operator_on_a_grid(state, tend, pole)
      call meridional_momentum_advection_operator_on_a_grid(state, tend, pole)
      call coriolis_operator_on_a_grid(state, tend, pole)
      call curvature_operator_on_a_grid(state, tend, pole)
      call zonal_pressure_gradient_force_operator_on_a_grid(state, tend, pole)
      call meridional_pressure_gradient_force_operator_on_a_grid(state, tend, pole)
      call zonal_mass_divergence_operator_on_a_grid(state, tend, pole)
      call meridional_mass_divergence_operator_on_a_grid(state, tend, pole)

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du_a(i,j) = - tend%u_adv_lon_a(i,j) - tend%u_adv_lat_a(i,j) + tend%fv_a(i,j) + tend%cv_a(i,j) - tend%u_pgf_a(i,j)
        end do
      end do

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv_a(i,j) = - tend%v_adv_lon_a(i,j) - tend%v_adv_lat_a(i,j) - tend%fu_a(i,j) - tend%cu_a(i,j) - tend%v_pgf_a(i,j)
        end do
      end do

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
    case (slow_pass)
      call zonal_momentum_advection_operator_on_a_grid(state, tend, pole)
      call meridional_momentum_advection_operator_on_a_grid(state, tend, pole)

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du_a(i,j) = - tend%u_adv_lon_a(i,j) - tend%u_adv_lat_a(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du_a(i1:i2,j) * state%iap%u_a(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du_a(:,j))
            s2 = sum(tend%du_a(i1:i2,j) * state%iap%u_a(i1:i2,j))
            tend%du_a(i1:i2,j) = tend%du_a(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv_a(i,j) = - tend%v_adv_lon_a(i,j) - tend%v_adv_lat_a(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%dv_a(i1:i2,j) * state%iap%v_a(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dv_a(:,j))
            s2 = sum(tend%dv_a(i1:i2,j) * state%iap%v_a(i1:i2,j))
            tend%dv_a(i1:i2,j) = tend%dv_a(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      tend%dgd = 0.0
    case (fast_pass)
      call coriolis_operator_on_a_grid(state, tend, pole)
      call curvature_operator_on_a_grid(state, tend, pole)
      call zonal_pressure_gradient_force_operator_on_a_grid(state, tend, pole)
      call meridional_pressure_gradient_force_operator_on_a_grid(state, tend, pole)
      call zonal_mass_divergence_operator_on_a_grid(state, tend, pole)
      call meridional_mass_divergence_operator_on_a_grid(state, tend, pole)

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%mass_div_lon(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dgd(:,j))
            s2 = sum(tend%dgd(i1:i2,j) * (state%gd(i1:i2,j) + static%ghs(i1:i2,j)))
            tend%dgd(i1:i2,j) = tend%dgd(i1:i2,j) * s1 / s2
          end if
        end if
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%du_a(i,j) = tend%fv_a(i,j) + tend%cv_a(i,j) - tend%u_pgf_a(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du_a(i1:i2,j) * state%iap%u_a(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du_a(:,j))
            s2 = sum(tend%du_a(i1:i2,j) * state%iap%u_a(i1:i2,j))
            tend%du_a(i1:i2,j) = tend%du_a(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = full_lat_start_idx_a_grid(pole), full_lat_end_idx_a_grid(pole)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv_a(i,j) = - tend%fu_a(i,j) - tend%cu_a(i,j) - tend%v_pgf_a(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%dv_a(i1:i2,j) * state%iap%v_a(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%dv_a(:,j))
            s2 = sum(tend%dv_a(i1:i2,j) * state%iap%v_a(i1:i2,j))
            tend%dv_a(i1:i2,j) = tend%dv_a(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
    end select

  end subroutine space_operators_run_on_a_grids

  subroutine space_operators_run_on_c_grids(state, tend, dt, pass)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real, intent(in) :: dt
    integer, intent(in) :: pass

    integer i, j, i1, i2
    real s1, s2

    ! call reduce_run(state, static)

    ! Allow me to use i1 and i2 as shorthands.
    i1 = parallel%full_lon_start_idx
    i2 = parallel%full_lon_end_idx

    select case (pass)
    case (all_pass)
      call zonal_momentum_advection_operator_on_c_grid(state, tend)
      call meridional_momentum_advection_operator_on_c_grid(state, tend)
      call coriolis_operator_on_c_grid(state, tend)
      call curvature_operator_on_c_grid(state, tend)
      call zonal_pressure_gradient_force_operator_on_c_grid(state, tend)
      call meridional_pressure_gradient_force_operator_on_c_grid(state, tend)
      call zonal_mass_divergence_operator_on_c_grid(state, tend)
      call meridional_mass_divergence_operator_on_c_grid(state, tend)

      do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du_c(i,j) = - tend%u_adv_lon_c(i,j) - tend%u_adv_lat_c(i,j) + tend%fv_c(i,j) + tend%cv_c(i,j) - tend%u_pgf_c(i,j)
        end do
      end do

      do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv_c(i,j) = - tend%v_adv_lon_c(i,j) - tend%v_adv_lat_c(i,j) - tend%fu_c(i,j) - tend%cu_c(i,j) - tend%v_pgf_c(i,j)
        end do
      end do

      do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
    case (slow_pass)
      call zonal_momentum_advection_operator_on_c_grid(state, tend)
      call meridional_momentum_advection_operator_on_c_grid(state, tend)

      do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du_c(i,j) = - tend%u_adv_lon_c(i,j) - tend%u_adv_lat_c(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du_c(:,j))
            s2 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
            tend%du_c(i1:i2,j) = tend%du_c(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv_c(i,j) = - tend%v_adv_lon_c(i,j) - tend%v_adv_lat_c(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
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
      call coriolis_operator_on_c_grid(state, tend)
      call curvature_operator_on_c_grid(state, tend)
      call zonal_pressure_gradient_force_operator_on_c_grid(state, tend)
      call meridional_pressure_gradient_force_operator_on_c_grid(state, tend)
      call zonal_mass_divergence_operator_on_c_grid(state, tend)
      call meridional_mass_divergence_operator_on_c_grid(state, tend)

      do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
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

      do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du_c(i,j) = tend%fv_c(i,j) + tend%cv_c(i,j) - tend%u_pgf_c(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_full_zonal_tend(j)) then
          s1 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
          if (abs(s1) > 1.0e-16) then
            call filter_array_at_full_lat(j, tend%du_c(:,j))
            s2 = sum(tend%du_c(i1:i2,j) * state%iap%u_c(i1:i2,j))
            tend%du_c(i1:i2,j) = tend%du_c(i1:i2,j) * s1 / s2
          end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do

      do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv_c(i,j) = - tend%fu_c(i,j) - tend%cu_c(i,j) - tend%v_pgf_c(i,j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SMOOTHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (filter_half_zonal_tend(j)) then
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
    !  call history_write(tend, tag)
    ! end if

    ! call check_antisymmetry(tend, state)

  end subroutine space_operators_run_on_c_grids

  subroutine space_operators_run(state, tend, dt, pass)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real, intent(in) :: dt
    integer, intent(in) :: pass

    call fill_ac_interface(state)
    call space_operators_run_on_a_grids(state, tend, dt, pass, south_pole)
    call space_operators_run_on_a_grids(state, tend, dt, pass, north_pole)
    call space_operators_run_on_c_grids(state, tend, dt, pass)

  end subroutine space_operators_run

  subroutine fill_ac_interface(state)

    type(state_type), intent(inout) :: state

  end subroutine fill_ac_interface

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

    do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ip_u_adv_lon = ip_u_adv_lon + tend%u_adv_lon_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_u_adv_lat = ip_u_adv_lat + tend%u_adv_lat_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_fv = ip_fv + tend%fv_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_cv = ip_cv + tend%cv_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
        ip_u_pgf = ip_u_pgf + tend%u_pgf_c(i,j) * state%iap%u_c(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_v_adv_lon = ip_v_adv_lon + tend%v_adv_lon_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_v_adv_lat = ip_v_adv_lat + tend%v_adv_lat_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_fu = ip_fu + tend%fu_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_cu = ip_cu + tend%cu_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
        ip_v_pgf = ip_v_pgf + tend%v_pgf_c(i,j) * state%iap%v_c(i,j) * mesh%half_cos_lat(j)
      end do
    end do

    do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_mass_div_lon = ip_mass_div_lon + tend%mass_div_lon(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%full_cos_lat(j)
        ip_mass_div_lat = ip_mass_div_lat + tend%mass_div_lat(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%full_cos_lat(j)
      end do
    end do

    print *, &
      ip_u_adv_lon + ip_v_adv_lon + ip_u_adv_lat + ip_v_adv_lat, &
      ip_u_pgf + ip_mass_div_lon, &
      ip_fv - ip_fu, &
      ip_cv - ip_cu, &
      ip_v_pgf + ip_mass_div_lat

  end subroutine check_antisymmetry

end module space_operators_mod