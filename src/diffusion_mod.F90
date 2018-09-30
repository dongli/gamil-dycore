module diffusion_mod

  use parallel_mod
  use params_mod
  use mesh_mod
  use types_mod
  use data_mod
  use diag_mod
  use filter_mod

  implicit none

  private

  public diffusion_init
  public divergence_diffusion
  public normal_diffusion
  public diffusion_final

  real, allocatable :: ud(:,:)
  real, allocatable :: vd(:,:)
  real, allocatable :: gdd(:,:)

contains

  subroutine diffusion_init()

    if (.not. allocated(ud))  call parallel_allocate(ud, half_lon=.true.)
    if (.not. allocated(vd))  call parallel_allocate(vd, half_lat=.true.)
    if (.not. allocated(gdd)) call parallel_allocate(gdd)

  end subroutine diffusion_init

  subroutine divergence_diffusion(dt, diag, state)

    real, intent(in) :: dt
    type(diag_type), intent(in) :: diag
    type(state_type), intent(inout) :: state

    real, parameter :: vd = 1.0e6
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      ! if (filter_full_zonal_tend(j)) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          state%u(i,j) = state%u(i,j) + vd * (diag%div(i+1,j) - diag%div(i,j)) / coef%full_dlon(j)
        end do
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          state%iap%u(i,j) = state%u(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j))
        end do
      ! end if
    end do
    do j = parallel%half_lon_start_idx, parallel%half_lat_end_idx
      ! if (filter_half_zonal_tend(j)) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          state%v(i,j) = state%v(i,j) + vd * (diag%div(i,j+1) - diag%div(i,j)) / radius / mesh%dlat
        end do
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          state%iap%v(i,j) = state%v(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1))
        end do
      ! end if
    end do

    call parallel_fill_halo(state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(state%v(:,:), all_halo=.true.)

  end subroutine divergence_diffusion

  subroutine normal_diffusion(dt, state)

    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    real sp, np
    integer i, j, k

    !
    ! ∂ F         1      ∂² F        1      ∂         ∂ F
    ! --- = 𝞶 ---------- ---- + 𝞶 --------- --- cos(φ) ---
    ! ∂ t     a² cos²(φ) ∂ λ²     a² cos(φ) ∂ φ        ∂ φ
    !

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ud(i,j) = (state%u(i+1,j) - 2 * state%u(i,j) + state%u(i-1,j)) / coef%full_dlon(j)**2 + &
                   ((state%u(i,j+1) - state%u(i,j)) * mesh%half_cos_lat(j) - &
                    (state%u(i,j) - state%u(i,j-1)) * mesh%half_cos_lat(j-1)) / &
                    coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
      end do
    end do

    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_array_at_full_lat(j, ud(:,j))
      end if
    end do

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + dt * diffusion_coef * ud(i,j)
      end do
    end do

    call parallel_fill_halo(state%u, all_halo=.true.)

    ! V
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = (state%v(i+1,j) - 2 * state%v(i,j) + state%v(i-1,j)) / coef%half_dlon(j)**2
        end do
    end do
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd(i,j) = vd(i,j) + ((state%v(i,j+1) - state%v(i,j)) * mesh%full_cos_lat(j+1) - &
                             (state%v(i,j) - state%v(i,j-1)) * mesh%full_cos_lat(j)) / &
                             coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
      end do
    end do
    if (parallel%has_south_pole) then
      j = parallel%half_lat_start_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd(i,j) = vd(i,j) + (state%v(i,j+1) - state%v(i,j)) * mesh%full_cos_lat(j+1) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
      end do
    end if
    if (parallel%has_north_pole) then
      j = parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd(i,j) = vd(i,j) - (state%v(i,j) - state%v(i,j-1)) * mesh%full_cos_lat(j) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
      end do
    end if

    ! Do FFT filter.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (filter_half_zonal_tend(j)) then
        call filter_array_at_half_lat(j, vd(:,j))
      end if
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + dt * diffusion_coef * vd(i,j)
      end do
    end do

    call parallel_fill_halo(state%v, all_halo=.true.)

    ! H
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd(i,j) = (state%gd(i+1,j) - 2 * state%gd(i,j) + state%gd(i-1,j)) / coef%full_dlon(j)**2 + &
                    ((state%gd(i,j+1) - state%gd(i,j)) * mesh%half_cos_lat(j) - &
                     (state%gd(i,j) - state%gd(i,j-1)) * mesh%half_cos_lat(j-1)) / &
                     coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
      end do
    end do
    if (parallel%has_south_pole) then
      j = parallel%full_lat_start_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + state%gd(i,j+1) - state%gd(i,j)
      end do
      call parallel_zonal_sum(sp)
      sp = sp * mesh%half_cos_lat(j) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd(i,j) = sp
      end do
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_end_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - (state%gd(i,j) - state%gd(i,j-1))
      end do
      call parallel_zonal_sum(np)
      np = np * mesh%half_cos_lat(j-1) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd(i,j) = np
      end do
    end if

    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_array_at_full_lat(j, gdd(:,j))
      end if
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%gd(i,j) = state%gd(i,j) + dt * diffusion_coef * gdd(i,j)
      end do
    end do

    call parallel_fill_halo(state%gd, all_halo=.true.)

    call iap_transform(state)

  end subroutine normal_diffusion

  subroutine diffusion_final()

    if (allocated(ud))  deallocate(ud)
    if (allocated(vd))  deallocate(vd)
    if (allocated(gdd)) deallocate(gdd)
  end subroutine diffusion_final

end module diffusion_mod
